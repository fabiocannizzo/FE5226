#include "MCEngine.h"
#include "Models.h"
#include <algorithm>
#include <random>

namespace mc {

MCEngine::MCEngine(std::vector<pmcpricer_t> pricers,
                   Date                     today,
                   size_t                   n_paths,
                   uint64_t                 seed,
                   bool                     antithetic)
    : m_pricers(std::move(pricers))
    , m_today(today)
    , m_n_paths(n_paths)
    , m_seed(seed)
    , m_antithetic(antithetic) {
    MYASSERT(!m_pricers.empty(), "MCEngine requires at least one pricer");
    MYASSERT(m_n_paths > 0, "Number of Monte Carlo paths must be positive");
    buildMasterTimelineAndModels();
}

void MCEngine::buildMasterTimelineAndModels() {
    m_dates.clear();
    m_models.clear();
    m_dates.push_back(m_today); // t_0 is always today

    for (const auto& pricer : m_pricers) {
        auto dts = pricer->requiredDates();
        m_dates.insert(m_dates.end(), dts.begin(), dts.end());
        auto mdls = pricer->modelNeeded();
        m_models.insert(m_models.end(), mdls.begin(), mdls.end());
    }

    // Sort and deduplicate master timeline
    std::sort(m_dates.begin(), m_dates.end());
    m_dates.erase(std::unique(m_dates.begin(), m_dates.end()), m_dates.end());
    MYASSERT(m_dates.front() == m_today, "Payoff required date cannot precede today");

    // Sort and deduplicate required models
    std::sort(m_models.begin(), m_models.end());
    m_models.erase(std::unique(m_models.begin(), m_models.end()), m_models.end());
}

std::vector<double> MCEngine::price(const MarketDataServer& mkt) const {
    const size_t n_models = m_models.size();

    // 1. Query Model Factory with collated models needed and prepopulated MarketDataServer
    std::vector<std::shared_ptr<IModelMC>> models = createModelsMC(m_models, mkt);

    // 2. Add intermediate dates where required by models (e.g. Euler <= 30d monthly steps)
    for (size_t m = 0; m < n_models; ++m) {
        m_dates = models[m]->mergeDates(m_dates);
    }
    std::sort(m_dates.begin(), m_dates.end());
    m_dates.erase(std::unique(m_dates.begin(), m_dates.end()), m_dates.end());
    const size_t n_dates = m_dates.size();
    const size_t n_steps = n_dates - 1;

    // 3. Configure model timelines and bind indices on pricers
    std::vector<size_t>      brownian_offset(n_models + 1, 0);
    std::vector<std::string> factor_names;

    for (size_t m = 0; m < n_models; ++m) {
        models[m]->configureTimeline(m_dates);
        size_t dW              = models[m]->numberOfBrownians();
        brownian_offset[m + 1] = brownian_offset[m] + dW;
        for (size_t b = 0; b < dW; ++b) {
            factor_names.push_back(m_models[m] + "." + std::to_string(b));
        }
    }
    const size_t total_brownians = brownian_offset.back();

    std::vector<const IModelMC*> model_ptrs;
    for (const auto& m_ptr : models) {
        model_ptrs.push_back(m_ptr.get());
    }
    for (auto& pricer : m_pricers) {
        pricer->bindIndices(m_dates, m_models, model_ptrs);
    }

    // Pre-fetch settlement discount factors for each pricer
    std::vector<double> df_settle(m_pricers.size(), 1.0);
    for (size_t pr = 0; pr < m_pricers.size(); ++pr) {
        std::string ccy = m_pricers[pr]->baseCurrency();
        auto        it  = mkt.discount_curves.find(ccy);
        if (it != mkt.discount_curves.end() && it->second != nullptr) {
            df_settle[pr] = it->second->df(m_pricers[pr]->settlementDate());
        }
    }

    // 4. Assemble global M x M correlation matrix Sigma and compute Cholesky factor L
    Matrix Sigma(total_brownians, total_brownians, 1.0);
    for (size_t i = 0; i < total_brownians; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double rho  = mkt.getCorrelation(factor_names[i], factor_names[j]);
            Sigma(i, j) = Sigma(j, i) = rho;
        }
    }
    Matrix L = cholesky(Sigma);

    // 5. Pre-allocate reusable path workspace outside the Monte Carlo loop
    Matrix              rndNumbers(total_brownians, n_steps);
    Matrix              Z_iid(total_brownians, n_steps);
    std::vector<Matrix> stateVariables(n_models);
    for (size_t m = 0; m < n_models; ++m) {
        stateVariables[m] = Matrix(models[m]->numberOfStateVariables(), n_dates);
    }

    std::mt19937_64                  rng(m_seed);
    std::normal_distribution<double> norm(0.0, 1.0);
    std::vector<double>              pv_sum(m_pricers.size(), 0.0);

    // Helper lambda to diffuse all models given correlated rndNumbers and evaluate payoffs
    auto diffuseAndPrice = [&]() {
        for (size_t m = 0; m < n_models; ++m) {
            ConstMatrixView rnd_slice = rndNumbers.rowSlice(brownian_offset[m], brownian_offset[m + 1]);
            MatrixView      state_view(stateVariables[m]);
            models[m]->diffuse(rnd_slice, state_view);
        }
        for (size_t pr = 0; pr < m_pricers.size(); ++pr) {
            pv_sum[pr] += m_pricers[pr]->pricePath(stateVariables, df_settle[pr], /*fx*/ 1.0);
        }
    };

    // 6. Main Monte Carlo Path Loop
    size_t p = 0;
    while (p < m_n_paths) {
        for (size_t i = 0; i < total_brownians; ++i) {
            for (size_t k = 0; k < n_steps; ++k) {
                Z_iid(i, k) = norm(rng);
            }
        }

        // Correlate random numbers at each step: rndNumbers(:, k) = L * Z_iid(:, k)
        for (size_t k = 0; k < n_steps; ++k) {
            for (size_t i = 0; i < total_brownians; ++i) {
                double acc = 0.0;
                for (size_t j = 0; j <= i; ++j) {
                    acc += L(i, j) * Z_iid(j, k);
                }
                rndNumbers(i, k) = acc;
            }
        }

        diffuseAndPrice();
        ++p;

        // Antithetic variate path (-rndNumbers)
        if (m_antithetic && p < m_n_paths) {
            for (size_t i = 0; i < total_brownians; ++i) {
                for (size_t k = 0; k < n_steps; ++k) {
                    rndNumbers(i, k) = -rndNumbers(i, k);
                }
            }
            diffuseAndPrice();
            ++p;
        }
    }

    // 7. Average across all simulated paths
    std::vector<double> pv_results(m_pricers.size(), 0.0);
    for (size_t pr = 0; pr < m_pricers.size(); ++pr) {
        pv_results[pr] = pv_sum[pr] / static_cast<double>(m_n_paths);
    }
    return pv_results;
}

std::vector<double> MCEngine::price(const std::map<std::string, IMarketObject*>& mkt) const {
    const size_t n_models = m_models.size();

    // 1. Retrieve active IModelMC* objects, let models inject intermediate dates
    std::vector<IModelMC*> models(n_models, nullptr);
    for (size_t m = 0; m < n_models; ++m) {
        auto it = mkt.find(m_models[m]);
        MYASSERT(it != mkt.end(), "Missing model in market: " << m_models[m]);
        models[m] = dynamic_cast<IModelMC*>(it->second);
        MYASSERT(models[m] != nullptr, "Market object is not an IModelMC: " << m_models[m]);
        m_dates = models[m]->mergeDates(m_dates);
    }
    std::sort(m_dates.begin(), m_dates.end());
    m_dates.erase(std::unique(m_dates.begin(), m_dates.end()), m_dates.end());
    const size_t n_dates = m_dates.size();
    const size_t n_steps = n_dates - 1;

    // 2. Configure model timelines and bind indices on pricers
    std::vector<size_t>      brownian_offset(n_models + 1, 0);
    std::vector<std::string> factor_names;

    for (size_t m = 0; m < n_models; ++m) {
        models[m]->configureTimeline(m_dates);
        size_t dW              = models[m]->numberOfBrownians();
        brownian_offset[m + 1] = brownian_offset[m] + dW;
        for (size_t b = 0; b < dW; ++b) {
            factor_names.push_back(m_models[m] + "." + std::to_string(b));
        }
    }
    const size_t total_brownians = brownian_offset.back();

    std::vector<const IModelMC*> model_ptrs(models.begin(), models.end());
    for (auto& pricer : m_pricers) {
        pricer->bindIndices(m_dates, m_models, model_ptrs);
    }

    // Pre-fetch settlement discount factors for each pricer
    std::vector<double> df_settle(m_pricers.size(), 1.0);
    for (size_t pr = 0; pr < m_pricers.size(); ++pr) {
        std::string curve_key = "IR.DISCOUNT." + m_pricers[pr]->baseCurrency();
        auto        it        = mkt.find(curve_key);
        if (it != mkt.end()) {
            auto* curve = dynamic_cast<ICurve*>(it->second);
            MYASSERT(curve != nullptr, "Object is not an ICurve: " << curve_key);
            df_settle[pr] = curve->df(m_pricers[pr]->settlementDate());
        }
    }

    // 3. Assemble global M x M correlation matrix Sigma and compute Cholesky factor L
    Matrix Sigma(total_brownians, total_brownians, 1.0);
    for (size_t i = 0; i < total_brownians; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double rho  = getCorrelation(mkt, factor_names[i], factor_names[j]);
            Sigma(i, j) = Sigma(j, i) = rho;
        }
    }
    Matrix L = cholesky(Sigma);

    // 4. Pre-allocate reusable path workspace outside the Monte Carlo loop
    Matrix              rndNumbers(total_brownians, n_steps);
    Matrix              Z_iid(total_brownians, n_steps);
    std::vector<Matrix> stateVariables(n_models);
    for (size_t m = 0; m < n_models; ++m) {
        stateVariables[m] = Matrix(models[m]->numberOfStateVariables(), n_dates);
    }

    std::mt19937_64                  rng(m_seed);
    std::normal_distribution<double> norm(0.0, 1.0);
    std::vector<double>              pv_sum(m_pricers.size(), 0.0);

    auto diffuseAndPrice = [&]() {
        for (size_t m = 0; m < n_models; ++m) {
            ConstMatrixView rnd_slice = rndNumbers.rowSlice(brownian_offset[m], brownian_offset[m + 1]);
            MatrixView      state_view(stateVariables[m]);
            models[m]->diffuse(rnd_slice, state_view);
        }
        for (size_t pr = 0; pr < m_pricers.size(); ++pr) {
            pv_sum[pr] += m_pricers[pr]->pricePath(stateVariables, df_settle[pr], /*fx*/ 1.0);
        }
    };

    size_t p = 0;
    while (p < m_n_paths) {
        for (size_t i = 0; i < total_brownians; ++i) {
            for (size_t k = 0; k < n_steps; ++k) {
                Z_iid(i, k) = norm(rng);
            }
        }
        for (size_t k = 0; k < n_steps; ++k) {
            for (size_t i = 0; i < total_brownians; ++i) {
                double acc = 0.0;
                for (size_t j = 0; j <= i; ++j) {
                    acc += L(i, j) * Z_iid(j, k);
                }
                rndNumbers(i, k) = acc;
            }
        }
        diffuseAndPrice();
        ++p;

        if (m_antithetic && p < m_n_paths) {
            for (size_t i = 0; i < total_brownians; ++i) {
                for (size_t k = 0; k < n_steps; ++k) {
                    rndNumbers(i, k) = -rndNumbers(i, k);
                }
            }
            diffuseAndPrice();
            ++p;
        }
    }

    std::vector<double> pv_results(m_pricers.size(), 0.0);
    for (size_t pr = 0; pr < m_pricers.size(); ++pr) {
        pv_results[pr] = pv_sum[pr] / static_cast<double>(m_n_paths);
    }
    return pv_results;
}

} // namespace mc
