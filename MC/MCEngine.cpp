#include "MCEngine.h"
#include "Models.h"
#include <algorithm>
#include <cmath>
#include <random>

namespace mc {

// Constructor: only structural work. Composite pricers are flattened into leaves here; nothing
// touches the market until price(). That keeps the engine cheap to build and lets the same
// engine be re-priced under several markets (bumps, scenarios).
MCEngine::MCEngine(std::vector<pmcpricer_t> pricers, PricerConfig config)
    : m_config(std::move(config)) {
    MYASSERT(!pricers.empty(), "MCEngine requires at least one pricer");
    MYASSERT(m_config.n_paths > 0, "Number of Monte Carlo paths must be positive");
    m_terms.resize(pricers.size());
    for (size_t i = 0; i < pricers.size(); ++i) {
        MYASSERT(pricers[i] != nullptr, "MCEngine: null pricer");
        // Remember which requested results are composites: only those get a per-component breakdown.
        m_composite.push_back(!pricers[i]->components().empty());
        flatten(pricers[i], 1.0, m_terms[i]);
    }
}

// Depth-first expansion of the pricer tree. A leaf is stored once in m_leaf_pricers even if it appears
// in several results (so it is evaluated once per path); its weight is the product of the weights
// along the path from the root.
void MCEngine::flatten(const pmcpricer_t& p, double weight, std::vector<Term>& terms) {
    auto comps = p->components();
    if (comps.empty()) {
        auto it = std::find(m_leaf_pricers.begin(), m_leaf_pricers.end(), p);
        if (it == m_leaf_pricers.end()) {
            m_leaf_pricers.push_back(p);
            it = m_leaf_pricers.end() - 1;
        }
        terms.push_back({ static_cast<size_t>(it - m_leaf_pricers.begin()), weight });
        return;
    }
    auto w = p->componentWeights();
    MYASSERT(w.size() == comps.size(), "Composite pricer: one weight per component required");
    for (size_t i = 0; i < comps.size(); ++i)
        flatten(comps[i], weight * w[i], terms);
}

namespace {

// Running statistics of one quantity: exact path mean, and standard error computed on
// "draws" (a path and its antithetic counterpart are averaged, as they are not independent)
// Running statistics of one quantity (a leaf, or a requested result).
//   path_sum : sum over ALL paths -> the estimate is path_sum / n_paths (exact path count)
//   cur      : sum over the paths of the draw being processed
//   draw_sum, draw_sq : sum and sum of squares of the DRAW averages
// A "draw" is one set of random numbers: with antithetic variates it produces 2 paths (Z and -Z)
// which are negatively correlated, NOT independent. The standard error must therefore be computed
// on the draw averages, which are independent; computing it on single paths would overstate the
// benefit of antithetic sampling.
struct Acc {
    double path_sum = 0.0;
    double cur      = 0.0; // sum over the paths of the current draw
    double draw_sum = 0.0;
    double draw_sq  = 0.0;

    void endDraw(int n_paths_in_draw) {
        double s = cur / n_paths_in_draw;
        draw_sum += s;
        draw_sq += s * s;
        cur = 0.0;
    }
    double mean(size_t n_paths) const { return path_sum / static_cast<double>(n_paths); }
    // Standard error of the mean of the draw averages: sqrt(sample variance / n_draws),
    // with Bessel's correction n/(n-1) on the variance. The max(0, ...) guards against a tiny negative
    // variance from rounding when the quantity is almost constant.
    double stdErr(size_t n_draws) const {
        if (n_draws < 2)
            return 0.0;
        double n   = static_cast<double>(n_draws);
        double m   = draw_sum / n;
        double var = std::max(0.0, draw_sq / n - m * m) * n / (n - 1.0);
        return std::sqrt(var / n);
    }
};

} // namespace

// The valuation. Read it in eight steps; everything before the path loop is "setup", done once.
std::vector<PriceResult> MCEngine::price(const Market& mkt) const {
    const Date   today    = mkt.today();
    const size_t n_leaves = m_leaf_pricers.size();

    // STEP 1 - Ask every leaf pricer what it needs.
    // requiredDates() are the FUTURE observation dates (the pricer has already removed what is in the
    // events); modelNeeded() are model tags such as "BS_CLOSED.EUR.USD", including the FX model of a
    // foreign payment. Both are collected in plain vectors and then made sorted and unique with
    // sort + erase(unique): a date or a model needed by ten trades appears once, which is what lets all
    // trades share one simulation. `today` is always the first date (the initial state).
    // 1. Master timeline and models required by all leaves
    std::vector<Date>        dates{ today };
    std::vector<std::string> model_tags;
    for (const auto& leaf : m_leaf_pricers) {
        auto d = leaf->requiredDates();
        dates.insert(dates.end(), d.begin(), d.end());
        auto m = leaf->modelNeeded();
        model_tags.insert(model_tags.end(), m.begin(), m.end());
    }
    std::sort(dates.begin(), dates.end());
    dates.erase(std::unique(dates.begin(), dates.end()), dates.end());
    MYASSERT(dates.front() == today, "Payoff required date cannot precede today");
    std::sort(model_tags.begin(), model_tags.end());
    model_tags.erase(std::unique(model_tags.begin(), model_tags.end()), model_tags.end());

    // Tags are sorted "FAMILY.ASSET": two different families for the same asset (e.g. BS_CLOSED.EUR.USD
    // and BS_EULER.EUR.USD) would be two different processes for one asset in a single run, with
    // independent noise: refuse. Mode 1 (one run per trade) avoids this by construction.
    // One simulation scheme per asset
    for (size_t i = 1; i < model_tags.size(); ++i) {
        MYASSERT(assetOfModelTag(model_tags[i]) != assetOfModelTag(model_tags[i - 1]),
                 "Asset simulated with conflicting models in one simulation: "
                     << model_tags[i - 1] << " and " << model_tags[i]
                     << ". Use a consistent PricerConfig, or price the trades separately.");
    }
    const size_t n_models = model_tags.size();

    // STEP 2 - Build the models and finalize the timeline.
    // The market creates each model from its data and the config. A model may then ADD dates it needs
    // for accuracy (an Euler scheme inserts steps of at most 30 days), so the final timeline is the
    // union of the dates required by the payoffs and by the models.
    // 2. Market builds the models; models may add intermediate dates (e.g. Euler steps)
    std::vector<std::shared_ptr<IModelMC>> models;
    for (const auto& tag : model_tags)
        models.push_back(mkt.createModel(tag, m_config));
    for (const auto& m : models)
        dates = m->mergeDates(dates);
    std::sort(dates.begin(), dates.end());
    dates.erase(std::unique(dates.begin(), dates.end()), dates.end());
    const size_t n_dates = dates.size();
    const size_t n_steps = n_dates - 1;

    // STEP 3 - Configure and bind.
    // Each model precomputes its step constants on the final timeline (configureTimeline) and reports
    // how many Brownians it consumes; brownian_offset[m] is the first row of the random matrix owned by
    // model m. factor_names ("EUR.USD.0", ...) label the rows for the correlation lookup.
    // Then every leaf resolves its dates and models into integer positions (bindIndices).
    // 3. Configure model timelines, bind pricers
    std::vector<size_t>      brownian_offset(n_models + 1, 0);
    std::vector<std::string> factor_names; // "<asset>.<brownian index>"
    for (size_t m = 0; m < n_models; ++m) {
        models[m]->configureTimeline(dates);
        size_t dW              = models[m]->numberOfBrownians();
        brownian_offset[m + 1] = brownian_offset[m] + dW;
        for (size_t b = 0; b < dW; ++b)
            factor_names.push_back(assetOfModelTag(model_tags[m]) + "." + std::to_string(b));
    }
    const size_t total_brownians = brownian_offset.back();

    std::vector<const IModelMC*> model_ptrs;
    for (const auto& m : models)
        model_ptrs.push_back(m.get());
    for (const auto& leaf : m_leaf_pricers)
        leaf->bindIndices(dates, model_tags, model_ptrs);

    // STEP 4 - Discount factors (path independent, so computed once, outside the loop).
    // 4. Deterministic part of the discounting of each cash flow, as of today:
    //    scale = P_num(0,Tf) * P_ccy(Tf,Ts) = P_num(0,Tf) * P_ccy(0,Ts) / P_ccy(0,Tf)
    //    A cash flow already paid (Ts < today) is worth nothing. A fixing in the past is
    //    treated as if fixed today.
    std::vector<double> scale(n_leaves, 0.0);
    for (size_t l = 0; l < n_leaves; ++l) {
        // Fixing date, clamped to today: a fixing in the past is known, so there is nothing to discount
        // from before today.
        const Date Tf = std::max(m_leaf_pricers[l]->fixingDate(), today);
        const Date Ts = m_leaf_pricers[l]->settlementDate();
        // The cash flow has already been paid: it is not part of the value any more.
        if (Ts < today)
            continue;
        const std::string ccy = m_leaf_pricers[l]->paymentCurrency();
        // P_USD(0,Tf) discounts the USD value at Tf to today (numeraire).
        // P_ccy(0,Ts) / P_ccy(0,Tf) = P_ccy(Tf,Ts) is the forward discount factor from payment to fixing in
        // the payment currency: the conditional expectation of the cash flow on its fixing date.
        // The pricer will multiply by the contingent FX spot (ccy.USD at Tf) when ccy != USD.
        scale[l]              = mkt.df(NUMERAIRE_CCY, Tf) * mkt.df(ccy, Ts) / mkt.df(ccy, Tf);
    }

    // STEP 5 - Correlation. One matrix Sigma over ALL Brownians of ALL models (market lookup by
    // risk-factor name, default 0), then Sigma = L L^T by Cholesky. Correlated normals are L * Z with
    // Z independent standard normals; the same matrix therefore covers any mix of models.
    // 5. Global correlation matrix and Cholesky factor
    Matrix Sigma(total_brownians, total_brownians, 1.0);
    for (size_t i = 0; i < total_brownians; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double rho  = mkt.correlation(factor_names[i], factor_names[j]);
            Sigma(i, j) = Sigma(j, i) = rho;
        }
    }
    Matrix L = cholesky(Sigma);

    // STEP 6 - Allocate the work space once. Re-using it on every path avoids millions of allocations.
    //   rndNumbers, Z_iid : [total_brownians x n_steps], the correlated / independent normals of ONE
    //                       path; row = one Brownian, column = one step. The rows of model m are
    //                       brownian_offset[m] .. brownian_offset[m+1]-1.
    //   stateVariables[m] : [numberOfStateVariables(m) x n_dates], the state history of model m on the
    //                       current path (column 0 = today). Pricers read it through observable().
    //   prev_states[m], next_states[m] : the states of model m at the current / next date, each
    //                       numberOfStateVariables(m) values (the arguments of generateStates).
    //   step_rnd[m]       : the normals of the current step for model m, numberOfBrownians(m) values,
    //                       copied from a column of rndNumbers because that column is not contiguous.
    // 6. Workspace
    Matrix              rndNumbers(total_brownians, n_steps);
    Matrix              Z_iid(total_brownians, n_steps);
    std::vector<Matrix> stateVariables(n_models);
    std::vector<std::vector<double>> prev_states(n_models), next_states(n_models), step_rnd(n_models);
    for (size_t m = 0; m < n_models; ++m) {
        stateVariables[m] = Matrix(models[m]->numberOfStateVariables(), n_dates);
        prev_states[m].assign(models[m]->numberOfStateVariables(), 0.0);
        next_states[m].assign(models[m]->numberOfStateVariables(), 0.0);
        step_rnd[m].assign(models[m]->numberOfBrownians(), 0.0);
    }

    std::vector<Acc>    leaf_acc(n_leaves);
    std::vector<Acc>    root_acc(m_terms.size());
    std::vector<double> leaf_val(n_leaves, 0.0);
    size_t              n_draws = 0;

    // One path: advance every model step by step, price every leaf on the resulting states, then
    // recombine the leaves into each requested result.
    // Evaluating all leaves on the same states is what makes a portfolio a single simulation.
    auto diffuseAndPrice = [&]() {
        // Date 0 (today): every model sets its initial states, which are also stored in column 0.
        for (size_t m = 0; m < n_models; ++m) {
            models[m]->initializeStates(prev_states[m]);
            for (size_t s = 0; s < prev_states[m].size(); ++s)
                stateVariables[m](s, 0) = prev_states[m][s];
        }
        // Step k goes from date k to date k+1. The loop over steps is OUTSIDE the loop over models:
        // all models move together, one step at a time, so a model could later depend on the
        // states of another one at the same date (e.g. stochastic rates) without changing the engine.
        for (size_t k = 0; k < n_steps; ++k) {
            for (size_t m = 0; m < n_models; ++m) {
                // The model's normals of step k: column k of its rows of rndNumbers.
                for (size_t b = 0; b < step_rnd[m].size(); ++b)
                    step_rnd[m][b] = rndNumbers(brownian_offset[m] + b, k);
                // previous states + random numbers -> new states
                models[m]->generateStates(k, step_rnd[m], prev_states[m], next_states[m]);
                for (size_t s = 0; s < next_states[m].size(); ++s)
                    stateVariables[m](s, k + 1) = next_states[m][s];
                std::swap(prev_states[m], next_states[m]);   // the new states are the next step's previous ones
            }
        }
        for (size_t l = 0; l < n_leaves; ++l) {
            // The pricer returns its value in USD as of today (it applies `scale` and the contingent FX).
            leaf_val[l] = m_leaf_pricers[l]->pricePath(stateVariables, scale[l]);
            leaf_acc[l].path_sum += leaf_val[l];
            leaf_acc[l].cur += leaf_val[l];
        }
        for (size_t r = 0; r < m_terms.size(); ++r) {
            double x = 0.0;
            for (const auto& t : m_terms[r])
                // A requested result is a weighted sum of leaves (weight 1 for a plain trade).
                x += t.weight * leaf_val[t.leaf];
            root_acc[r].path_sum += x;
            root_acc[r].cur += x;
        }
    };

    // STEP 7 - The path loop.
    // The generator is seeded from the config at EVERY call: identical inputs give identical outputs,
    // and two runs that differ only by a market bump see the same random numbers (Common Random
    // Numbers), so the bump is not drowned by sampling noise.
    // 7. Monte Carlo path loop
    std::mt19937_64                  rng(m_config.seed);
    std::normal_distribution<double> norm(0.0, 1.0);
    size_t                           p = 0;
    while (p < m_config.n_paths) {
        for (size_t i = 0; i < total_brownians; ++i)
            for (size_t k = 0; k < n_steps; ++k)
                Z_iid(i, k) = norm(rng);

        // L is lower triangular, so row i only needs the first i+1 independent normals.
        // Correlate: rndNumbers(:, k) = L * Z_iid(:, k)
        for (size_t k = 0; k < n_steps; ++k) {
            for (size_t i = 0; i < total_brownians; ++i) {
                double acc = 0.0;
                for (size_t j = 0; j <= i; ++j)
                    acc += L(i, j) * Z_iid(j, k);
                rndNumbers(i, k) = acc;
            }
        }

        diffuseAndPrice();
        ++p;
        int in_draw = 1;

        // Antithetic path: same draw with the sign flipped. For payoffs monotone in the noise the two
        // payoffs are negatively correlated, which lowers the variance of their average. The guard on
        // `p` keeps the total at exactly n_paths also when n_paths is odd.
        if (m_config.antithetic && p < m_config.n_paths) {
            for (size_t i = 0; i < total_brownians; ++i)
                for (size_t k = 0; k < n_steps; ++k)
                    rndNumbers(i, k) = -rndNumbers(i, k);
            diffuseAndPrice();
            ++p;
            in_draw = 2;
        }
        for (auto& a : leaf_acc) a.endDraw(in_draw);
        for (auto& a : root_acc) a.endDraw(in_draw);
        ++n_draws;
    }

    // STEP 8 - Results. PV = mean over all paths; standard error from the draw averages.
    // For a composite the breakdown lists each leaf's weighted PV and |weight| * its standard error;
    // the components add up to the total exactly, while their errors do not (the leaves are
    // correlated, so the total's standard error is smaller than the sum of theirs).
    // 8. Results: one per input pricer; composites carry the breakdown
    std::vector<PriceResult> results(m_terms.size());
    for (size_t r = 0; r < m_terms.size(); ++r) {
        results[r].pv      = root_acc[r].mean(m_config.n_paths);
        results[r].std_err = root_acc[r].stdErr(n_draws);
        if (m_composite[r]) {
            for (const auto& t : m_terms[r]) {
                PriceResult c;
                c.pv      = t.weight * leaf_acc[t.leaf].mean(m_config.n_paths);
                c.std_err = std::abs(t.weight) * leaf_acc[t.leaf].stdErr(n_draws);
                results[r].components.push_back(c);
            }
        }
    }
    return results;
}

} // namespace mc
