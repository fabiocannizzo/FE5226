#pragma once

#include "Interfaces.h"
#include <algorithm>
#include <cmath>
#include <memory>

namespace mc {

// Pairwise Risk-Factor Correlation Object stored in the Market Dictionary
class CorrelationObject : public ICorrelation {
public:
    CorrelationObject(std::string name, double rho)
        : m_name(std::move(name)), m_rho(rho) {
        MYASSERT(rho >= -1.0 && rho <= 1.0, "Correlation must be in [-1, 1]");
    }

    std::string name() const override  { return m_name; }
    double      value() const override { return m_rho; }

private:
    std::string m_name;
    double      m_rho;
};

// Flat Continuously-Compounded Discount Curve: DF(T) = exp(-r * (T - T_0) / 365)
class FlatDiscountCurve : public ICurve {
public:
    FlatDiscountCurve(std::string name, Date today, double rate)
        : m_name(std::move(name)), m_today(today), m_rate(rate) {}

    std::string name() const override { return m_name; }

    double df(Date d) const override {
        double tau = static_cast<double>(d - m_today) / 365.0;
        return std::exp(-m_rate * tau);
    }

private:
    std::string m_name;
    Date        m_today;
    double      m_rate;
};

// ----------------------------------------------------------------------------
// Model 1: Black-Scholes Closed-Form Model ("BS_CLOSED.<PAIR>")
//   - numberOfBrownians()      = 1
//   - numberOfStateVariables() = 1 (Brownian motion W(t))
//   - Input rndNumbers         = correlated normal increments
//   - Output stateVariables    = standard Brownian motion W(t_k)
//   - Reconstruction           = recovers spot S(t) = S_0 * exp(A(t) + sigma * W(t))
//   - mergeDates()             = returns timeline as-is (no discretization error)
// ----------------------------------------------------------------------------
class ModelBSClosedForm : public IModelMC {
public:
    ModelBSClosedForm(std::string name, double spot0, double mu, double sigma)
        : m_name(std::move(name))
        , m_spot0(spot0)
        , m_mu(mu)
        , m_sigma(sigma) {}

    std::string name() const override { return m_name; }

    size_t numberOfBrownians() const override      { return 1; }
    size_t numberOfStateVariables() const override { return 1; } // State 0 is Brownian W(t)

    double spot0() const { return m_spot0; }
    void   setSpot0(double s) { m_spot0 = s; }
    double sigma() const { return m_sigma; }
    double mu() const    { return m_mu; }

    // Exact transition: no intermediate discretization dates needed
    std::vector<Date> mergeDates(const std::vector<Date>& dates) const override {
        return dates;
    }

    void configureTimeline(const std::vector<Date>& dates) override {
        MYASSERT(dates.size() >= 2, "Timeline must have at least 2 dates");
        m_dates = dates;
        size_t n_steps = dates.size() - 1;
        m_sqrt_dt.resize(n_steps);
        m_drift_accum.resize(dates.size());

        m_drift_accum[0] = 0.0;
        for (size_t k = 0; k < n_steps; ++k) {
            double dt            = static_cast<double>(dates[k + 1] - dates[k]) / 365.0;
            m_sqrt_dt[k]         = std::sqrt(dt);
            double t             = static_cast<double>(dates[k + 1] - dates[0]) / 365.0;
            m_drift_accum[k + 1] = (m_mu - 0.5 * m_sigma * m_sigma) * t;
        }
    }

    // Diffusion: propagates Brownian motion W(t_k)
    void diffuse(const ConstMatrixView& rndNumbers,
                 MatrixView&            stateVariables) const override {
        double W = 0.0;
        stateVariables(0, 0) = W; // W(0) = 0
        for (size_t k = 0; k < m_sqrt_dt.size(); ++k) {
            W += m_sqrt_dt[k] * rndNumbers(0, k);
            stateVariables(0, k + 1) = W;
        }
    }

    // Observable function: queries model passing stateVariables to reconstruct Spot S(t)
    double observable(size_t date_idx, const Matrix& stateVariables) const override {
        double W = stateVariables(0, date_idx);
        return m_spot0 * std::exp(m_drift_accum[date_idx] + m_sigma * W);
    }

private:
    std::string         m_name;
    double              m_spot0;
    double              m_mu;
    double              m_sigma;
    std::vector<Date>   m_dates;
    std::vector<double> m_sqrt_dt;
    std::vector<double> m_drift_accum;
};

using ModelBSFX = ModelBSClosedForm;

// ----------------------------------------------------------------------------
// Model 2: Black-Scholes Euler Discretization ("BS_EULER.<PAIR>")
//   - numberOfBrownians()      = 1
//   - numberOfStateVariables() = 1 (Spot S(t))
//   - Input rndNumbers         = correlated normal increments
//   - Output stateVariables    = spot price S(t_k)
//   - Reconstruction           = state variable IS Spot!
//   - mergeDates()             = injects time steps so step size <= max_step_days (1 month = 30d)
// ----------------------------------------------------------------------------
class ModelBSEuler : public IModelMC {
public:
    ModelBSEuler(std::string name, double spot0, double mu, double sigma, int max_step_days = 30)
        : m_name(std::move(name))
        , m_spot0(spot0)
        , m_mu(mu)
        , m_sigma(sigma)
        , m_max_step_days(max_step_days) {}

    std::string name() const override { return m_name; }

    size_t numberOfBrownians() const override      { return 1; }
    size_t numberOfStateVariables() const override { return 1; } // State 0 is Spot S(t)

    double spot0() const { return m_spot0; }
    void   setSpot0(double s) { m_spot0 = s; }
    double sigma() const { return m_sigma; }
    double mu() const    { return m_mu; }

    // Inject intermediate dates so that no step is larger than 1 month (max_step_days = 30)
    std::vector<Date> mergeDates(const std::vector<Date>& dates) const override {
        if (dates.size() < 2)
            return dates;
        std::vector<Date> merged;
        merged.push_back(dates.front());
        for (size_t i = 1; i < dates.size(); ++i) {
            Date prev   = merged.back();
            Date target = dates[i];
            while (target - prev > m_max_step_days) {
                prev += m_max_step_days;
                merged.push_back(prev);
            }
            if (merged.back() != target) {
                merged.push_back(target);
            }
        }
        return merged;
    }

    void configureTimeline(const std::vector<Date>& dates) override {
        MYASSERT(dates.size() >= 2, "Timeline must have at least 2 dates");
        size_t n_steps = dates.size() - 1;
        m_dt.resize(n_steps);
        m_sqrt_dt.resize(n_steps);
        for (size_t k = 0; k < n_steps; ++k) {
            double dt      = static_cast<double>(dates[k + 1] - dates[k]) / 365.0;
            m_dt[k]        = dt;
            m_sqrt_dt[k]   = std::sqrt(dt);
        }
    }

    // Diffusion: propagates Spot S(t_k) directly using Euler-Maruyama
    void diffuse(const ConstMatrixView& rndNumbers,
                 MatrixView&            stateVariables) const override {
        double S = m_spot0;
        stateVariables(0, 0) = S;
        for (size_t k = 0; k < m_dt.size(); ++k) {
            S += m_mu * S * m_dt[k] + m_sigma * S * m_sqrt_dt[k] * rndNumbers(0, k);
            if (S < 0.0)
                S = 0.0; // Truncation
            stateVariables(0, k + 1) = S;
        }
    }

    // Observable function: state variable IS Spot!
    double observable(size_t date_idx, const Matrix& stateVariables) const override {
        return stateVariables(0, date_idx);
    }

private:
    std::string         m_name;
    double              m_spot0;
    double              m_mu;
    double              m_sigma;
    int                 m_max_step_days;
    std::vector<double> m_dt;
    std::vector<double> m_sqrt_dt;
};

// ----------------------------------------------------------------------------
// Model 3: Heston Stochastic Volatility Model ("HESTON.<PAIR>")
//   - numberOfBrownians()      = 2 (dW^S and dW^v)
//   - numberOfStateVariables() = 2 (Spot S_t and Variance v_t)
// ----------------------------------------------------------------------------
class ModelHestonFX : public IModelMC {
public:
    ModelHestonFX(std::string name,
                  double      spot0,
                  double      v0,
                  double      kappa,
                  double      theta,
                  double      xi,
                  double      mu)
        : m_name(std::move(name))
        , m_spot0(spot0)
        , m_v0(v0)
        , m_kappa(kappa)
        , m_theta(theta)
        , m_xi(xi)
        , m_mu(mu) {}

    std::string name() const override { return m_name; }

    size_t numberOfBrownians() const override      { return 2; }
    size_t numberOfStateVariables() const override { return 2; } // State 0: Spot, State 1: Variance

    void configureTimeline(const std::vector<Date>& dates) override {
        MYASSERT(dates.size() >= 2, "Timeline must have at least 2 dates");
        size_t n_steps = dates.size() - 1;
        m_dt.resize(n_steps);
        m_sqrt_dt.resize(n_steps);
        for (size_t k = 0; k < n_steps; ++k) {
            double dt    = static_cast<double>(dates[k + 1] - dates[k]) / 365.0;
            m_dt[k]      = dt;
            m_sqrt_dt[k] = std::sqrt(dt);
        }
    }

    void diffuse(const ConstMatrixView& rndNumbers,
                 MatrixView&            stateVariables) const override {
        double S = m_spot0;
        double v = m_v0;
        stateVariables(0, 0) = S;
        stateVariables(1, 0) = v;

        for (size_t k = 0; k < m_dt.size(); ++k) {
            double dt      = m_dt[k];
            double sdt     = m_sqrt_dt[k];
            double v_pos   = std::max(v, 0.0);
            double sqrt_v  = std::sqrt(v_pos);

            // Euler-Full Truncation for variance
            v += m_kappa * (m_theta - v_pos) * dt + m_xi * sqrt_v * sdt * rndNumbers(1, k);
            // Log-Euler for spot
            S *= std::exp((m_mu - 0.5 * v_pos) * dt + sqrt_v * sdt * rndNumbers(0, k));

            stateVariables(0, k + 1) = S;
            stateVariables(1, k + 1) = v;
        }
    }

    double observable(size_t date_idx, const Matrix& stateVariables) const override {
        return stateVariables(0, date_idx); // State 0 is Spot
    }

private:
    std::string         m_name;
    double              m_spot0;
    double              m_v0;
    double              m_kappa;
    double              m_theta;
    double              m_xi;
    double              m_mu;
    std::vector<double> m_dt;
    std::vector<double> m_sqrt_dt;
};

// ----------------------------------------------------------------------------
// Model Factory: Simple switch() statement constructing IModelMC instances
// given model tag (e.g. "BS_CLOSED.EUR.USD") and MarketDataServer
// ----------------------------------------------------------------------------
enum class ModelFamily {
    BS_CLOSED,
    BS_EULER,
    HESTON,
    UNKNOWN
};

inline ModelFamily parseModelFamily(const std::string& family) {
    if (family == "BS_CLOSED" || family == "BS")
        return ModelFamily::BS_CLOSED;
    if (family == "BS_EULER")
        return ModelFamily::BS_EULER;
    if (family == "HESTON")
        return ModelFamily::HESTON;
    return ModelFamily::UNKNOWN;
}

inline std::shared_ptr<IModelMC> createModelMC(const std::string&      model_tag,
                                               const MarketDataServer& mkt) {
    auto dot = model_tag.find('.');
    MYASSERT(dot != std::string::npos, "Invalid model tag format: " << model_tag);
    std::string family = model_tag.substr(0, dot);
    std::string asset  = model_tag.substr(dot + 1);

    double S0    = mkt.getSpot(asset);
    double mu    = mkt.getDrift(asset);
    double sigma = mkt.getVol(asset);

    switch (parseModelFamily(family)) {
    case ModelFamily::BS_CLOSED:
        return std::make_shared<ModelBSClosedForm>(model_tag, S0, mu, sigma);
    case ModelFamily::BS_EULER:
        return std::make_shared<ModelBSEuler>(model_tag, S0, mu, sigma, /*max_step_days*/ 30);
    case ModelFamily::HESTON:
        return std::make_shared<ModelHestonFX>(model_tag, S0, sigma * sigma, 2.0, sigma * sigma, 0.25, mu);
    default:
        MYASSERT(false, "Unknown model family in factory: " << family);
        return nullptr;
    }
}

inline std::vector<std::shared_ptr<IModelMC>> createModelsMC(
    const std::vector<std::string>& models_needed,
    const MarketDataServer&         mkt) {
    std::vector<std::shared_ptr<IModelMC>> models;
    models.reserve(models_needed.size());
    for (const auto& tag : models_needed) {
        models.push_back(createModelMC(tag, mkt));
    }
    return models;
}

} // namespace mc
