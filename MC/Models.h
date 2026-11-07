#pragma once

#include "Interfaces.h"
#include <algorithm>
#include <cmath>
#include <memory>

namespace mc {

// Pairwise Risk-Factor Correlation Object stored in the Market Dictionary
// Value object holding one correlation coefficient. (Market keeps plain doubles; this class
// shows how a correlation could be a first-class, named market object.)
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
// The simplest curve: one continuously compounded rate r for all maturities,
// P(0, d) = exp(-r * tau),  tau = (d - today) / 365.
// A real curve (bootstrapped from deposits / swaps) would implement the same ICurve interface
// and nothing else in the library would change.
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
//   - Input rndNumbers (step)  = correlated normal of the step
//   - States                   = standard Brownian motion W(t_k)
//   - Reconstruction           = recovers spot S(t) = S_0 * exp(A(t) + sigma * W(t))
//   - mergeDates()             = returns timeline as-is (no discretization error)
// ----------------------------------------------------------------------------
// Geometric Brownian motion, simulated EXACTLY:
//     dS = mu S dt + sigma S dW   =>   S(t) = S0 exp((mu - sigma^2/2) t + sigma W(t))
// The state variable is the Brownian motion W(t), not the spot. Since W has independent normal
// increments, W(t_k) is built exactly on ANY set of dates: no discretization error, and no extra
// dates are needed. The spot is recovered on demand by observable().
// `mu` is the risk neutral drift (for FX: domestic minus foreign rate); `sigma` the flat vol.
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
    // Closed form: any set of dates is fine, nothing to add.
    std::vector<Date> mergeDates(const std::vector<Date>& dates) const override {
        return dates;
    }

    // Called once by the engine, after the master timeline is final. Precompute everything that does
    // not depend on the random numbers, so the per-path work is a few multiplications.
    // Here: sqrt(dt) of each step, and the deterministic part of the exponent at each date.
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
            // Deterministic part of log(S/S0): (mu - sigma^2/2) * t, t in years since the first date.
            // The -sigma^2/2 "Ito correction" makes E[S(t)] = S0 exp(mu t), i.e. the forward.
            m_drift_accum[k + 1] = (m_mu - 0.5 * m_sigma * m_sigma) * t;
        }
    }

    // The only state is the Brownian motion W (not the spot), and W(today) = 0.
    void initializeStates(std::span<double> states) const override {
        states[0] = 0.0;
    }

    // One step of the Brownian motion:  W(t_{k+1}) = W(t_k) + sqrt(dt_k) * Z_k
    //   rndNumbers[0]    = Z_k, the correlated standard normal of step k (this model has 1 Brownian)
    //   previousStates[0]= W(t_k)
    //   states[0]        = W(t_{k+1})   (output)
    // Exact for any step length, hence no discretization error.
    void generateStates(size_t                  timeIndex,
                        std::span<const double> rndNumbers,
                        std::span<const double> previousStates,
                        std::span<double>       states) const override {
        states[0] = previousStates[0] + m_sqrt_dt[timeIndex] * rndNumbers[0];
    }

    // Observable function: queries model passing stateVariables to reconstruct Spot S(t)
    // Reconstruction: S(t) = S0 * exp(drift_accum(t) + sigma * W(t)). Computed on demand, only
    // for the dates and trades that need it.
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

// Historical alias used in earlier versions of the example.
using ModelBSFX = ModelBSClosedForm;

// ----------------------------------------------------------------------------
// Model 2: Black-Scholes Euler Discretization ("BS_EULER.<PAIR>")
//   - numberOfBrownians()      = 1
//   - numberOfStateVariables() = 1 (Spot S(t))
//   - Input rndNumbers (step)  = correlated normal of the step
//   - States                   = spot price S(t_k)
//   - Reconstruction           = state variable IS Spot!
//   - mergeDates()             = injects time steps so step size <= max_step_days (1 month = 30d)
// ----------------------------------------------------------------------------
// Same dynamics as ModelBSClosedForm but discretized with the Euler-Maruyama scheme,
//     S_{k+1} = S_k + mu S_k dt + sigma S_k sqrt(dt) Z_k,
// which is only an approximation: the error shrinks with the step size, hence mergeDates() adds
// intermediate dates so that no step exceeds max_step_days. The state variable IS the spot.
// It is here to show a model that is NOT exact (as local vol, ... would be), next to the
// exact one, and to measure the discretization error against the closed form.
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
    // Insert intermediate dates so that no step is longer than m_max_step_days (default 30 days).
    // Example with the default: dates {0, 100} become {0, 30, 60, 90, 100}.
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

    // Precompute dt and sqrt(dt) for every step of the (already merged) timeline.
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

    // The only state is the spot itself, starting from today's market spot.
    void initializeStates(std::span<double> states) const override {
        states[0] = m_spot0;
    }

    // One Euler-Maruyama step:  S_{k+1} = S_k + mu S_k dt + sigma S_k sqrt(dt) Z_k
    //   rndNumbers[0]    = Z_k, the correlated standard normal of step k (1 Brownian)
    //   previousStates[0]= S_k
    //   states[0]        = S_{k+1}   (output)
    // The step length comes from the timeline, which mergeDates() has made short enough.
    void generateStates(size_t                  timeIndex,
                        std::span<const double> rndNumbers,
                        std::span<const double> previousStates,
                        std::span<double>       states) const override {
        const double S0 = previousStates[0];
        // Deterministic drift plus volatility times the normal draw.
        double S = S0 + (m_mu * S0 * m_dt[timeIndex] + m_sigma * S0 * m_sqrt_dt[timeIndex] * rndNumbers[0]);
        // The scheme can overshoot below zero for large steps or vols; truncate to keep the spot valid.
        if (S < 0.0)
            S = 0.0; // Truncation
        states[0] = S;
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
// Model families. Models are built by Market::createModel (see Market.h) from a
// model tag (e.g. "BS_CLOSED.EUR.USD") and the PricerConfig
// ----------------------------------------------------------------------------
// Model families selectable by name in PricerConfig. The name is the first part of a model tag
// "<FAMILY>.<ASSET>" such as "BS_EULER.EUR.USD"; Market::createModel() builds the model.
enum class ModelFamily {
    BS_CLOSED,
    BS_EULER,
    UNKNOWN
};

// Text -> enum. "BS" is accepted as an alias of "BS_CLOSED". Unknown names give UNKNOWN, which the
// factory reports as an error.
inline ModelFamily parseModelFamily(const std::string& family) {
    if (family == "BS_CLOSED" || family == "BS")
        return ModelFamily::BS_CLOSED;
    if (family == "BS_EULER")
        return ModelFamily::BS_EULER;
    return ModelFamily::UNKNOWN;
}

} // namespace mc
