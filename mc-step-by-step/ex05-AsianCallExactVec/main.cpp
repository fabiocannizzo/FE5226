// ex05 - Asian call: MULTI-STEP paths, exact simulation, vectorized.
//
// The payoff now depends on the spot at several dates t_1 < ... < t_m (the fixing dates):
//     arithmetic Asian call :  max( (S(t_1)+...+S(t_m))/m - K , 0 )
//
// Same architecture as ex04 (model object, payoff object), generalized to several dates:
//   - the payoff TELLS the engine the dates it needs (dates()); they are the simulation timeline
//   - the engine gives the timeline to the model (configureTimeline), then advances the model
//     step by step: the states of date j are one Matrix (state variables x paths), and the states
//     of all dates are a vector of matrices
//   - the model still returns the spot at ONE date, as a vector over the paths: spot(j, ...).
//     The payoff queries the model once per date it needs, referring to the dates BY INDEX j
//     (index j = the j-th date of the payoff)
//
// Exact simulation on any grid, through the Brownian motion W (independent increments):
//     W(t_j) = W(t_{j-1}) + sqrt(t_j - t_{j-1}) Z_j
//     S(t_j) = S0 exp( (r - q - sigma^2/2) t_j + sigma W(t_j) )
//
// Validation: the GEOMETRIC Asian call has a closed form, and a single fixing date gives back
// the European call of ex01-ex04.
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <span>
#include <utility>
#include <vector>

struct Market { double S0, r, q, sigma; };
struct McResult { double price, stdErr; };

// Row-major matrix: element (i,j) is in data[i*cols + j]
struct Matrix {
    Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), data(rows * cols) {}
    double&              operator()(size_t i, size_t j)       { return data[i * cols + j]; }
    double               operator()(size_t i, size_t j) const { return data[i * cols + j]; }
    std::span<double>       row(size_t i)       { return {data.data() + i * cols, cols}; }
    std::span<const double> row(size_t i) const { return {data.data() + i * cols, cols}; }
    size_t rows, cols;
    std::vector<double> data;
};

// States of all the dates: states[j] holds (state variables x paths) at the j-th date
using States = std::vector<Matrix>;

struct IModelMC {
    virtual ~IModelMC() = default;

    virtual size_t numberOfBrownians() const = 0;
    virtual size_t numberOfStateVariables() const = 0;

    // The engine tells the model the dates of the simulation (years from today, increasing,
    // all > 0): date j is dates[j]. The model precomputes whatever it needs.
    virtual void configureTimeline(const std::vector<double>& dates) = 0;

    // states at today, for the first n paths
    virtual void initializeStates(Matrix& states, size_t n) const = 0;

    // One step, from the date j-1 (today if j == 0) to the date j, for the first n paths:
    //   rnd : Brownians x paths          independent N(0,1) draws of the step
    //   prev: states x paths             states at the date j-1
    //   next: states x paths             states at the date j (output)
    virtual void generateStates(size_t j, const Matrix& rnd, const Matrix& prev,
                                Matrix& next, size_t n) const = 0;

    // Spot at the date j, for the first n paths, rebuilt from states[j]: out[i] = S(t_j) on path i
    virtual void spot(size_t j, const States& states, size_t n, std::span<double> out) const = 0;
};

// Black-Scholes, exact on any grid. State variable: the Brownian motion W.
struct BSExact : IModelMC {
    explicit BSExact(const Market& m) : m_spot0(m.S0), m_mu(m.r - m.q), m_sigma(m.sigma) {}

    size_t numberOfBrownians() const override      { return 1; }
    size_t numberOfStateVariables() const override { return 1; }

    void configureTimeline(const std::vector<double>& dates) override {
        m_sqrtDt.resize(dates.size());
        m_drift.resize(dates.size());
        for (size_t j = 0; j < dates.size(); ++j) {
            const double tPrev = (j == 0) ? 0.0 : dates[j - 1];
            m_sqrtDt[j] = std::sqrt(dates[j] - tPrev);
            m_drift[j]  = (m_mu - 0.5 * m_sigma * m_sigma) * dates[j];
        }
    }

    void initializeStates(Matrix& states, size_t n) const override {
        for (size_t i = 0; i < n; ++i)
            states(0, i) = 0.0;                          // W(0) = 0
    }

    void generateStates(size_t j, const Matrix& rnd, const Matrix& prev,
                        Matrix& next, size_t n) const override {
        for (size_t i = 0; i < n; ++i)
            next(0, i) = prev(0, i) + m_sqrtDt[j] * rnd(0, i);
    }

    void spot(size_t j, const States& states, size_t n, std::span<double> out) const override {
        for (size_t i = 0; i < n; ++i)
            out[i] = m_spot0 * std::exp(m_drift[j] + m_sigma * states[j](0, i));
    }

    double m_spot0, m_mu, m_sigma;
    std::vector<double> m_sqrtDt, m_drift;   // per step, computed once by configureTimeline
};

struct IPayoff {
    virtual ~IPayoff() = default;
    // Dates (years from today, increasing) at which the payoff needs the spot; paid at the last
    virtual const std::vector<double>& dates() const = 0;
    // out[i] = payoff on path i, for the first n paths. The payoff queries the model for the
    // spot at each of its dates, by index: model.spot(j, states, n, ...) is the j-th date.
    virtual void payoff(const IModelMC& model, const States& states, size_t n,
                        std::span<double> out) const = 0;
};

McResult mcPrice(const Market& m, IModelMC& model, const IPayoff& payoff,
                 size_t nPaths, size_t batchSize, unsigned seed)
{
    const std::vector<double>& dates = payoff.dates();   // asked to the payoff
    const size_t nDates = dates.size();
    const double T      = dates.back();
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    model.configureTimeline(dates);                      // the model learns the timeline

    // Work space allocated once, re-used by every batch
    const size_t nS = model.numberOfStateVariables(), nB = model.numberOfBrownians();
    Matrix rnd(nB, batchSize), initial(nS, batchSize);
    States states(nDates, Matrix(nS, batchSize));        // states[j]: the states at the date j
    std::vector<double> pay(batchSize);

    double sum = 0.0, sumSq = 0.0;
    for (size_t done = 0; done < nPaths; done += batchSize) {
        const size_t n = std::min(batchSize, nPaths - done);   // the last batch may be shorter

        model.initializeStates(initial, n);
        for (size_t j = 0; j < nDates; ++j) {                   // advance all paths by one step
            for (size_t b = 0; b < nB; ++b)
                for (size_t i = 0; i < n; ++i)
                    rnd(b, i) = normal(rng);
            model.generateStates(j, rnd, j == 0 ? initial : states[j - 1], states[j], n);
        }

        payoff.payoff(model, states, n, std::span<double>(pay).first(n));
        for (size_t i = 0; i < n; ++i) {
            sum   += pay[i];
            sumSq += pay[i] * pay[i];
        }
    }

    const double N    = static_cast<double>(nPaths);
    const double mean = sum / N;
    const double var  = (sumSq / N - mean * mean) * N / (N - 1.0);
    const double df   = std::exp(-m.r * T);
    return {df * mean, df * std::sqrt(var / N)};
}

// ---- payoffs ----
// Arithmetic Asian call: max( average of the spots at the fixing dates - K, 0 )
struct PayoffAsianCall : IPayoff {
    PayoffAsianCall(double strike, std::vector<double> fixingDates)
        : m_strike(strike), m_dates(std::move(fixingDates)) {}

    const std::vector<double>& dates() const override { return m_dates; }

    void payoff(const IModelMC& model, const States& states, size_t n,
                std::span<double> out) const override {
        std::vector<double> spot(n);
        std::fill(out.begin(), out.end(), 0.0);
        for (size_t j = 0; j < m_dates.size(); ++j) {   // one query to the model per date
            model.spot(j, states, n, spot);
            for (size_t i = 0; i < n; ++i)
                out[i] += spot[i];
        }
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(out[i] / static_cast<double>(m_dates.size()) - m_strike, 0.0);
    }

    double m_strike;
    std::vector<double> m_dates;
};

// Geometric Asian call: max( geometric average of the spots - K, 0 ), has a closed form
struct PayoffGeometricAsianCall : IPayoff {
    PayoffGeometricAsianCall(double strike, std::vector<double> fixingDates)
        : m_strike(strike), m_dates(std::move(fixingDates)) {}

    const std::vector<double>& dates() const override { return m_dates; }

    void payoff(const IModelMC& model, const States& states, size_t n,
                std::span<double> out) const override {
        std::vector<double> spot(n);
        std::fill(out.begin(), out.end(), 0.0);          // sum of log-spots
        for (size_t j = 0; j < m_dates.size(); ++j) {
            model.spot(j, states, n, spot);
            for (size_t i = 0; i < n; ++i)
                out[i] += std::log(spot[i]);
        }
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(std::exp(out[i] / static_cast<double>(m_dates.size())) - m_strike, 0.0);
    }

    double m_strike;
    std::vector<double> m_dates;
};

// ---- closed forms ----
double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// E[max(X-K,0)] for X lognormal with E[log X] = M and Var[log X] = V, discounted
double lognormalCall(double M, double V, double K, double df)
{
    const double sd = std::sqrt(V);
    const double d1 = (M + V - std::log(K)) / sd;
    return df * (std::exp(M + 0.5 * V) * normCdf(d1) - K * normCdf(d1 - sd));
}

// Geometric average G = (prod S(t_j))^(1/m): log G is normal.
double geometricAsianCallExact(const Market& m, const std::vector<double>& t, double K)
{
    const size_t n = t.size();
    double mean = 0.0, var = 0.0;
    for (size_t i = 0; i < n; ++i) {
        mean += t[i];
        for (size_t j = 0; j < n; ++j)
            var += std::min(t[i], t[j]);                // Cov(W(t_i), W(t_j)) = min(t_i, t_j)
    }
    const double nn = static_cast<double>(n);
    const double M  = std::log(m.S0) + (m.r - m.q - 0.5 * m.sigma * m.sigma) * mean / nn;
    const double V  = m.sigma * m.sigma * var / (nn * nn);
    return lognormalCall(M, V, K, std::exp(-m.r * t.back()));
}

int main()
{
    const Market m{100.0, 0.03, 0.01, 0.20};
    const double K = 100.0;
    BSExact model(m);
    const size_t nPaths = 1'000'000, batch = 4096;
    int failures = 0;

    auto check = [&](const char* what, const McResult& mc, double exact) {
        const double z = (mc.price - exact) / mc.stdErr;
        std::cout << what << ": MC = " << mc.price << " +/- " << mc.stdErr
                  << "  exact = " << exact << "  (" << z << " std errors)\n";
        if (std::abs(z) > 4.0) ++failures;
    };

    // Monthly fixings over one year
    std::vector<double> monthly;
    for (int j = 1; j <= 12; ++j) monthly.push_back(j / 12.0);

    check("geometric Asian, 12 fixings", mcPrice(m, model, PayoffGeometricAsianCall(K, monthly), nPaths, batch, 42),
          geometricAsianCallExact(m, monthly, K));

    // One fixing at T: the Asian call is a European call
    check("geometric Asian, 1 fixing  ", mcPrice(m, model, PayoffGeometricAsianCall(K, {1.0}), nPaths, batch, 42),
          geometricAsianCallExact(m, {1.0}, K));

    // The arithmetic Asian call has no closed form: this is what Monte Carlo is for
    const McResult asian = mcPrice(m, model, PayoffAsianCall(K, monthly), nPaths, batch, 42);
    std::cout << "arithmetic Asian, 12 fixings: MC = " << asian.price << " +/- " << asian.stdErr << "\n";

    return failures;
}
