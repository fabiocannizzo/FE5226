// ex06 - Asian call, EULER discretization: the model is no longer exact.
//
// With a model that cannot be simulated exactly (local volatility, stochastic volatility, ...)
// the paths are built by small time steps, and the result carries a DISCRETIZATION BIAS that
// shrinks with the step size. One new idea, visible in the code:
//   - the model decides the simulation grid: IModelMC gets a method mergeDates, which turns
//     the payoff dates into the (finer) simulation dates, IN PLACE.
//   - the model is then configured with the grid AND the payoff dates (spot(j, ...) still
//     refers to the j-th payoff date).
//   - the engine keeps the states only at the PAYOFF dates (the cube has one matrix per payoff
//     date): a payoff date is generated straight into the cube, the other steps alternate
//     between two auxiliary matrices.
// The payoff is untouched: it is the one of ex05, and works with the exact and the Euler model.
//
// Euler-Maruyama (state variable: the spot itself):
//     S_{k+1} = S_k (1 + (r-q) dt + sigma sqrt(dt) Z_k), truncated at 0.
//
// The bias is measured against the closed form of the geometric Asian call.
#include <algorithm>
#include <cmath>
#include <iomanip>
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

// States at the payoff dates: states[j] holds (state variables x paths) at the j-th payoff date
using States = std::vector<Matrix>;

struct IModelMC {
    virtual ~IModelMC() = default;

    virtual size_t numberOfBrownians() const = 0;
    virtual size_t numberOfStateVariables() const = 0;

    // In: the payoff dates (years from today, increasing, all > 0).
    // Out: the dates of the simulation, which must include the payoff dates.
    virtual void mergeDates(std::vector<double>& payoffDates) const = 0;

    // The engine tells the model the simulation grid (as returned by mergeDates) and the payoff
    // dates. The model precomputes whatever it needs.
    virtual void configureTimeline(const std::vector<double>& grid,
                                   const std::vector<double>& payoffDates) = 0;

    // states at today, for the first n paths
    virtual void initializeStates(Matrix& states, size_t n) const = 0;

    // One step, from the grid date k-1 (today if k == 0) to the grid date k, for the first n paths
    virtual void generateStates(size_t k, const Matrix& rnd, const Matrix& prev,
                                Matrix& next, size_t n) const = 0;

    // Spot at the PAYOFF date j, for the first n paths: out[i] = S(payoffDates[j]) on path i,
    // rebuilt from states[j]
    virtual void spot(size_t j, const States& states, size_t n, std::span<double> out) const = 0;
};

// Black-Scholes, exact on any grid. State variable: the Brownian motion W.
struct BSExact : IModelMC {
    explicit BSExact(const Market& m) : m_spot0(m.S0), m_mu(m.r - m.q), m_sigma(m.sigma) {}

    size_t numberOfBrownians() const override      { return 1; }
    size_t numberOfStateVariables() const override { return 1; }

    void mergeDates(std::vector<double>&) const override {}   // the payoff dates are enough

    void configureTimeline(const std::vector<double>& grid,
                           const std::vector<double>& payoffDates) override {
        m_sqrtDt.resize(grid.size());
        for (size_t k = 0; k < grid.size(); ++k) {
            const double tPrev = (k == 0) ? 0.0 : grid[k - 1];
            m_sqrtDt[k] = std::sqrt(grid[k] - tPrev);
        }
        m_drift.resize(payoffDates.size());              // needed at the payoff dates only
        for (size_t j = 0; j < payoffDates.size(); ++j)
            m_drift[j] = (m_mu - 0.5 * m_sigma * m_sigma) * payoffDates[j];
    }

    void initializeStates(Matrix& states, size_t n) const override {
        for (size_t i = 0; i < n; ++i)
            states(0, i) = 0.0;                          // W(0) = 0
    }

    void generateStates(size_t k, const Matrix& rnd, const Matrix& prev,
                        Matrix& next, size_t n) const override {
        for (size_t i = 0; i < n; ++i)
            next(0, i) = prev(0, i) + m_sqrtDt[k] * rnd(0, i);
    }

    void spot(size_t j, const States& states, size_t n, std::span<double> out) const override {
        for (size_t i = 0; i < n; ++i)
            out[i] = m_spot0 * std::exp(m_drift[j] + m_sigma * states[j](0, i));
    }

    double m_spot0, m_mu, m_sigma;
    std::vector<double> m_sqrtDt;                        // per grid step
    std::vector<double> m_drift;                         // per payoff date
};

// Black-Scholes, Euler scheme. State variable: the spot S itself.
struct BSEuler : IModelMC {
    BSEuler(const Market& m, double maxStep)
        : m_spot0(m.S0), m_mu(m.r - m.q), m_sigma(m.sigma), m_maxStep(maxStep) {}

    size_t numberOfBrownians() const override      { return 1; }
    size_t numberOfStateVariables() const override { return 1; }

    // Splits every interval between payoff dates in equal steps no longer than maxStep
    void mergeDates(std::vector<double>& dates) const override {
        std::vector<double> grid;
        double prev = 0.0;
        for (double t : dates) {
            const int nSub = static_cast<int>(std::ceil((t - prev) / m_maxStep - 1e-12));
            for (int k = 1; k < nSub; ++k)
                grid.push_back(prev + (t - prev) * k / nSub);
            grid.push_back(t);                           // the payoff date itself, exactly
            prev = t;
        }
        dates = std::move(grid);
    }

    void configureTimeline(const std::vector<double>& grid,
                           const std::vector<double>&) override {
        m_a.resize(grid.size());
        m_b.resize(grid.size());
        for (size_t k = 0; k < grid.size(); ++k) {
            const double dt = grid[k] - ((k == 0) ? 0.0 : grid[k - 1]);
            m_a[k] = 1.0 + m_mu * dt;
            m_b[k] = m_sigma * std::sqrt(dt);
        }
    }

    void initializeStates(Matrix& states, size_t n) const override {
        for (size_t i = 0; i < n; ++i)
            states(0, i) = m_spot0;                      // S(0) = S0
    }

    void generateStates(size_t k, const Matrix& rnd, const Matrix& prev,
                        Matrix& next, size_t n) const override {
        for (size_t i = 0; i < n; ++i)
            next(0, i) = std::max(prev(0, i) * (m_a[k] + m_b[k] * rnd(0, i)), 0.0);
    }

    void spot(size_t j, const States& states, size_t n, std::span<double> out) const override {
        for (size_t i = 0; i < n; ++i)
            out[i] = states[j](0, i);                    // the state is the spot
    }

    double m_spot0, m_mu, m_sigma, m_maxStep;
    std::vector<double> m_a, m_b;
};

struct IPayoff {
    virtual ~IPayoff() = default;
    // Dates (years from today, increasing) at which the payoff needs the spot; paid at the last
    virtual const std::vector<double>& dates() const = 0;
    // out[i] = payoff on path i, for the first n paths. The payoff queries the model for the
    // spot at each of ITS dates, by index: model.spot(j, states, n, ...) is the j-th date.
    virtual void payoff(const IModelMC& model, const States& states, size_t n,
                        std::span<double> out) const = 0;
};

// Position of every payoff date in the (sorted) grid: grid[obs[j]] is the j-th payoff date
std::vector<size_t> positions(const std::vector<double>& grid, const std::vector<double>& payoffDates)
{
    std::vector<size_t> obs;
    for (double t : payoffDates)
        obs.push_back(static_cast<size_t>(std::lower_bound(grid.begin(), grid.end(), t - 1e-12) - grid.begin()));
    return obs;
}

McResult mcPrice(const Market& m, IModelMC& model, const IPayoff& payoff,
                 size_t nPaths, size_t batchSize, unsigned seed)
{
    std::vector<double> grid = payoff.dates();           // a copy ...
    model.mergeDates(grid);                              // ... which the model turns into the grid
    const size_t nSteps = grid.size(), nDates = payoff.dates().size();
    const double T      = grid.back();
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    model.configureTimeline(grid, payoff.dates());       // the model learns the grid
    const std::vector<size_t> obs = positions(grid, payoff.dates());

    // Work space allocated once, re-used by every batch
    const size_t nS = model.numberOfStateVariables(), nB = model.numberOfBrownians();
    Matrix rnd(nB, batchSize);
    Matrix bufferA(nS, batchSize), bufferB(nS, batchSize);   // two auxiliary matrices, alternating
    Matrix *prev = &bufferA;
    States states(nDates, Matrix(nS, batchSize));        // stored only at the payoff dates
    std::vector<double> pay(batchSize);

    double sum = 0.0, sumSq = 0.0;
    for (size_t done = 0; done < nPaths; done += batchSize) {
        const size_t n = std::min(batchSize, nPaths - done);   // the last batch may be shorter

        prev = &bufferA;
        model.initializeStates(*prev, n);
        size_t j = 0;                                           // next payoff date
        for (size_t k = 0; k < nSteps; ++k) {                   // advance all paths by one step
            for (size_t b = 0; b < nB; ++b)
                for (size_t i = 0; i < n; ++i)
                    rnd(b, i) = normal(rng);
            // A payoff date is written straight into the cube, any other step in the buffer
            // which does not hold the previous states (so the buffers alternate)
            Matrix* spare = (prev == &bufferA) ? &bufferB : &bufferA;
            Matrix* next  = (j < nDates && k == obs[j]) ? &states[j++] : spare;
            model.generateStates(k, rnd, *prev, *next, n);
            prev = next;
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
    const Market m{100.0, 0.03, 0.01, 0.40};            // high vol: the bias is easier to see
    const double K = 100.0;
    const size_t nPaths = 2'000'000, batch = 4096;

    std::vector<double> quarterly{0.25, 0.5, 0.75, 1.0};
    const PayoffGeometricAsianCall geometric(K, quarterly);
    const double exact = geometricAsianCallExact(m, quarterly, K);
    std::cout << "geometric Asian call, 4 fixings, exact = " << exact << "\n\n"
              << std::fixed << std::setprecision(4);

    int failures = 0;
    {
        BSExact model(m);
        const McResult r = mcPrice(m, model, geometric, nPaths, batch, 42);
        const double z = (r.price - exact) / r.stdErr;
        std::cout << "exact model          : " << r.price << " +/- " << r.stdErr << "  bias = " << z << " std errors\n";
        if (std::abs(z) > 4.0) ++failures;
    }
    for (double maxStep : {0.25, 0.10, 0.05, 0.01}) {
        BSEuler model(m, maxStep);
        const McResult r = mcPrice(m, model, geometric, nPaths, batch, 42);
        const double z = (r.price - exact) / r.stdErr;
        std::cout << "Euler, maxStep = " << std::setw(4) << maxStep << " : " << r.price << " +/- " << r.stdErr
                  << "  bias = " << z << " std errors\n";
        if (maxStep <= 0.01 && std::abs(z) > 4.0) ++failures;
    }

    BSEuler model(m, 0.01);
    const McResult asian = mcPrice(m, model, PayoffAsianCall(K, quarterly), nPaths, batch, 42);
    std::cout << "\narithmetic Asian call (Euler, maxStep = 0.01): " << asian.price << " +/- " << asian.stdErr << "\n";
    return failures;
}
