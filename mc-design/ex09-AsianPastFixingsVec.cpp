// ex09 - ASIAN options with PAST FIXINGS.
//
// An Asian option is observed at fixing dates. Valuing it in the middle of its life, some fixings
// are already in the PAST: they are known numbers, not random variables. They are not simulated
// and the model knows nothing about them. The PAYOFF object takes them in its constructor:
//
//     average = ( sum of the past fixings  +  sum of the simulated spots at the future dates ) / N
//
// where N = number of past + number of future fixings. The payoff precomputes the scalar "sum of
// the past fixings" once, and combines it with the simulated vector (one value per path) when it is
// evaluated. dates() contains only the FUTURE fixing dates: they are the simulation timeline.
//
// Everything else is as in ex08 (dictionary of models, views, correlation list), with the exact
// Black-Scholes model of ex07.
//
// Validation: the GEOMETRIC Asian call with past fixings has a closed form: log G is still normal,
// with the log of the past fixings shifting the mean. With no past fixing it is the ex05 formula.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct Market { double S0, r, q, sigma; };
struct McResult { double price, stdErr; };

// ---- Matrix views: a window on a block of an existing matrix (no copy, no ownership) ----
// Like std::span, but 2-D. Element (i,j) is in data[i*stride + j]; `stride` is the number of
// columns of the ORIGINAL matrix, so a block of rows/columns is still addressed correctly.
struct ConstMatrixView {
    ConstMatrixView(const double* data, size_t rows, size_t cols, size_t stride)
        : data(data), rows(rows), cols(cols), stride(stride) {}
    double operator()(size_t i, size_t j) const { return data[i * stride + j]; }
    std::span<const double> getRow(size_t i) const { return {data + i * stride, cols}; }
    // sub-block of nr rows and nc columns, starting at (r0, c0)
    ConstMatrixView block(size_t r0, size_t nr, size_t c0, size_t nc) const
        { return {data + r0 * stride + c0, nr, nc, stride}; }
    const double* data;
    size_t rows, cols, stride;
};

struct MatrixView {
    MatrixView(double* data, size_t rows, size_t cols, size_t stride)
        : data(data), rows(rows), cols(cols), stride(stride) {}
    double& operator()(size_t i, size_t j) const { return data[i * stride + j]; }
    std::span<double> getRow(size_t i) const { return {data + i * stride, cols}; }
    MatrixView block(size_t r0, size_t nr, size_t c0, size_t nc) const
        { return {data + r0 * stride + c0, nr, nc, stride}; }
    operator ConstMatrixView() const { return {data, rows, cols, stride}; }   // writable -> read-only
    double* data;
    size_t rows, cols, stride;
};

// Row-major matrix that OWNS its data: element (i,j) is in data[i*cols + j]
struct Matrix {
    Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), data(rows * cols) {}
    double&              operator()(size_t i, size_t j)       { return data[i * cols + j]; }
    double               operator()(size_t i, size_t j) const { return data[i * cols + j]; }
    std::span<double>       getRow(size_t i)       { return {data.data() + i * cols, cols}; }
    std::span<const double> getRow(size_t i) const { return {data.data() + i * cols, cols}; }
    MatrixView      view()       { return {data.data(), rows, cols, cols}; }
    ConstMatrixView view() const { return {data.data(), rows, cols, cols}; }
    // view of the nr rows starting at r0 (all the columns)
    MatrixView      rowsView(size_t r0, size_t nr)       { return view().block(r0, nr, 0, cols); }
    ConstMatrixView rowsView(size_t r0, size_t nr) const { return view().block(r0, nr, 0, cols); }
    size_t rows, cols;
    std::vector<double> data;
};

// States at the payoff dates: states[j] holds (state variables x paths) at the j-th payoff date.
// The state variables of ALL the models are stacked in the rows.
using States = std::vector<Matrix>;
// What ONE model sees of the states: ModelStates[j] is the view of its own rows at the date j
using ModelStates = std::vector<ConstMatrixView>;

// Correlation between two Brownians. Between models, a Brownian is referred to by "model.index":
// the model name and the position of the Brownian inside the model, e.g. "AAPL.0". The index is
// always required. Inside a model, the names are just the indices: "0", "1".
struct Correlation { std::string first, second; double rho; };


struct IModelMC {
    virtual ~IModelMC() = default;

    virtual size_t numberOfBrownians() const = 0;
    virtual size_t numberOfStateVariables() const = 0;

    // Correlations between the model's OWN Brownians, named by their index inside the model:
    // {"0", "1", rho}. Same format as the correlations between models, but without the model name
    // (the model does not know it): the engine adds it. All the pairs must be given.
    virtual std::vector<Correlation> intraModelCorrelation() const = 0;

    virtual void mergeDates(std::vector<double>& payoffDates) const = 0;
    virtual void configureTimeline(const std::vector<double>& grid,
                                   const std::vector<double>& payoffDates) = 0;

    // The matrices below are views of the model's OWN rows: row 0 is its first state variable
    // (or Brownian), whatever the position of the model in the engine.
    virtual void initializeStates(MatrixView states, size_t n) const = 0;
    virtual void generateStates(size_t k, ConstMatrixView rnd, ConstMatrixView prev,
                                MatrixView next, size_t n) const = 0;
    virtual void spot(size_t j, const ModelStates& states, size_t n, std::span<double> out) const = 0;
};

using Models = std::map<std::string, std::unique_ptr<IModelMC>>;   // name -> model

// Black-Scholes, exact on any grid. State variable: the Brownian motion W (one row).
// Written as if it were alone: it knows nothing of the other models.
struct BSExact : IModelMC {
    explicit BSExact(const Market& m) : m_spot0(m.S0), m_mu(m.r - m.q), m_sigma(m.sigma) {}

    size_t numberOfBrownians() const override      { return 1; }
    size_t numberOfStateVariables() const override { return 1; }

    std::vector<Correlation> intraModelCorrelation() const override { return {}; }   // one Brownian: none

    void mergeDates(std::vector<double>&) const override {}   // the payoff dates are enough

    void configureTimeline(const std::vector<double>& grid,
                           const std::vector<double>& payoffDates) override {
        m_sqrtDt.resize(grid.size());
        for (size_t k = 0; k < grid.size(); ++k)
            m_sqrtDt[k] = std::sqrt(grid[k] - ((k == 0) ? 0.0 : grid[k - 1]));
        m_drift.resize(payoffDates.size());              // needed at the payoff dates only
        for (size_t j = 0; j < payoffDates.size(); ++j)
            m_drift[j] = (m_mu - 0.5 * m_sigma * m_sigma) * payoffDates[j];
    }

    void initializeStates(MatrixView states, size_t n) const override {
        std::span<double> w = states.getRow(0);             // span on my only state variable
        for (size_t i = 0; i < n; ++i)
            w[i] = 0.0;                                  // W(0) = 0
    }

    void generateStates(size_t k, ConstMatrixView rnd, ConstMatrixView prev,
                        MatrixView next, size_t n) const override {
        std::span<const double> z  = rnd.getRow(0);         // spans on the relevant rows
        std::span<const double> w0 = prev.getRow(0);
        std::span<double>       w1 = next.getRow(0);
        for (size_t i = 0; i < n; ++i)
            w1[i] = w0[i] + m_sqrtDt[k] * z[i];
    }

    void spot(size_t j, const ModelStates& states, size_t n, std::span<double> out) const override {
        std::span<const double> w = states[j].getRow(0);
        for (size_t i = 0; i < n; ++i)
            out[i] = m_spot0 * std::exp(m_drift[j] + m_sigma * w[i]);
    }

    double m_spot0, m_mu, m_sigma;
    std::vector<double> m_sqrtDt;                        // per grid step
    std::vector<double> m_drift;                         // per payoff date
};

struct IPayoff {
    virtual ~IPayoff() = default;
    virtual const std::vector<double>& dates() const = 0;
    // The engine tells the payoff the position of every model (by name) in the array of models
    virtual void bind(const std::map<std::string, size_t>& modelIndex) = 0;
    // The payoff asks the models it needs by INDEX: models[i]->spot(j, states[i], ...), where
    // states[i] is what the model i sees of the states
    virtual void payoff(const std::vector<IModelMC*>& models, const std::vector<ModelStates>& states,
                        size_t n, std::span<double> out) const = 0;
};

// "AAPL.1" -> {"AAPL", 1}
std::pair<std::string, size_t> parseBrownian(const std::string& name)
{
    const size_t dot = name.rfind('.');
    if (dot == std::string::npos)
        throw std::invalid_argument("Brownian '" + name + "': expected 'model.index', e.g. 'AAPL.0'");
    return {name.substr(0, dot), static_cast<size_t>(std::stoul(name.substr(dot + 1)))};
}

// Lower triangular L with L L^T = R (R symmetric positive definite)
std::vector<std::vector<double>> cholesky(const std::vector<std::vector<double>>& R)
{
    const size_t n = R.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j <= i; ++j) {
            double x = R[i][j];
            for (size_t k = 0; k < j; ++k) x -= L[i][k] * L[j][k];
            if (i == j && x <= 0.0)
                throw std::invalid_argument("the correlation matrix is not positive definite");
            L[i][j] = (i == j) ? std::sqrt(x) : x / L[j][j];
        }
    return L;
}

// Position of every payoff date in the (sorted) grid: grid[obs[j]] is the j-th payoff date
std::vector<size_t> positions(const std::vector<double>& grid, const std::vector<double>& payoffDates)
{
    std::vector<size_t> obs;
    for (double t : payoffDates)
        obs.push_back(static_cast<size_t>(std::lower_bound(grid.begin(), grid.end(), t - 1e-12) - grid.begin()));
    return obs;
}

McResult mcPrice(double r, Models& modelDict, const std::vector<Correlation>& correlations,
                 IPayoff& payoff, size_t nPaths, size_t batchSize, unsigned seed)
{
    // ---- setup: dictionary -> array, and the coordinate maps ----
    std::vector<IModelMC*> models;                       // the array
    std::map<std::string, size_t> modelIndex;            // name -> position in the array
    std::vector<std::string> names;
    std::vector<size_t> firstState, firstBrownian;
    size_t nS = 0, nB = 0;
    for (auto& [name, model] : modelDict) {
        modelIndex[name] = models.size();
        names.push_back(name);
        models.push_back(model.get());
        firstState.push_back(nS);       nS += model->numberOfStateVariables();
        firstBrownian.push_back(nB);    nB += model->numberOfBrownians();
    }
    payoff.bind(modelIndex);                             // the payoff learns the model indices

    // ALL the correlations, in one list of the same format: the user's ones (between models) ...
    std::vector<Correlation> all = correlations;
    for (size_t m = 0; m < models.size(); ++m)           // ... and the models' ones (inside a model),
        for (const Correlation& c : models[m]->intraModelCorrelation())   // qualified with the name
            all.push_back({names[m] + "." + c.first, names[m] + "." + c.second, c.rho});

    // Correlation matrix of ALL the Brownians: 1 on the diagonal, the rest from the list
    std::vector<std::vector<double>> R(nB, std::vector<double>(nB, 0.0));
    std::vector<std::vector<bool>> given(nB, std::vector<bool>(nB, false));
    for (size_t b = 0; b < nB; ++b) { R[b][b] = 1.0; given[b][b] = true; }
    for (const Correlation& c : all) {
        const auto [nameA, idxA] = parseBrownian(c.first);     // "AAPL.0" -> model, Brownian
        const auto [nameB, idxB] = parseBrownian(c.second);
        const size_t a = firstBrownian[modelIndex.at(nameA)] + idxA;   // row of R
        const size_t b = firstBrownian[modelIndex.at(nameB)] + idxB;
        R[a][b] = R[b][a] = c.rho;
        given[a][b] = given[b][a] = true;
    }
    for (size_t a = 0; a < nB; ++a)                      // no defaults: a missing one is an error
        for (size_t b = a + 1; b < nB; ++b)
            if (!given[a][b])
                throw std::invalid_argument("missing correlation between Brownians " + std::to_string(a)
                                            + " and " + std::to_string(b));
    const auto L = cholesky(R);                          // R = L L^T, once

    // The simulation grid: each model may add dates
    std::vector<double> grid = payoff.dates();
    for (IModelMC* model : models) model->mergeDates(grid);
    const size_t nSteps = grid.size(), nDates = payoff.dates().size();
    const double T      = grid.back();
    for (IModelMC* model : models)
        model->configureTimeline(grid, payoff.dates());
    const std::vector<size_t> obs = positions(grid, payoff.dates());

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Work space allocated once, re-used by every batch
    Matrix z(nB, batchSize), rnd(nB, batchSize);         // independent / correlated normals
    Matrix bufferA(nS, batchSize), bufferB(nS, batchSize);
    Matrix *prev = &bufferA;
    States states(nDates, Matrix(nS, batchSize));
    // What each model sees: views of its own rows (the matrices are allocated once, so the views
    // stay valid for the whole simulation)
    std::vector<ModelStates> modelStates(models.size());
    for (size_t m = 0; m < models.size(); ++m)
        for (size_t j = 0; j < nDates; ++j)
            modelStates[m].push_back(states[j].rowsView(firstState[m], models[m]->numberOfStateVariables()));
    std::vector<double> pay(batchSize);

    double sum = 0.0, sumSq = 0.0;
    for (size_t done = 0; done < nPaths; done += batchSize) {
        const size_t n = std::min(batchSize, nPaths - done);

        prev = &bufferA;
        for (size_t m = 0; m < models.size(); ++m)
            models[m]->initializeStates(prev->rowsView(firstState[m], models[m]->numberOfStateVariables()), n);
        size_t j = 0;                                    // next payoff date
        for (size_t k = 0; k < nSteps; ++k) {
            for (size_t b = 0; b < nB; ++b)              // independent normals
                for (size_t i = 0; i < n; ++i)
                    z(b, i) = normal(rng);
            for (size_t a = 0; a < nB; ++a) {            // correlate: Zc_a = sum_{b<=a} L_ab Z_b
                for (size_t i = 0; i < n; ++i)
                    rnd(a, i) = 0.0;
                for (size_t b = 0; b <= a; ++b)
                    for (size_t i = 0; i < n; ++i)
                        rnd(a, i) += L[a][b] * z(b, i);
            }
            Matrix* spare = (prev == &bufferA) ? &bufferB : &bufferA;
            Matrix* next  = (j < nDates && k == obs[j]) ? &states[j++] : spare;
            for (size_t m = 0; m < models.size(); ++m) { // every model fills its own rows
                const size_t ns = models[m]->numberOfStateVariables();
                models[m]->generateStates(k, rnd.rowsView(firstBrownian[m], models[m]->numberOfBrownians()),
                                          prev->rowsView(firstState[m], ns),
                                          next->rowsView(firstState[m], ns), n);
            }
            prev = next;
        }

        payoff.payoff(models, modelStates, n, std::span<double>(pay).first(n));
        for (size_t i = 0; i < n; ++i) {
            sum   += pay[i];
            sumSq += pay[i] * pay[i];
        }
    }

    const double N    = static_cast<double>(nPaths);
    const double mean = sum / N;
    const double var  = (sumSq / N - mean * mean) * N / (N - 1.0);
    const double df   = std::exp(-r * T);
    return {df * mean, df * std::sqrt(var / N)};
}

// ---- payoffs ----
// Arithmetic Asian call with past fixings: max( average - K, 0 )
//   average = ( sum(past) + sum of the spots at the future dates ) / ( nPast + nFuture )
struct PayoffAsianCall : IPayoff {
    PayoffAsianCall(double strike, std::vector<double> futureDates, const std::vector<double>& pastFixings,
                    std::string model)
        : m_strike(strike), m_dates(std::move(futureDates)), m_name(std::move(model)),
          m_count(pastFixings.size() + m_dates.size()) {
        for (double f : pastFixings) m_pastSum += f;     // the scalar of the past, computed once
    }

    const std::vector<double>& dates() const override { return m_dates; }

    void bind(const std::map<std::string, size_t>& modelIndex) override { m_model = modelIndex.at(m_name); }

    void payoff(const std::vector<IModelMC*>& models, const std::vector<ModelStates>& states,
                size_t n, std::span<double> out) const override {
        std::vector<double> spot(n);
        std::fill(out.begin(), out.end(), m_pastSum);    // start from the past: a scalar, for all paths
        for (size_t j = 0; j < m_dates.size(); ++j) {    // add the simulated spots, one date at a time
            models[m_model]->spot(j, states[m_model], n, spot);
            for (size_t i = 0; i < n; ++i)
                out[i] += spot[i];
        }
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(out[i] / static_cast<double>(m_count) - m_strike, 0.0);
    }

    double m_strike;
    std::vector<double> m_dates;
    std::string m_name;
    size_t m_count, m_model = 0;
    double m_pastSum = 0.0;
};

// Geometric Asian call with past fixings: max( (prod of all the fixings)^(1/N) - K, 0 )
struct PayoffGeometricAsianCall : IPayoff {
    PayoffGeometricAsianCall(double strike, std::vector<double> futureDates, const std::vector<double>& pastFixings,
                             std::string model)
        : m_strike(strike), m_dates(std::move(futureDates)), m_name(std::move(model)),
          m_count(pastFixings.size() + m_dates.size()) {
        for (double f : pastFixings) m_pastLogSum += std::log(f);   // the scalar of the past
    }

    const std::vector<double>& dates() const override { return m_dates; }

    void bind(const std::map<std::string, size_t>& modelIndex) override { m_model = modelIndex.at(m_name); }

    void payoff(const std::vector<IModelMC*>& models, const std::vector<ModelStates>& states,
                size_t n, std::span<double> out) const override {
        std::vector<double> spot(n);
        std::fill(out.begin(), out.end(), m_pastLogSum);
        for (size_t j = 0; j < m_dates.size(); ++j) {
            models[m_model]->spot(j, states[m_model], n, spot);
            for (size_t i = 0; i < n; ++i)
                out[i] += std::log(spot[i]);
        }
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(std::exp(out[i] / static_cast<double>(m_count)) - m_strike, 0.0);
    }

    double m_strike;
    std::vector<double> m_dates;
    std::string m_name;
    size_t m_count, m_model = 0;
    double m_pastLogSum = 0.0;
};

// ---- closed form ----
double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// E[max(X-K,0)] for X lognormal with E[log X] = M and Var[log X] = V, discounted
double lognormalCall(double M, double V, double K, double df)
{
    const double sd = std::sqrt(V);
    const double d1 = (M + V - std::log(K)) / sd;
    return df * (std::exp(M + 0.5 * V) * normCdf(d1) - K * normCdf(d1 - sd));
}

// Geometric Asian with past fixings. t: the FUTURE fixing times (years from today).
//   log G = ( sum log(past) + sum_j log S(t_j) ) / N,   N = nPast + nFuture
double geometricAsianCallExact(const Market& m, const std::vector<double>& t,
                               const std::vector<double>& past, double K)
{
    const size_t nF = t.size(), N = nF + past.size();
    double sumT = 0.0, var = 0.0, sumLogPast = 0.0;
    for (size_t i = 0; i < nF; ++i) {
        sumT += t[i];
        for (size_t j = 0; j < nF; ++j)
            var += std::min(t[i], t[j]);                 // Cov(W(t_i), W(t_j)) = min(t_i, t_j)
    }
    for (double f : past) sumLogPast += std::log(f);
    const double NN = static_cast<double>(N);
    const double M  = (sumLogPast + nF * std::log(m.S0) + (m.r - m.q - 0.5 * m.sigma * m.sigma) * sumT) / NN;
    const double V  = m.sigma * m.sigma * var / (NN * NN);
    return lognormalCall(M, V, K, std::exp(-m.r * t.back()));
}

int main()
{
    const Market m{100.0, 0.03, 0.01, 0.20};
    const double K = 100.0;
    const size_t nPaths = 1'000'000, batch = 4096;
    int failures = 0;

    // A one-year Asian with 12 monthly fixings, valued after 4 months: 4 fixings are in the past
    const std::vector<double> past{98.0, 101.0, 103.0, 99.0};
    std::vector<double> future;                          // the 8 remaining fixings, from today
    for (int j = 1; j <= 8; ++j) future.push_back(j / 12.0);

    auto run = [&](const char* what, const McResult& mc, double exact) {
        const double z = (mc.price - exact) / mc.stdErr;
        std::cout << what << ": MC = " << mc.price << " +/- " << mc.stdErr
                  << "  exact = " << exact << "  (" << z << " std errors)\n";
        if (std::abs(z) > 4.0) ++failures;
    };

    auto price = [&](IPayoff& payoff) {
        Models models;                                   // the dictionary: one model
        models["AAPL"] = std::make_unique<BSExact>(m);
        return mcPrice(m.r, models, {}, payoff, nPaths, batch, 42);   // one model: no correlation
    };

    PayoffGeometricAsianCall geoPast(K, future, past, "AAPL");
    run("geometric Asian, 4 past + 8 future fixings", price(geoPast),
        geometricAsianCallExact(m, future, past, K));

    // No past fixing: back to the Asian of ex05
    PayoffGeometricAsianCall geoNone(K, future, {}, "AAPL");
    run("geometric Asian, 0 past + 8 future fixings", price(geoNone),
        geometricAsianCallExact(m, future, {}, K));

    // The arithmetic Asian has no closed form: this is what Monte Carlo is for
    PayoffAsianCall asianPast(K, future, past, "AAPL");
    const McResult asian = price(asianPast);
    std::cout << "arithmetic Asian, 4 past + 8 future fixings: MC = " << asian.price
              << " +/- " << asian.stdErr << "\n";

    return failures;
}
