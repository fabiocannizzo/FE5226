// ex10 - CURRENCY TRIANGLES.
//
// Three currencies, EUR, USD, GBP, make three exchange rates but only TWO are independent:
//     EURGBP = EURUSD / GBPUSD
// The engine simulates in ONE currency, USD: the two USD pairs EURUSD and GBPUSD are models
// (Black-Scholes; domestic rate = USD, foreign rate = EUR or GBP). The cross EURGBP is NOT
// simulated: it is a DUMMY MODEL, with zero Brownians and zero state variables, which takes no part
// in the diffusion. It DEPENDS on the models "EURUSD" and "GBPUSD" (dependencies(), by name); the
// engine gives them to it (bindDependencies). Its spot() asks the two models for their spots and
// divides.
//
// The payoff asks for "EURGBP" like for any other model: it does not know how it is obtained.
//
// Validation:
//   1. a call on EURGBP, strike K in GBP, paid in GBP at T, converted into USD at GBPUSD(T): the
//      cross is lognormal with volatility sX^2 = sE^2 + sG^2 - 2 rho sE sG, Garman-Kohlhagen
//      in GBP, converted at GBPUSD(0).
//   2. a call on the sum of the performances EURUSD(T)/EURUSD(0) + EURGBP(T)/EURGBP(0), paid in
//      USD: with K = 0 the price is a sum of two forwards (closed form).
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
    // Spot at the payoff date j. The states are the model's own (none for a model without states).
    virtual void spot(size_t j, const ModelStates& states, size_t n, std::span<double> out) const = 0;

    // A model may DEPEND on other models, by name (default: none). The engine resolves the names
    // and gives the model those models and their states (bindDependencies), once, before the simulation.
    virtual std::vector<std::string> dependencies() const { return {}; }
    virtual void bindDependencies(const std::vector<const IModelMC*>&, const std::vector<const ModelStates*>&) {}
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

// Triangle: a model WITHOUT Brownians and WITHOUT state variables. It reconstructs the spot of a
// cross from two other models:  cross = base1 / base2   (e.g. EURGBP = EURUSD / GBPUSD).
struct TriangleModel : IModelMC {
    TriangleModel(std::string numerator, std::string denominator)
        : m_numeratorName(std::move(numerator)), m_denominatorName(std::move(denominator)) {}

    size_t numberOfBrownians() const override      { return 0; }
    size_t numberOfStateVariables() const override { return 0; }
    std::vector<Correlation> intraModelCorrelation() const override { return {}; }

    void mergeDates(std::vector<double>&) const override {}
    void configureTimeline(const std::vector<double>&, const std::vector<double>&) override {}
    void initializeStates(MatrixView, size_t) const override {}                  // nothing to do
    void generateStates(size_t, ConstMatrixView, ConstMatrixView, MatrixView, size_t) const override {}

    std::vector<std::string> dependencies() const override { return {m_numeratorName, m_denominatorName}; }
    void bindDependencies(const std::vector<const IModelMC*>& models,
                          const std::vector<const ModelStates*>& states) override {
        m_numerator = models[0];   m_numeratorStates = states[0];
        m_denominator = models[1]; m_denominatorStates = states[1];
    }

    // the spots of the two models, and the division
    void spot(size_t j, const ModelStates&, size_t n, std::span<double> out) const override {
        std::vector<double> denominator(n);
        m_numerator->spot(j, *m_numeratorStates, n, out);
        m_denominator->spot(j, *m_denominatorStates, n, denominator);
        for (size_t i = 0; i < n; ++i)
            out[i] /= denominator[i];
    }

    std::string m_numeratorName, m_denominatorName;
    const IModelMC *m_numerator = nullptr, *m_denominator = nullptr;
    const ModelStates *m_numeratorStates = nullptr, *m_denominatorStates = nullptr;
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
    // Models that depend on other models (by name) learn them now. The views are all built,
    // so their addresses do not change any more.
    for (size_t m = 0; m < models.size(); ++m) {
        std::vector<const IModelMC*> deps;
        std::vector<const ModelStates*> depStates;
        for (const std::string& name : models[m]->dependencies()) {
            deps.push_back(models[modelIndex.at(name)]);
            depStates.push_back(&modelStates[modelIndex.at(name)]);
        }
        models[m]->bindDependencies(deps, depStates);
    }
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
// Call on a cross, paid in the currency of the denominator, converted into USD at the payment date:
//     fxToUsd(T) * max( cross(T) - K , 0 )         (e.g. GBPUSD(T) * max(EURGBP(T) - K, 0))
struct PayoffFxCall : IPayoff {
    PayoffFxCall(double strike, double T, std::string cross, std::string fxToUsd)
        : m_strike(strike), m_dates{T}, m_crossName(std::move(cross)), m_fxName(std::move(fxToUsd)) {}

    const std::vector<double>& dates() const override { return m_dates; }

    void bind(const std::map<std::string, size_t>& modelIndex) override {
        m_cross = modelIndex.at(m_crossName);
        m_fx    = modelIndex.at(m_fxName);
    }

    void payoff(const std::vector<IModelMC*>& models, const std::vector<ModelStates>& states,
                size_t n, std::span<double> out) const override {
        std::vector<double> fx(n);
        models[m_cross]->spot(0, states[m_cross], n, out);   // EURGBP(T): the payoff does not know HOW
        models[m_fx]->spot(0, states[m_fx], n, fx);          // GBPUSD(T)
        for (size_t i = 0; i < n; ++i)
            out[i] = fx[i] * std::max(out[i] - m_strike, 0.0);
    }

    double m_strike;
    std::vector<double> m_dates;
    std::string m_crossName, m_fxName;
    size_t m_cross = 0, m_fx = 0;
};

// Call on the sum of the performances of two models, paid in USD:
//     max( first(T)/first0 + second(T)/second0 - K , 0 )
struct PayoffFxPerformanceCall : IPayoff {
    PayoffFxPerformanceCall(double strike, double T, std::string first, double first0,
                            std::string second, double second0)
        : m_strike(strike), m_dates{T}, m_firstName(std::move(first)), m_first0(first0),
          m_secondName(std::move(second)), m_second0(second0) {}

    const std::vector<double>& dates() const override { return m_dates; }

    void bind(const std::map<std::string, size_t>& modelIndex) override {
        m_first  = modelIndex.at(m_firstName);
        m_second = modelIndex.at(m_secondName);
    }

    void payoff(const std::vector<IModelMC*>& models, const std::vector<ModelStates>& states,
                size_t n, std::span<double> out) const override {
        std::vector<double> s2(n);
        models[m_first]->spot(0, states[m_first], n, out);
        models[m_second]->spot(0, states[m_second], n, s2);
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(out[i] / m_first0 + s2[i] / m_second0 - m_strike, 0.0);
    }

    double m_strike;
    std::vector<double> m_dates;
    std::string m_firstName; double m_first0;
    std::string m_secondName; double m_second0;
    size_t m_first = 0, m_second = 0;
};

double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Garman-Kohlhagen: call on an exchange rate S (units of domestic currency per foreign one);
// rd, rf: domestic and foreign rates. Price in the domestic currency.
double garmanKohlhagen(double S, double K, double rd, double rf, double sigma, double T)
{
    const double sd = sigma * std::sqrt(T);
    const double d1 = (std::log(S / K) + (rd - rf) * T) / sd + 0.5 * sd;
    return S * std::exp(-rf * T) * normCdf(d1) - K * std::exp(-rd * T) * normCdf(d1 - sd);
}

int main()
{
    // Rates of the three currencies, and the USD pairs: EURUSD (USD per EUR), GBPUSD (USD per GBP)
    const double rU = 0.04, rE = 0.025, rG = 0.035;
    const double eurusd0 = 1.10, sE = 0.08, gbpusd0 = 1.30, sG = 0.09, rho = 0.6;   // rho(EURUSD, GBPUSD)
    const double T = 1.0;
    const size_t nPaths = 1'000'000, batch = 4096;
    int failures = 0;

    // The simulated models: the two USD pairs. The cross is a dummy model built from them.
    auto makeModels = [&]() {
        Models models;
        models["EURUSD"] = std::make_unique<BSExact>(Market{eurusd0, rU, rE, sE});
        models["GBPUSD"] = std::make_unique<BSExact>(Market{gbpusd0, rU, rG, sG});
        models["EURGBP"] = std::make_unique<TriangleModel>("EURUSD", "GBPUSD");   // no Brownian, no state
        return models;
    };
    // The user gives the correlation between the two USD pairs: EURGBP has no Brownian
    const std::vector<Correlation> correlations{{"EURUSD.0", "GBPUSD.0", rho}};

    const double eurgbp0 = eurusd0 / gbpusd0;
    const double sX      = std::sqrt(sE * sE + sG * sG - 2.0 * rho * sE * sG);   // volatility of EURGBP

    // ---- 1. call on EURGBP, strike in GBP: Garman-Kohlhagen in GBP, converted at GBPUSD(0)
    {
        const double K = 0.85;
        Models models = makeModels();
        PayoffFxCall payoff(K, T, "EURGBP", "GBPUSD");
        const McResult res = mcPrice(rU, models, correlations, payoff, nPaths, batch, 42);
        const double exact = gbpusd0 * garmanKohlhagen(eurgbp0, K, rG, rE, sX, T);
        const double z = (res.price - exact) / res.stdErr;
        std::cout << "EURGBP call, K = " << K << " GBP, in USD: MC = " << res.price << " +/- " << res.stdErr
                  << "  Garman-Kohlhagen = " << exact << "  (" << z << " std errors)\n";
        if (std::abs(z) > 4.0) ++failures;
    }

    // ---- 2. call on the sum of the performances of EURUSD and EURGBP, paid in USD
    {
        // K = 0: the sum of two forwards. E[EURGBP(T)] / EURGBP(0) = exp( (rG - rE + sG^2 - rho sE sG) T )
        Models models = makeModels();
        PayoffFxPerformanceCall payoff0(0.0, T, "EURUSD", eurusd0, "EURGBP", eurgbp0);
        const McResult res0 = mcPrice(rU, models, correlations, payoff0, nPaths, batch, 42);
        const double exact0 = std::exp(-rU * T) * (std::exp((rU - rE) * T) + std::exp((rG - rE + sG * sG - rho * sE * sG) * T));
        const double z0 = (res0.price - exact0) / res0.stdErr;
        std::cout << "\nperformances, K = 0: MC = " << res0.price << " +/- " << res0.stdErr
                  << "  exact = " << exact0 << "  (" << z0 << " std errors)\n";
        if (std::abs(z0) > 4.0) ++failures;

        PayoffFxPerformanceCall payoff2(2.0, T, "EURUSD", eurusd0, "EURGBP", eurgbp0);
        const McResult res2 = mcPrice(rU, models, correlations, payoff2, nPaths, batch, 42);
        std::cout << "performances, K = 2: MC = " << res2.price << " +/- " << res2.stdErr << "\n";
    }
    return failures;
}
