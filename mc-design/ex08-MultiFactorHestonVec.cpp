// ex08 - MULTIFACTOR models: AAPL and MSFT both follow a HESTON model.
//
//     dS = (r-q) S dt + sqrt(v) S dW1                 spot
//     dv = kappa (theta - v) dt + xi sqrt(v) dW2      variance,    dW1 dW2 = rho dt
//
// Incremental to ex06 (Euler steps, mergeDates) and ex07 (a dictionary of models, matrix views):
//   - a model with TWO risk factors: numberOfBrownians() = numberOfStateVariables() = 2
//       state variables : x = ln S  and  v          (rows 0 and 1 of the model's states)
//       Brownians       : Z1 (spot) and Z2 (variance)
//   - the CORRELATION MATRIX of all the Brownians is now 4 x 4 (AAPL.0, AAPL.1, MSFT.0, MSFT.1):
//       * INSIDE a model (AAPL.0 vs AAPL.1) : provided by the MODEL, intraModelCorrelation(), as a
//         list of Correlation, with the same format as the user's ones but names without the model
//         ("0", "1"). The engine adds the model name: "0","1" of AAPL -> "AAPL.0","AAPL.1".
//       * BETWEEN models (AAPL.x vs MSFT.y) : provided by the user, by name.
//     All together they must cover every pair of Brownians: a missing correlation is an error.
//
// Euler scheme with FULL TRUNCATION (v+ = max(v,0) wherever v is used), on x = ln S:
//     x' = x + (r - q - v+/2) dt + sqrt(v+ dt) Z1
//     v' = v + kappa (theta - v+) dt + xi sqrt(v+ dt) Z2
//
// Validation:
//   1. xi = 0, v0 = theta: the variance is constant and Heston is Black-Scholes with sigma = sqrt(v0):
//      the exchange option has Margrabe's closed form (this tests the cross-model correlation).
//   2. the European call on AAPL has a semi-closed form (Heston characteristic function), which
//      tests the dynamics and the correlation inside the model.
#include <algorithm>
#include <cmath>
#include <complex>
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

struct HestonParams { double S0, r, q, v0, kappa, theta, xi, rho; };   // rho: spot / variance

// Heston, Euler with full truncation. State variables: x = ln S (row 0) and v (row 1).
// Brownians: Z1 (row 0), Z2 (row 1). Written as if alone: it knows nothing of the other models.
struct HestonEuler : IModelMC {
    HestonEuler(const HestonParams& p, double maxStep) : m_p(p), m_maxStep(maxStep) {}

    size_t numberOfBrownians() const override      { return 2; }
    size_t numberOfStateVariables() const override { return 2; }

    // the model's own Brownians "0" (spot) and "1" (variance)
    std::vector<Correlation> intraModelCorrelation() const override { return {{"0", "1", m_p.rho}}; }

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

    void configureTimeline(const std::vector<double>& grid, const std::vector<double>&) override {
        m_dt.resize(grid.size());
        m_sqrtDt.resize(grid.size());
        for (size_t k = 0; k < grid.size(); ++k) {
            m_dt[k]     = grid[k] - ((k == 0) ? 0.0 : grid[k - 1]);
            m_sqrtDt[k] = std::sqrt(m_dt[k]);
        }
    }

    void initializeStates(MatrixView states, size_t n) const override {
        std::span<double> x = states.getRow(0), v = states.getRow(1);
        for (size_t i = 0; i < n; ++i) {
            x[i] = std::log(m_p.S0);                     // x(0) = ln S0
            v[i] = m_p.v0;                               // v(0) = v0
        }
    }

    void generateStates(size_t k, ConstMatrixView rnd, ConstMatrixView prev,
                        MatrixView next, size_t n) const override {
        std::span<const double> z1 = rnd.getRow(0),  z2 = rnd.getRow(1);    // the Brownians
        std::span<const double> x0 = prev.getRow(0), v0 = prev.getRow(1);   // states before
        std::span<double>       x1 = next.getRow(0), v1 = next.getRow(1);   // states after
        const double dt = m_dt[k], sqrtDt = m_sqrtDt[k], mu = m_p.r - m_p.q;
        for (size_t i = 0; i < n; ++i) {
            const double vp = std::max(v0[i], 0.0);      // full truncation
            const double sv = std::sqrt(vp) * sqrtDt;
            x1[i] = x0[i] + (mu - 0.5 * vp) * dt + sv * z1[i];
            v1[i] = v0[i] + m_p.kappa * (m_p.theta - vp) * dt + m_p.xi * sv * z2[i];
        }
    }

    void spot(size_t j, const ModelStates& states, size_t n, std::span<double> out) const override {
        std::span<const double> x = states[j].getRow(0);
        for (size_t i = 0; i < n; ++i)
            out[i] = std::exp(x[i]);
    }

    HestonParams m_p;
    double m_maxStep;
    std::vector<double> m_dt, m_sqrtDt;                  // per grid step
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

// European call on ONE model, given by name: max( S(T) - K , 0 ), paid at T
struct PayoffCall : IPayoff {
    PayoffCall(double K, double T, std::string name) : m_strike(K), m_dates{T}, m_name(std::move(name)) {}

    const std::vector<double>& dates() const override { return m_dates; }

    void bind(const std::map<std::string, size_t>& modelIndex) override { m_model = modelIndex.at(m_name); }

    void payoff(const std::vector<IModelMC*>& models, const std::vector<ModelStates>& states,
                size_t n, std::span<double> out) const override {
        models[m_model]->spot(0, states[m_model], n, out);
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(out[i] - m_strike, 0.0);
    }

    double m_strike;
    std::vector<double> m_dates;
    std::string m_name;
    size_t m_model = 0;
};

// Exchange option: max( S1(T) - S2(T) , 0 ), paid at T. The two models are given BY NAME
struct PayoffExchange : IPayoff {
    PayoffExchange(double T, std::string first, std::string second)
        : m_dates{T}, m_firstName(std::move(first)), m_secondName(std::move(second)) {}

    const std::vector<double>& dates() const override { return m_dates; }

    void bind(const std::map<std::string, size_t>& modelIndex) override {
        m_first  = modelIndex.at(m_firstName);           // names -> indices, once
        m_second = modelIndex.at(m_secondName);
    }

    void payoff(const std::vector<IModelMC*>& models, const std::vector<ModelStates>& states,
                size_t n, std::span<double> out) const override {
        std::vector<double> s2(n);
        models[m_first]->spot(0, states[m_first], n, out);      // date 0 (= T): S1(T)
        models[m_second]->spot(0, states[m_second], n, s2);     // date 0 (= T): S2(T)
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(out[i] - s2[i], 0.0);
    }

    std::vector<double> m_dates;
    std::string m_firstName, m_secondName;
    size_t m_first = 0, m_second = 0;                    // positions in the array of models
};

double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Heston European call (semi-closed form): C = S e^{-qT} P1 - K e^{-rT} P2, with
//   Pj = 1/2 + 1/pi Int_0^inf Re[ e^{-iu ln K} phi_j(u) / (iu) ] du
// phi_j: characteristic function of ln S(T), in the "little trap" formulation (stable for long T)
double hestonCall(const HestonParams& p, double K, double T)
{
    using C = std::complex<double>;
    const C I(0.0, 1.0);
    auto P = [&](int j) {
        const double uj = (j == 1) ? 0.5 : -0.5;
        const double b  = (j == 1) ? p.kappa - p.rho * p.xi : p.kappa;
        const double a  = p.kappa * p.theta, xi2 = p.xi * p.xi;
        auto integrand = [&](double u) {
            const C iu = I * u;
            const C d  = std::sqrt((p.rho * p.xi * iu - b) * (p.rho * p.xi * iu - b) - xi2 * (2.0 * uj * iu - u * u));
            const C g  = (b - p.rho * p.xi * iu - d) / (b - p.rho * p.xi * iu + d);
            const C e  = std::exp(-d * T);
            const C Cj = (p.r - p.q) * iu * T + a / xi2 * ((b - p.rho * p.xi * iu - d) * T - 2.0 * std::log((1.0 - g * e) / (1.0 - g)));
            const C Dj = (b - p.rho * p.xi * iu - d) / xi2 * (1.0 - e) / (1.0 - g * e);
            const C phi = std::exp(Cj + Dj * p.v0 + iu * std::log(p.S0));
            return std::real(std::exp(-iu * std::log(K)) * phi / iu);
        };
        const double du = 0.01, uMax = 200.0;            // trapezoid rule
        double sum = 0.5 * integrand(1e-8);
        for (double u = du; u < uMax; u += du) sum += integrand(u);
        return 0.5 + sum * du / M_PI;
    };
    return p.S0 * std::exp(-p.q * T) * P(1) - K * std::exp(-p.r * T) * P(2);
}

// Margrabe: price of max(S1(T) - S2(T), 0)
double margrabe(const Market& a, const Market& b, double rho, double T)
{
    const double var = a.sigma * a.sigma + b.sigma * b.sigma - 2.0 * rho * a.sigma * b.sigma;
    const double sd  = std::sqrt(var * T);
    const double F1  = a.S0 * std::exp(-a.q * T), F2 = b.S0 * std::exp(-b.q * T);
    const double d1  = (std::log(F1 / F2) + 0.5 * sd * sd) / sd;
    return F1 * normCdf(d1) - F2 * normCdf(d1 - sd);
}

int main()
{
    const double T = 1.0, K = 100.0;
    const size_t batch = 4096;
    int failures = 0;

    // ---- 1. xi = 0, v0 = theta: constant variance, Heston = Black-Scholes. Cross-model correlation.
    {
        HestonParams a{100.0, 0.03, 0.01, 0.04, 2.0, 0.04, 0.0, -0.7};     // sigma = 20%
        HestonParams b{95.0, 0.03, 0.02, 0.09, 2.0, 0.09, 0.0, -0.6};      // sigma = 30%
        const double rhoS = 0.6;
        Models models;
        models["AAPL"] = std::make_unique<HestonEuler>(a, 1.0);          // exact in this case: 1 step
        models["MSFT"] = std::make_unique<HestonEuler>(b, 1.0);
        PayoffExchange payoff(T, "AAPL", "MSFT");
        // The user gives the correlations BETWEEN the models: 4 pairs (the 2 inside are the models')
        const std::vector<Correlation> cross{{"AAPL.0", "MSFT.0", rhoS}, {"AAPL.1", "MSFT.1", 0.5},
                                             {"AAPL.0", "MSFT.1", -0.35}, {"AAPL.1", "MSFT.0", -0.4}};
        const McResult res = mcPrice(a.r, models, cross, payoff, 1'000'000, batch, 42);
        const double exact = margrabe({a.S0, a.r, a.q, std::sqrt(a.v0)}, {b.S0, b.r, b.q, std::sqrt(b.v0)}, rhoS, T);
        const double z = (res.price - exact) / res.stdErr;
        std::cout << "xi = 0, exchange option: MC = " << res.price << " +/- " << res.stdErr
                  << "  Margrabe = " << exact << "  (" << z << " std errors)\n";
        if (std::abs(z) > 4.0) ++failures;
    }

    // ---- 2. European call on AAPL, full Heston: dynamics and correlation inside the model
    {
        const HestonParams a{100.0, 0.03, 0.01, 0.04, 2.0, 0.04, 0.4, -0.7};
        const double exact = hestonCall(a, K, T);
        std::cout << "\nEuropean call on AAPL (Heston), exact = " << exact << "\n";
        for (double maxStep : {0.1, 0.02, 0.005}) {
            Models models;
            models["AAPL"] = std::make_unique<HestonEuler>(a, maxStep);
            PayoffCall payoff(K, T, "AAPL");
            const McResult res = mcPrice(a.r, models, {}, payoff, 1'000'000, batch, 42);   // no user correlation: one model
            const double z = (res.price - exact) / res.stdErr;
            std::cout << "  maxStep = " << maxStep << ": MC = " << res.price << " +/- " << res.stdErr
                      << "  bias = " << z << " std errors\n";
            if (maxStep <= 0.005 && std::abs(z) > 4.0) ++failures;
        }
    }

    // ---- 3. exchange option, both Heston: the correlation between the assets is a price driver
    {
        const HestonParams a{100.0, 0.03, 0.01, 0.04, 2.0, 0.04, 0.4, -0.7};
        const HestonParams b{95.0, 0.03, 0.02, 0.09, 1.5, 0.09, 0.5, -0.6};
        std::cout << "\nExchange option AAPL / MSFT (Heston):\n";
        for (double rhoS : {0.6, 0.9}) {
            Models models;
            models["AAPL"] = std::make_unique<HestonEuler>(a, 0.02);
            models["MSFT"] = std::make_unique<HestonEuler>(b, 0.02);
            PayoffExchange payoff(T, "AAPL", "MSFT");
            // spot-spot rhoS, variance-variance 0.5, spot-variance across assets: rhoS times the own rho
            const std::vector<Correlation> cross{{"AAPL.0", "MSFT.0", rhoS}, {"AAPL.1", "MSFT.1", 0.5},
                                                 {"AAPL.0", "MSFT.1", rhoS * b.rho}, {"AAPL.1", "MSFT.0", rhoS * a.rho}};
            const McResult res = mcPrice(a.r, models, cross, payoff, 500'000, batch, 42);
            std::cout << "  rho(S1,S2) = " << rhoS << ": MC = " << res.price << " +/- " << res.stdErr << "\n";
        }
    }
    return failures;
}
