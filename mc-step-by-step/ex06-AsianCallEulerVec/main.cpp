// ex06 - Asian call, EULER discretization: the model is no longer exact.
//
// With a model that cannot be simulated exactly (local volatility, stochastic volatility, ...)
// the paths are built by small time steps, and the result carries a DISCRETIZATION BIAS that
// shrinks with the step size. Two new ideas, both visible in the code:
//   1. A Model object bundles what depends on the dynamics:
//        mergeDates : adds intermediate dates to the payoff dates (steps no longer than maxStep)
//        simulate   : builds the spots on that simulation grid
//   2. The payoff is untouched: it still sees only the spots at ITS dates. The engine simulates
//      on the fine grid and hands the payoff only the rows of the fixing dates.
// Scheme (Euler-Maruyama):  S_{k+1} = S_k (1 + (r-q) dt + sigma sqrt(dt) Z_k), truncated at 0.
//
// The bias is measured against the closed form of the geometric Asian call.
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <random>
#include <span>
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

using Payoff = std::function<void(const Matrix& S, std::span<double> out)>;

// The model: how to choose the simulation grid, and how to simulate on it.
// grid[k] are the simulation dates (> 0, increasing); S has one row per grid date, one column per path.
struct Model {
    std::function<std::vector<double>(const std::vector<double>& payoffDates)> mergeDates;
    std::function<void(const std::vector<double>& grid, size_t n, std::mt19937_64& rng, Matrix& S)> simulate;
};

// Exact model: any grid is fine, no date is added (as in ex05)
Model exactModel(const Market& m)
{
    return Model{
        [](const std::vector<double>& dates) { return dates; },
        [m](const std::vector<double>& grid, size_t n, std::mt19937_64& rng, Matrix& S) {
            std::normal_distribution<double> normal(0.0, 1.0);
            std::vector<double> W(n, 0.0);
            double tPrev = 0.0;
            for (size_t j = 0; j < grid.size(); ++j) {
                const double sq    = std::sqrt(grid[j] - tPrev);
                const double drift = (m.r - m.q - 0.5 * m.sigma * m.sigma) * grid[j];
                std::span<double> Sj = S.row(j);
                for (size_t i = 0; i < n; ++i) {
                    W[i]  += sq * normal(rng);
                    Sj[i]  = m.S0 * std::exp(drift + m.sigma * W[i]);
                }
                tPrev = grid[j];
            }
        }};
}

// Euler model: splits every interval of the payoff dates in equal steps no longer than maxStep
Model eulerModel(const Market& m, double maxStep)
{
    return Model{
        [maxStep](const std::vector<double>& dates) {
            std::vector<double> grid;
            double prev = 0.0;
            for (double t : dates) {
                const int nSub = static_cast<int>(std::ceil((t - prev) / maxStep - 1e-12));
                for (int k = 1; k < nSub; ++k)
                    grid.push_back(prev + (t - prev) * k / nSub);
                grid.push_back(t);                      // the payoff date itself, exactly
                prev = t;
            }
            return grid;
        },
        [m](const std::vector<double>& grid, size_t n, std::mt19937_64& rng, Matrix& S) {
            std::normal_distribution<double> normal(0.0, 1.0);
            double tPrev = 0.0;
            for (size_t j = 0; j < grid.size(); ++j) {
                const double dt = grid[j] - tPrev;
                const double a  = 1.0 + (m.r - m.q) * dt;
                const double b  = m.sigma * std::sqrt(dt);
                std::span<double> Sj = S.row(j);
                for (size_t i = 0; i < n; ++i) {
                    const double Sprev = (j == 0) ? m.S0 : S(j - 1, i);
                    Sj[i] = std::max(Sprev * (a + b * normal(rng)), 0.0);
                }
                tPrev = grid[j];
            }
        }};
}

McResult mcPrice(const Market& m, const Model& model, const std::vector<double>& payoffDates,
                 const Payoff& payoff, size_t nPaths, size_t batchSize, unsigned seed)
{
    const double T = payoffDates.back();
    std::mt19937_64 rng(seed);

    // The model decides the simulation grid; locate the payoff dates inside it
    const std::vector<double> grid = model.mergeDates(payoffDates);
    std::vector<size_t> obsRow;                         // obsRow[j] = row of grid for payoff date j
    for (double t : payoffDates)
        obsRow.push_back(static_cast<size_t>(std::lower_bound(grid.begin(), grid.end(), t - 1e-12) - grid.begin()));

    Matrix Sgrid(grid.size(), batchSize);               // spots on the simulation grid
    Matrix Sobs(payoffDates.size(), batchSize);         // spots on the payoff dates only
    std::vector<double> pay(batchSize);

    double sum = 0.0, sumSq = 0.0;
    for (size_t done = 0; done < nPaths; done += batchSize) {
        const size_t n = std::min(batchSize, nPaths - done);

        model.simulate(grid, n, rng, Sgrid);
        for (size_t j = 0; j < obsRow.size(); ++j)      // keep only what the payoff asks for
            std::copy_n(Sgrid.row(obsRow[j]).begin(), n, Sobs.row(j).begin());

        payoff(Sobs, std::span<double>(pay).first(n));
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
Payoff asianCall(double K)
{
    return [K](const Matrix& S, std::span<double> out) {
        const size_t n = out.size();
        std::vector<double> avg(n, 0.0);
        for (size_t j = 0; j < S.rows; ++j)             // sum over dates, inner loop over paths
            for (size_t i = 0; i < n; ++i)
                avg[i] += S(j, i);
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(avg[i] / static_cast<double>(S.rows) - K, 0.0);
    };
}

Payoff geometricAsianCall(double K)
{
    return [K](const Matrix& S, std::span<double> out) {
        const size_t n = out.size();
        std::vector<double> logSum(n, 0.0);
        for (size_t j = 0; j < S.rows; ++j)
            for (size_t i = 0; i < n; ++i)
                logSum[i] += std::log(S(j, i));
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(std::exp(logSum[i] / static_cast<double>(S.rows)) - K, 0.0);
    };
}

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
    const double exact = geometricAsianCallExact(m, quarterly, K);
    std::cout << "geometric Asian call, 4 fixings, exact = " << exact << "\n\n"
              << std::fixed << std::setprecision(4);

    int failures = 0;
    {
        const McResult r = mcPrice(m, exactModel(m), quarterly, geometricAsianCall(K), nPaths, batch, 42);
        const double z = (r.price - exact) / r.stdErr;
        std::cout << "exact model          : " << r.price << " +/- " << r.stdErr << "  bias = " << z << " std errors\n";
        if (std::abs(z) > 4.0) ++failures;
    }
    for (double maxStep : {0.25, 0.10, 0.05, 0.01}) {
        const McResult r = mcPrice(m, eulerModel(m, maxStep), quarterly, geometricAsianCall(K), nPaths, batch, 42);
        const double z = (r.price - exact) / r.stdErr;
        std::cout << "Euler, maxStep = " << std::setw(4) << maxStep << " : " << r.price << " +/- " << r.stdErr
                  << "  bias = " << z << " std errors\n";
        if (maxStep <= 0.01 && std::abs(z) > 4.0) ++failures;
    }

    const McResult asian = mcPrice(m, eulerModel(m, 0.01), quarterly, asianCall(K), nPaths, batch, 42);
    std::cout << "\narithmetic Asian call (Euler, maxStep = 0.01): " << asian.price << " +/- " << asian.stdErr << "\n";
    return failures;
}
