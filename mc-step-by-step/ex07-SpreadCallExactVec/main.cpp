// ex07 - Spread call on TWO correlated assets, exact simulation, vectorized.
//
//     payoff = max( S1(T) - S2(T) - K , 0 )
//
// New ingredient: several risk factors. Each asset has its own Brownian motion; the Brownian
// motions are correlated with a matrix R. Correlated normals are obtained from independent ones:
//     R = L L^T  (Cholesky),       Zc = L Z,      Z independent N(0,1)
// A path is now a LIST of matrices, one per asset: S[a](j,i) = spot of asset a, date j, path i.
// The payoff takes the list; it can combine assets at any date (here only the last one matters).
//
// Validation: with K = 0 the spread call is an exchange option, with Margrabe's closed form.
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <span>
#include <vector>

struct Asset { double S0, q, sigma; };                  // one underlying; the rate r is common
struct Market {
    std::vector<Asset>               assets;
    std::vector<std::vector<double>> corr;              // correlation matrix of the Brownians
    double                           r;
};
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

// S[a] is the matrix (dates x paths) of asset a
using Payoff = std::function<void(const std::vector<Matrix>& S, std::span<double> out)>;

// Lower triangular L with L L^T = R (R symmetric positive definite)
std::vector<std::vector<double>> cholesky(const std::vector<std::vector<double>>& R)
{
    const size_t n = R.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j <= i; ++j) {
            double x = R[i][j];
            for (size_t k = 0; k < j; ++k) x -= L[i][k] * L[j][k];
            L[i][j] = (i == j) ? std::sqrt(x) : x / L[j][j];
        }
    return L;
}

McResult mcPrice(const Market& m, const std::vector<double>& dates, const Payoff& payoff,
                 size_t nPaths, size_t batchSize, unsigned seed)
{
    const size_t nA = m.assets.size(), nDates = dates.size();
    const double T  = dates.back();
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    const auto L = cholesky(m.corr);

    std::vector<Matrix> S(nA, Matrix(nDates, batchSize));
    std::vector<std::vector<double>> W(nA, std::vector<double>(batchSize));   // Brownian of each asset
    std::vector<std::vector<double>> Z(nA, std::vector<double>(batchSize));   // independent normals
    std::vector<double> pay(batchSize);

    double sum = 0.0, sumSq = 0.0;
    for (size_t done = 0; done < nPaths; done += batchSize) {
        const size_t n = std::min(batchSize, nPaths - done);
        for (auto& w : W) std::fill(w.begin(), w.end(), 0.0);

        double tPrev = 0.0;
        for (size_t j = 0; j < nDates; ++j) {
            const double sq = std::sqrt(dates[j] - tPrev);
            for (size_t a = 0; a < nA; ++a)             // independent normals
                for (size_t i = 0; i < n; ++i) Z[a][i] = normal(rng);

            for (size_t a = 0; a < nA; ++a) {           // correlate, advance W, rebuild the spot
                const Asset& as = m.assets[a];
                const double drift = (m.r - as.q - 0.5 * as.sigma * as.sigma) * dates[j];
                std::span<double> Sj = S[a].row(j);
                for (size_t i = 0; i < n; ++i) {
                    double zc = 0.0;
                    for (size_t b = 0; b <= a; ++b) zc += L[a][b] * Z[b][i];   // L is triangular
                    W[a][i] += sq * zc;
                    Sj[i]    = as.S0 * std::exp(drift + as.sigma * W[a][i]);
                }
            }
            tPrev = dates[j];
        }

        payoff(S, std::span<double>(pay).first(n));
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

Payoff spreadCall(double K)
{
    return [K](const std::vector<Matrix>& S, std::span<double> out) {
        const size_t last = S[0].rows - 1;              // the last date
        for (size_t i = 0; i < out.size(); ++i)
            out[i] = std::max(S[0](last, i) - S[1](last, i) - K, 0.0);
    };
}

double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Margrabe: price of max(S1(T) - S2(T), 0)
double margrabe(const Market& m, double T)
{
    const Asset &a = m.assets[0], &b = m.assets[1];
    const double rho = m.corr[0][1];
    const double var = a.sigma * a.sigma + b.sigma * b.sigma - 2.0 * rho * a.sigma * b.sigma;
    const double sd  = std::sqrt(var * T);
    const double F1  = a.S0 * std::exp(-a.q * T), F2 = b.S0 * std::exp(-b.q * T);
    const double d1  = (std::log(F1 / F2) + 0.5 * sd * sd) / sd;
    return F1 * normCdf(d1) - F2 * normCdf(d1 - sd);
}

int main()
{
    const Market m{{{100.0, 0.01, 0.20}, {95.0, 0.02, 0.30}}, {{1.0, 0.6}, {0.6, 1.0}}, 0.03};
    const double T = 1.0;
    const std::vector<double> dates{T};
    const size_t nPaths = 1'000'000, batch = 4096;

    const McResult ex = mcPrice(m, dates, spreadCall(0.0), nPaths, batch, 42);
    const double exact = margrabe(m, T);
    const double z = (ex.price - exact) / ex.stdErr;
    std::cout << "exchange option: MC = " << ex.price << " +/- " << ex.stdErr
              << "  Margrabe = " << exact << "  (" << z << " std errors)\n";

    const McResult sp = mcPrice(m, dates, spreadCall(5.0), nPaths, batch, 42);
    std::cout << "spread call K=5: MC = " << sp.price << " +/- " << sp.stdErr << "\n";

    // The correlation matters: a higher correlation lowers the volatility of the spread
    Market m2 = m;
    m2.corr = {{1.0, 0.9}, {0.9, 1.0}};
    const McResult sp2 = mcPrice(m2, dates, spreadCall(5.0), nPaths, batch, 42);
    std::cout << "spread call K=5, rho = 0.9: MC = " << sp2.price << " +/- " << sp2.stdErr << "\n";

    return std::abs(z) < 4.0 ? 0 : 1;
}
