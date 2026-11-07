// ex01 - European call, Monte Carlo, exact simulation, one path at a time.
//
// The simplest possible pricer, everything hardcoded and passed as plain numbers:
//   inputs  : S0, r, q, sigma, K, T
//   model   : S(T) = S0 exp((r - q - sigma^2/2) T + sigma sqrt(T) Z),  Z ~ N(0,1)   (exact)
//   payoff  : max(S(T) - K, 0), written inside the path loop
//   price   : discounted average of the payoffs, with its standard error
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

struct McResult {
    double price;
    double stdErr;
};

McResult mcEuroCall(double S0, double r, double q, double sigma, double K, double T,
                    size_t nPaths, unsigned seed)
{
    std::mt19937_64 rng(seed); // seed the generator for repeatability
    std::normal_distribution<double> normal(0.0, 1.0);

    const double drift = (r - q - 0.5 * sigma * sigma) * T;  // computed once for efficiency
    const double stdev = sigma * std::sqrt(T);               // computed once for efficiency

    double sum = 0.0, sumSq = 0.0;
    for (size_t i = 0; i < nPaths; ++i) {
        const double z  = normal(rng);                              // random number N(0,1)
        const double ST = S0 * std::exp(drift + stdev * z);         // price simulated from the model
        const double p  = std::max(ST - K, 0.0);                    // call payoff
        sum += p;  sumSq += p * p;                                  // accumulate statistics
    }

    const double n    = static_cast<double>(nPaths);
    const double mean = sum / n;
    const double var  = (sumSq / n - mean * mean) * n / (n - 1.0);   // sample variance
    const double df   = std::exp(-r * T);
    return {df * mean, df * std::sqrt(var / n)};
}

// ---- closed form, only used to validate the Monte Carlo ----
double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

double bsCall(double S0, double r, double q, double sigma, double K, double T)
{
    const double F  = S0 * std::exp((r - q) * T);
    const double sd = sigma * std::sqrt(T);
    const double d1 = (std::log(F / K) + 0.5 * sd * sd) / sd;
    return std::exp(-r * T) * (F * normCdf(d1) - K * normCdf(d1 - sd));
}

int main()
{
    const double S0 = 100.0, r = 0.03, q = 0.01, sigma = 0.20;
    const double K = 105.0, T = 1.0;

    const McResult mc    = mcEuroCall(S0, r, q, sigma, K, T, 1'000'000, 42);
    const double   exact = bsCall(S0, r, q, sigma, K, T);

    std::cout << "MC    = " << mc.price << "  +/- " << mc.stdErr << "\n"
              << "exact = " << exact << "  (" << (mc.price - exact) / mc.stdErr << " std errors)\n";
    return std::abs(mc.price - exact) < 4.0 * mc.stdErr ? 0 : 1;
}
