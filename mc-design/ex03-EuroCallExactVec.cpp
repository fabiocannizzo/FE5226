// ex03 - European options, exact simulation, VECTORIZED in batches.
//
// Same inputs, model and payoff objects as ex02; the only change is the organisation of the work.
// Instead of one path at a time we process a BATCH of paths at a time, in contiguous arrays,
// one stage after the other:
//   - random numbers : one loop fills a whole array
//   - model          : a simple loop over an array (the compiler can use SIMD instructions)
//   - payoff         : ONE call per batch, with arrays in and out (the call overhead is shared
//                      by all the paths of the batch; the payoff body is a simple loop)
//   - memory         : bounded by the batch size, whatever the total number of paths
// The IPayoff interface therefore changes: payoff() takes and returns a batch.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

struct McResult {
    double price;
    double stdErr;
};

// ---------------------------------------------------------------------------------------------
// Payoff interface and concrete payoffs
// ---------------------------------------------------------------------------------------------
enum class CallPut { Call, Put };

struct IPayoff {
    virtual ~IPayoff() = default;
    // Maturity: when the payoff is paid, in years from today
    virtual double maturity() const = 0;
    // Payoffs paid at maturity, for a batch of spots at maturity: out[i] = payoff(ST[i])
    virtual void payoff(std::span<const double> ST, std::span<double> out) const = 0;
};

// Vanilla option: max(S - K, 0) for a call, max(K - S, 0) for a put
struct PayoffVanilla : IPayoff {
    PayoffVanilla(double strike, double T, CallPut callPut)
        : m_strike(strike), m_T(T), m_sign(callPut == CallPut::Call ? 1.0 : -1.0)
    {
        if (strike <= 0.0) throw std::invalid_argument("strike must be positive");
        if (T <= 0.0)      throw std::invalid_argument("maturity must be positive");
    }

    double maturity() const override { return m_T; }

    void payoff(std::span<const double> ST, std::span<double> out) const override {
        for (size_t i = 0; i < ST.size(); ++i)
            out[i] = std::max(m_sign * (ST[i] - m_strike), 0.0);
    }

    double m_strike;
    double m_T;
    double m_sign;   // +1 call, -1 put
};

// Digital (cash-or-nothing) option: pays 1 if S > K (call) or S < K (put), else 0
struct PayoffDigital : IPayoff {
    PayoffDigital(double strike, double T, CallPut callPut)
        : m_strike(strike), m_T(T), m_sign(callPut == CallPut::Call ? 1.0 : -1.0)
    {
        if (strike <= 0.0) throw std::invalid_argument("strike must be positive");
        if (T <= 0.0)      throw std::invalid_argument("maturity must be positive");
    }

    double maturity() const override { return m_T; }

    void payoff(std::span<const double> ST, std::span<double> out) const override {
        for (size_t i = 0; i < ST.size(); ++i)
            out[i] = m_sign * (ST[i] - m_strike) > 0.0 ? 1.0 : 0.0;
    }

    double m_strike;
    double m_T;
    double m_sign;
};

// ---------------------------------------------------------------------------------------------
// Pricer: knows the market and the model, not the product
// ---------------------------------------------------------------------------------------------
McResult mcPrice(double S0, double r, double q, double sigma,
                 const IPayoff& payoff, size_t nPaths, size_t batchSize, unsigned seed)
{
    const double T = payoff.maturity();
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    const double drift = (r - q - 0.5 * sigma * sigma) * T;
    const double stdev = sigma * std::sqrt(T);

    // Work space allocated once, re-used by every batch
    std::vector<double> z(batchSize), ST(batchSize), pay(batchSize);

    double sum = 0.0, sumSq = 0.0;
    for (size_t done = 0; done < nPaths; done += batchSize) {
        const size_t n = std::min(batchSize, nPaths - done);  // the last batch may be shorter

        // 1. random numbers
        for (size_t i = 0; i < n; ++i)
            z[i] = normal(rng);

        // 2. model: terminal spot of every path (no dependency between iterations)
        for (size_t i = 0; i < n; ++i)
            ST[i] = S0 * std::exp(drift + stdev * z[i]);

        // 3. payoff: one call for the whole batch
        payoff.payoff(std::span<const double>(ST).first(n), std::span<double>(pay).first(n));

        // 4. statistics
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

// ---- closed forms, only used to validate the Monte Carlo ----
double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

double bsVanilla(double S0, double r, double q, double sigma, double K, double T, CallPut cp)
{
    const double F  = S0 * std::exp((r - q) * T);
    const double sd = sigma * std::sqrt(T);
    const double d1 = (std::log(F / K) + 0.5 * sd * sd) / sd;
    const double d2 = d1 - sd;
    const double df = std::exp(-r * T);
    return cp == CallPut::Call ? df * (F * normCdf(d1) - K * normCdf(d2))
                               : df * (K * normCdf(-d2) - F * normCdf(-d1));
}

double bsDigital(double S0, double r, double q, double sigma, double K, double T, CallPut cp)
{
    const double F  = S0 * std::exp((r - q) * T);
    const double sd = sigma * std::sqrt(T);
    const double d2 = (std::log(F / K) - 0.5 * sd * sd) / sd;
    const double df = std::exp(-r * T);
    return cp == CallPut::Call ? df * normCdf(d2) : df * normCdf(-d2);
}

int main()
{
    const double S0 = 100.0, r = 0.03, q = 0.01, sigma = 0.20;
    const double K = 105.0, T = 1.0;

    struct Trade {
        std::string                name;
        std::unique_ptr<IPayoff>   payoff;
        double                     exact;
    };
    std::vector<Trade> trades;
    trades.push_back({"vanilla call", std::make_unique<PayoffVanilla>(K, T, CallPut::Call),
                      bsVanilla(S0, r, q, sigma, K, T, CallPut::Call)});
    trades.push_back({"vanilla put ", std::make_unique<PayoffVanilla>(K, T, CallPut::Put),
                      bsVanilla(S0, r, q, sigma, K, T, CallPut::Put)});
    trades.push_back({"digital call", std::make_unique<PayoffDigital>(K, T, CallPut::Call),
                      bsDigital(S0, r, q, sigma, K, T, CallPut::Call)});
    trades.push_back({"digital put ", std::make_unique<PayoffDigital>(K, T, CallPut::Put),
                      bsDigital(S0, r, q, sigma, K, T, CallPut::Put)});

    int failures = 0;
    for (const Trade& t : trades) {
        const McResult mc = mcPrice(S0, r, q, sigma, *t.payoff, 1'000'000, 4096, 42);
        const double z = (mc.price - t.exact) / mc.stdErr;
        std::cout << t.name << ": MC = " << mc.price << " +/- " << mc.stdErr
                  << "  exact = " << t.exact << "  (" << z << " std errors)\n";
        if (std::abs(z) > 4.0) ++failures;
    }
    return failures;
}
