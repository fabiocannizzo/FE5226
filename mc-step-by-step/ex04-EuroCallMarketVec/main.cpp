// ex04 - European options, vectorized: the MARKET and the MODEL become objects.
//
// Same vectorized engine and payoff objects as ex03. Two changes:
//   1. The market data (spot, rate, yield, volatility) are bundled in a Market struct.
//   2. The dynamics of the underlying move out of the pricer into a model object, behind the
//      interface IModelMC. A model says
//        - how many Brownian motions it needs                  numberOfBrownians()
//        - how many state variables it carries                 numberOfStateVariables()
//        - the states at the initial date                      initializeStates()
//        - how to advance the states by ONE time step          generateStates()
//        - how to rebuild the spot from the states             spot()
//      The concrete model here is BSExact (Black-Scholes, exact step) whose state variable is the
//      Brownian motion W, not the spot. The payoff asks the model for the spot (spot()); the
//      pricer only talks to IModelMC: another model (Euler scheme, local volatility, ...)
//      needs no change in it.
//
// Data layout: states and random numbers are matrices with one ROW per state variable (or per
// Brownian) and one COLUMN per path of the batch, so that each row is contiguous in memory.
// The payoff reads the SPOT, which the model reconstructs from its states.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

struct Market {
    double S0;     // spot today
    double r;      // continuously compounded risk free rate (flat)
    double q;      // continuous dividend yield / foreign rate (flat); the drift is r - q
    double sigma;  // Black-Scholes volatility (flat)
};

// Row-major matrix: element (i,j) is data[i*cols + j]; row i is contiguous
struct Matrix {
    Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), data(rows * cols) {}
    double&                 operator()(size_t i, size_t j)       { return data[i * cols + j]; }
    double                  operator()(size_t i, size_t j) const { return data[i * cols + j]; }
    std::span<double>       row(size_t i)       { return {data.data() + i * cols, cols}; }
    std::span<const double> row(size_t i) const { return {data.data() + i * cols, cols}; }
    size_t rows, cols;
    std::vector<double> data;
};

struct McResult {
    double price;
    double stdErr;
};

// ---------------------------------------------------------------------------------------------
// Model interface and the concrete Black-Scholes model
// ---------------------------------------------------------------------------------------------
struct IModelMC {
    virtual ~IModelMC() = default;

    virtual size_t numberOfBrownians() const = 0;       // rows of the random number matrix
    virtual size_t numberOfStateVariables() const = 0;  // rows of the state matrices

    // states at the initial date, for the first n paths of the batch
    virtual void initializeStates(Matrix& states, size_t n) const = 0;

    // One time step from time t0 to time t1, for the first n paths:
    //   rnd  : numberOfBrownians() x paths      independent N(0,1) draws of the step
    //   prev : numberOfStateVariables() x paths states at t0
    //   next : numberOfStateVariables() x paths states at t1  (output)
    virtual void generateStates(double t0, double t1, const Matrix& rnd,
                                const Matrix& prev, Matrix& next, size_t n) const = 0;

    // Spot at time t for the first n paths, rebuilt from the states at t: out[i] = S(t) on path i
    virtual void spot(double t, const Matrix& states, size_t n, std::span<double> out) const = 0;
};

// Geometric Brownian motion, simulated EXACTLY on any step. The state variable is the Brownian
// motion W (not the spot): W has independent normal increments, so
//     W(t1) = W(t0) + sqrt(t1 - t0) Z
// is exact on any grid. The spot is rebuilt on demand:
//     S(t)  = S0 exp( (r - q - sigma^2/2) t + sigma W(t) )
// One Brownian, one state variable.
struct BSExact : IModelMC {
    explicit BSExact(const Market& m) : m_spot0(m.S0), m_mu(m.r - m.q), m_sigma(m.sigma) {}

    size_t numberOfBrownians() const override      { return 1; }
    size_t numberOfStateVariables() const override { return 1; }

    void initializeStates(Matrix& states, size_t n) const override {
        for (size_t i = 0; i < n; ++i)
            states(0, i) = 0.0;                          // W(0) = 0
    }

    void generateStates(double t0, double t1, const Matrix& rnd,
                        const Matrix& prev, Matrix& next, size_t n) const override {
        const double sqrtDt = std::sqrt(t1 - t0);
        for (size_t i = 0; i < n; ++i)
            next(0, i) = prev(0, i) + sqrtDt * rnd(0, i);
    }

    void spot(double t, const Matrix& states, size_t n, std::span<double> out) const override {
        const double drift = (m_mu - 0.5 * m_sigma * m_sigma) * t;
        for (size_t i = 0; i < n; ++i)
            out[i] = m_spot0 * std::exp(drift + m_sigma * states(0, i));
    }

    double m_spot0, m_mu, m_sigma;
};

// ---------------------------------------------------------------------------------------------
// Payoff interface and concrete payoffs
// ---------------------------------------------------------------------------------------------
enum class CallPut { Call, Put };

struct IPayoff {
    virtual ~IPayoff() = default;
    // Maturity: when the payoff is paid, in years from today
    virtual double maturity() const = 0;
    // Payoffs paid at maturity, for a batch of n paths: out[i] = payoff on path i.
    // The payoff does not see the states: it asks the model for the spot, through spot().
    virtual void payoff(const IModelMC& model, const Matrix& states, size_t n,
                        std::span<double> out) const = 0;
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

    void payoff(const IModelMC& model, const Matrix& states, size_t n,
                std::span<double> out) const override {
        model.spot(m_T, states, n, out);                 // out[i] = S(T) on path i
        for (size_t i = 0; i < n; ++i)
            out[i] = std::max(m_sign * (out[i] - m_strike), 0.0);
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

    void payoff(const IModelMC& model, const Matrix& states, size_t n,
                std::span<double> out) const override {
        model.spot(m_T, states, n, out);                 // out[i] = S(T) on path i
        for (size_t i = 0; i < n; ++i)
            out[i] = m_sign * (out[i] - m_strike) > 0.0 ? 1.0 : 0.0;
    }

    double m_strike;
    double m_T;
    double m_sign;
};

// ---------------------------------------------------------------------------------------------
// Pricer: knows neither the product nor the dynamics, only the two interfaces
// ---------------------------------------------------------------------------------------------
McResult mcPrice(const Market& m, const IModelMC& model, const IPayoff& payoff,
                 size_t nPaths, size_t batchSize, unsigned seed)
{
    const double T = payoff.maturity();
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Work space allocated once, re-used by every batch
    Matrix rnd(model.numberOfBrownians(), batchSize);
    Matrix prev(model.numberOfStateVariables(), batchSize);
    Matrix next(model.numberOfStateVariables(), batchSize);
    std::vector<double> pay(batchSize);

    double sum = 0.0, sumSq = 0.0;
    for (size_t done = 0; done < nPaths; done += batchSize) {
        const size_t n = std::min(batchSize, nPaths - done);  // the last batch may be shorter

        // 1. random numbers: one row per Brownian
        for (size_t b = 0; b < rnd.rows; ++b)
            for (size_t i = 0; i < n; ++i)
                rnd(b, i) = normal(rng);

        // 2. model: initial states, then ONE step from today (0) to the maturity T
        model.initializeStates(prev, n);
        model.generateStates(0.0, T, rnd, prev, next, n);

        // 3. payoff: one call for the whole batch. It rebuilds the spot from the states
        //    through the model.
        payoff.payoff(model, next, n, std::span<double>(pay).first(n));

        // 4. statistics
        for (size_t i = 0; i < n; ++i) {
            sum   += pay[i];
            sumSq += pay[i] * pay[i];
        }
    }

    const double N    = static_cast<double>(nPaths);
    const double mean = sum / N;
    const double var  = (sumSq / N - mean * mean) * N / (N - 1.0);
    const double df   = std::exp(-m.r * T);     // the market provides the discounting
    return {df * mean, df * std::sqrt(var / N)};
}

// ---- closed forms, only used to validate the Monte Carlo ----
double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

double bsVanilla(const Market& m, double K, double T, CallPut cp)
{
    const double S0 = m.S0, r = m.r, q = m.q, sigma = m.sigma;
    const double F  = S0 * std::exp((r - q) * T);
    const double sd = sigma * std::sqrt(T);
    const double d1 = (std::log(F / K) + 0.5 * sd * sd) / sd;
    const double d2 = d1 - sd;
    const double df = std::exp(-r * T);
    return cp == CallPut::Call ? df * (F * normCdf(d1) - K * normCdf(d2))
                               : df * (K * normCdf(-d2) - F * normCdf(-d1));
}

double bsDigital(const Market& m, double K, double T, CallPut cp)
{
    const double S0 = m.S0, r = m.r, q = m.q, sigma = m.sigma;
    const double F  = S0 * std::exp((r - q) * T);
    const double sd = sigma * std::sqrt(T);
    const double d2 = (std::log(F / K) - 0.5 * sd * sd) / sd;
    const double df = std::exp(-r * T);
    return cp == CallPut::Call ? df * normCdf(d2) : df * normCdf(-d2);
}

int main()
{
    const Market m{100.0, 0.03, 0.01, 0.20};
    const double K = 105.0, T = 1.0;

    const BSExact model(m);

    struct Trade {
        std::string                name;
        std::unique_ptr<IPayoff>   payoff;
        double                     exact;
    };
    std::vector<Trade> trades;
    trades.push_back({"vanilla call", std::make_unique<PayoffVanilla>(K, T, CallPut::Call),
                      bsVanilla(m, K, T, CallPut::Call)});
    trades.push_back({"vanilla put ", std::make_unique<PayoffVanilla>(K, T, CallPut::Put),
                      bsVanilla(m, K, T, CallPut::Put)});
    trades.push_back({"digital call", std::make_unique<PayoffDigital>(K, T, CallPut::Call),
                      bsDigital(m, K, T, CallPut::Call)});
    trades.push_back({"digital put ", std::make_unique<PayoffDigital>(K, T, CallPut::Put),
                      bsDigital(m, K, T, CallPut::Put)});

    int failures = 0;
    for (const Trade& t : trades) {
        const McResult mc = mcPrice(m, model, *t.payoff, 1'000'000, 4096, 42);
        const double z = (mc.price - t.exact) / mc.stdErr;
        std::cout << t.name << ": MC = " << mc.price << " +/- " << mc.stdErr
                  << "  exact = " << t.exact << "  (" << z << " std errors)\n";
        if (std::abs(z) > 4.0) ++failures;
    }
    return failures;
}
