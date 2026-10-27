#include "MCEngine.h"
#include "Models.h"
#include "Payoffs.h"
#include "Pricers.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

using namespace mc;

namespace {

// Standard normal CDF via std::erfc
double normCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// Analytical Black-Scholes European Call / Put formula for benchmark validation
double blackScholesPrice(double S0,
                         double K,
                         double r_disc,
                         double mu,
                         double sigma,
                         double T,
                         int    omega,
                         double notional = 1.0) {
    double df     = std::exp(-r_disc * T);
    double F      = S0 * std::exp(mu * T);
    double vol_sq = sigma * std::sqrt(T);
    double d1     = (std::log(F / K) + 0.5 * sigma * sigma * T) / vol_sq;
    double d2     = d1 - vol_sq;
    return notional * df * (omega * F * normCDF(omega * d1) - omega * K * normCDF(omega * d2));
}

} // namespace

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "====================================================================\n";
    std::cout << "  Multi-Model MC Engine: Closed-Form vs. Euler Discretization\n";
    std::cout << "  Diffusion vs. Reconstruction & Model Factory (FE5226)\n";
    std::cout << "====================================================================\n\n";

    const Date   today   = 0;
    const size_t n_paths = 200000;

    // ------------------------------------------------------------------------
    // 1. Observed Historical Events (Fixing Events up to Today)
    // ------------------------------------------------------------------------
    Events events = {
        makeFixingEvent("GBP.USD", -90, 1.25), // Historical fixing at t = -90d
        makeFixingEvent("GBP.USD",   0, 1.28)  // Spot fixing at t = today (0d)
    };

    // ------------------------------------------------------------------------
    // 2. Payoffs (Pure Contract Specification & Pricer Factory)
    // ------------------------------------------------------------------------
    // Trade 0: 6M Vanilla Call on EUR.USD
    auto vanilla_call = std::make_shared<PayoffVanilla>(
        "EUR.USD", /*notional*/ 1000000.0, /*strike*/ 1.10, /*Call*/ +1,
        /*expiry*/ 180, /*settle*/ 180, "USD");

    // Trade 1: Seasoned Asian Call on GBP.USD
    // Schedule: -90d (past), 0d (past), 90d (future), 180d (future)
    std::vector<Date> asian_sched = { -90, 0, 90, 180 };
    auto asian_call = std::make_shared<PayoffAsian>(
        "GBP.USD", /*notional*/ 1000000.0, /*strike*/ 1.26, /*Call*/ +1,
        asian_sched, /*settle*/ 180, "USD");

    // Trade 2: Worst-Of Call on EUR.USD and GBP.USD
    auto worst_of_call = std::make_shared<PayoffWorstOfCall>(
        "EUR.USD", "GBP.USD", /*spot0*/ 1.10, 1.28,
        /*notional*/ 1000000.0, /*strike*/ 1.00,
        /*expiry*/ 180, /*settle*/ 180, "USD");

    // Portfolio of payoffs passed to the pricing function
    std::vector<ppayoff_t> portfolio = {
        vanilla_call,
        asian_call,
        worst_of_call
    };

    // Inspect IPayoff contract query methods
    std::cout << "[IPayoff Contract Query Methods with Events]\n";
    std::cout << "  Asian Call observationDates() : { ";
    auto obs = asian_call->observationDates(events);
    for (size_t i = 0; i < obs.size(); ++i) {
        std::cout << obs[i] << (i + 1 < obs.size() ? ", " : " }\n");
    }
    std::cout << "  Asian Call maturityDate()     : " << asian_call->maturityDate(events) << "\n";
    std::cout << "  Asian Call settlementDates()  : { " << asian_call->settlementDates(events)[0] << " }\n";
    std::cout << "  Asian Call cashflow ccy       : { " << asian_call->getCashFlowCurrencies(events)[0] << " }\n\n";

    // ------------------------------------------------------------------------
    // 3. Configurations
    // ------------------------------------------------------------------------
    PricerConfig config_closed{ "BS_CLOSED" };
    PricerConfig config_euler{ "BS_EULER" };

    // ------------------------------------------------------------------------
    // 4. Prepopulated Market Data Server
    // ------------------------------------------------------------------------
    const double r_usd     = 0.04;
    const double mu_eurusd = 0.015;
    const double vol_eur   = 0.12;
    const double mu_gbpusd = 0.010;
    const double vol_gbp   = 0.14;

    FlatDiscountCurve usd_curve("IR.DISCOUNT.USD", today, r_usd);

    MarketDataServer mds;
    mds.spots["EUR.USD"]  = 1.10;
    mds.spots["GBP.USD"]  = 1.28;
    mds.drifts["EUR.USD"] = mu_eurusd;
    mds.drifts["GBP.USD"] = mu_gbpusd;
    mds.vols["EUR.USD"]   = vol_eur;
    mds.vols["GBP.USD"]   = vol_gbp;

    std::string corr_closed_key = corrKey("BS_CLOSED.EUR.USD.0", "BS_CLOSED.GBP.USD.0");
    std::string corr_euler_key  = corrKey("BS_EULER.EUR.USD.0", "BS_EULER.GBP.USD.0");
    mds.correlations[corr_closed_key] = 0.50;
    mds.correlations[corr_euler_key]  = 0.50;
    mds.discount_curves["USD"]        = &usd_curve;

    // ------------------------------------------------------------------------
    // 5. Valuation: Pricing function taking a vector of payoffs directly!
    // ------------------------------------------------------------------------
    std::cout << "[Pricing Function taking a Vector of Payoffs]\n";
    auto pvs_closed = price(portfolio, events, config_closed, mds, today, n_paths, /*seed*/ 42, /*antithetic*/ true);
    auto pvs_euler  = price(portfolio, events, config_euler,  mds, today, n_paths, /*seed*/ 42, /*antithetic*/ true);

    double bs_exact = blackScholesPrice(
        1.10, 1.10, r_usd, mu_eurusd, vol_eur, 180.0 / 365.0, +1, 1000000.0);

    std::cout << "[Pricing Results (N = " << n_paths << " paths)]\n";
    std::cout << "  Trade 0: EUR Vanilla Call (BS Closed) : " << std::setw(12) << pvs_closed[0]
              << " USD  (Analytical BS: " << bs_exact << ", diff: " << std::abs(pvs_closed[0] - bs_exact) << ")\n";
    std::cout << "  Trade 0: EUR Vanilla Call (BS Euler)  : " << std::setw(12) << pvs_euler[0]
              << " USD  (Euler vs Closed diff: " << std::abs(pvs_euler[0] - pvs_closed[0]) << ")\n";
    std::cout << "  Trade 1: GBP Asian Call   (BS Closed) : " << std::setw(12) << pvs_closed[1] << " USD\n";
    std::cout << "  Trade 1: GBP Asian Call   (BS Euler)  : " << std::setw(12) << pvs_euler[1]
              << " USD  (Euler vs Closed diff: " << std::abs(pvs_euler[1] - pvs_closed[1]) << ")\n";
    std::cout << "  Trade 2: EUR/GBP Worst-Of (BS Closed) : " << std::setw(12) << pvs_closed[2] << " USD\n";
    std::cout << "  Trade 2: EUR/GBP Worst-Of (BS Euler)  : " << std::setw(12) << pvs_euler[2]
              << " USD  (Euler vs Closed diff: " << std::abs(pvs_euler[2] - pvs_closed[2]) << ")\n\n";

    // ------------------------------------------------------------------------
    // 6. Common Random Numbers (CRN) Finite-Difference Delta on EUR.USD
    // ------------------------------------------------------------------------
    const double S0   = 1.10;
    const double bump = 1e-4;

    mds.spots["EUR.USD"] = S0 + bump;
    auto pv_up           = price(portfolio, events, config_closed, mds, today, n_paths, 42, true);
    mds.spots["EUR.USD"] = S0 - bump;
    auto pv_dn           = price(portfolio, events, config_closed, mds, today, n_paths, 42, true);
    mds.spots["EUR.USD"] = S0;

    double mc_delta = (pv_up[0] - pv_dn[0]) / (2.0 * bump);

    double T      = 180.0 / 365.0;
    double d1     = (std::log(1.0) + (mu_eurusd + 0.5 * vol_eur * vol_eur) * T) / (vol_eur * std::sqrt(T));
    double bs_del = 1000000.0 * std::exp((mu_eurusd - r_usd) * T) * normCDF(d1);

    std::cout << "[Common Random Numbers (CRN) Finite-Difference Delta (dS = " << bump << ")]\n";
    std::cout << "  EUR.USD Vanilla Call MC Delta : " << std::setw(12) << mc_delta
              << "  (Analytical BS Delta: " << bs_del << ")\n";

    return 0;
}
