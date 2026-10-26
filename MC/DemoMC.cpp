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
    std::vector<Events> events = {
        Events("GBP.USD", -90, 1.25), // Historical fixing at t = -90d
        Events("GBP.USD",   0, 1.28)  // Spot fixing at t = today (0d)
    };

    // ------------------------------------------------------------------------
    // 2. Payoffs (Pure Contract Specification & Pricer Factory)
    //    Notice: IPayoff has no fixingDates() and no eval()!
    // ------------------------------------------------------------------------
    // Trade 0: 6M Vanilla Call on EUR.USD under BS Closed-Form (Brownian state variable)
    auto vanilla_call_closed = std::make_shared<PayoffVanilla>(
        "EUR.USD", /*notional*/ 1000000.0, /*strike*/ 1.10, /*Call*/ +1,
        /*expiry*/ 180, /*settle*/ 180, "USD", "BS_CLOSED");

    // Trade 1: 6M Vanilla Call on EUR.USD under BS Euler (Spot state variable + monthly step injection)
    auto vanilla_call_euler = std::make_shared<PayoffVanilla>(
        "EUR.USD", /*notional*/ 1000000.0, /*strike*/ 1.10, /*Call*/ +1,
        /*expiry*/ 180, /*settle*/ 180, "USD", "BS_EULER");

    // Trade 2: Seasoned Asian Call on GBP.USD under BS Closed-Form
    // Schedule: -90d (past), 0d (past), 90d (future), 180d (future)
    std::vector<Date> asian_sched = { -90, 0, 90, 180 };
    auto asian_call_closed = std::make_shared<PayoffAsian>(
        "GBP.USD", /*notional*/ 1000000.0, /*strike*/ 1.26, /*Call*/ +1,
        asian_sched, /*settle*/ 180, "USD", "BS_CLOSED");

    // Trade 3: Seasoned Asian Call on GBP.USD under BS Euler
    auto asian_call_euler = std::make_shared<PayoffAsian>(
        "GBP.USD", /*notional*/ 1000000.0, /*strike*/ 1.26, /*Call*/ +1,
        asian_sched, /*settle*/ 180, "USD", "BS_EULER");

    // Trade 4: Worst-Of Call on EUR.USD and GBP.USD
    auto worst_of_call = std::make_shared<PayoffWorstOfCall>(
        "EUR.USD", "GBP.USD", /*spot0*/ 1.10, 1.28,
        /*notional*/ 1000000.0, /*strike*/ 1.00,
        /*expiry*/ 180, /*settle*/ 180, "USD", "BS_CLOSED", "BS_CLOSED");

    // ------------------------------------------------------------------------
    // 3. Factory Generation: payoff->getPricer(events)
    // ------------------------------------------------------------------------
    std::vector<pmcpricer_t> pricers = {
        vanilla_call_closed->getPricer(events), // Pricer 0: Vanilla Call (BS Closed)
        vanilla_call_euler->getPricer(events),  // Pricer 1: Vanilla Call (BS Euler)
        asian_call_closed->getPricer(events),   // Pricer 2: Asian Call (BS Closed)
        asian_call_euler->getPricer(events),    // Pricer 3: Asian Call (BS Euler)
        worst_of_call->getPricer(events)        // Pricer 4: Worst-Of Call (BS Closed)
    };

    std::cout << "[Pricer Lifecycling via getPricer(events)]\n";
    std::cout << "  Asian Call contract schedule      : { -90, 0, 90, 180 }\n";
    std::cout << "  Asian Call requiredDates() (future): { ";
    for (size_t i = 0; i < pricers[2]->requiredDates().size(); ++i) {
        std::cout << pricers[2]->requiredDates()[i]
                  << (i + 1 < pricers[2]->requiredDates().size() ? ", " : " }\n\n");
    }

    // ------------------------------------------------------------------------
    // 4. MCEngine Master Setup
    // ------------------------------------------------------------------------
    MCEngine engine(pricers, today, n_paths, /*seed*/ 42, /*antithetic*/ true);

    std::cout << "[Engine Initial Setup (Before Model Date Injection)]\n";
    std::cout << "  Pricer Required Dates: { ";
    for (size_t i = 0; i < engine.masterDates().size(); ++i) {
        std::cout << engine.masterDates()[i] << (i + 1 < engine.masterDates().size() ? ", " : " }\n");
    }
    std::cout << "  Models Needed        : { ";
    for (size_t i = 0; i < engine.masterModels().size(); ++i) {
        std::cout << engine.masterModels()[i] << (i + 1 < engine.masterModels().size() ? ", " : " }\n\n");
    }

    // ------------------------------------------------------------------------
    // 5. Prepopulated Market Data Server
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
    // 6. Valuation: Engine Queries Model Factory & Injects Dates
    // ------------------------------------------------------------------------
    std::vector<double> pvs = engine.price(mds);

    std::cout << "[Timeline After Model Date Injection (Euler Injected <= 30d Steps)]\n";
    std::cout << "  Merged Timeline (" << engine.masterDates().size() << " dates): { ";
    for (size_t i = 0; i < engine.masterDates().size(); ++i) {
        std::cout << engine.masterDates()[i] << (i + 1 < engine.masterDates().size() ? ", " : " }\n\n");
    }

    double bs_exact = blackScholesPrice(
        1.10, 1.10, r_usd, mu_eurusd, vol_eur, 180.0 / 365.0, +1, 1000000.0);

    std::cout << "[Pricing Results (N = " << n_paths << " paths)]\n";
    std::cout << "  Pricer 0: EUR Vanilla Call (BS Closed, W state) : " << std::setw(12) << pvs[0]
              << " USD  (Analytical BS: " << bs_exact << ", diff: " << std::abs(pvs[0] - bs_exact) << ")\n";
    std::cout << "  Pricer 1: EUR Vanilla Call (BS Euler,  S state) : " << std::setw(12) << pvs[1]
              << " USD  (Euler vs Closed diff: " << std::abs(pvs[1] - pvs[0]) << ")\n";
    std::cout << "  Pricer 2: GBP Asian Call   (BS Closed, W state) : " << std::setw(12) << pvs[2] << " USD\n";
    std::cout << "  Pricer 3: GBP Asian Call   (BS Euler,  S state) : " << std::setw(12) << pvs[3]
              << " USD  (Euler vs Closed diff: " << std::abs(pvs[3] - pvs[2]) << ")\n";
    std::cout << "  Pricer 4: EUR/GBP Worst-Of Call (BS Closed)     : " << std::setw(12) << pvs[4] << " USD\n\n";

    // ------------------------------------------------------------------------
    // 7. Common Random Numbers (CRN) Finite-Difference Delta on EUR.USD
    // ------------------------------------------------------------------------
    const double S0   = 1.10;
    const double bump = 1e-4;

    mds.spots["EUR.USD"] = S0 + bump;
    auto pv_up           = engine.price(mds);
    mds.spots["EUR.USD"] = S0 - bump;
    auto pv_dn           = engine.price(mds);
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
