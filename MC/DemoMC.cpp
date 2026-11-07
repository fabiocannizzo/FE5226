#include "MCEngine.h"
#include "Market.h"
#include "Models.h"
#include "Payoffs.h"
#include "Pricers.h"
#include "SuperPayoff.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace mc;

namespace {

// Standard normal CDF via std::erfc
double normCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// Analytical Black-Scholes European Call / Put formula for benchmark validation.
// 'df' is the discount factor to the payment date, 'mu' the drift of the underlying.
double blackScholesPrice(double S0,
                         double K,
                         double df,
                         double mu,
                         double sigma,
                         double T,
                         int    omega,
                         double notional = 1.0) {
    double F      = S0 * std::exp(mu * T);
    double vol_sq = sigma * std::sqrt(T);
    double d1     = (std::log(F / K) + 0.5 * sigma * sigma * T) / vol_sq;
    double d2     = d1 - vol_sq;
    return notional * df * (omega * F * normCDF(omega * d1) - omega * K * normCDF(omega * d2));
}

// The demo is also a regression test: every check() that fails increments this counter and
// main() returns non-zero, so it can be run from a script or a CI job.
int failures = 0;

// Statistical check: |mc - exact| within 4 standard errors
// Statistical check: a Monte Carlo value is accepted if it lies within 4 standard errors of the
// exact value (probability of a false alarm ~ 6e-5). Never compare MC results with a fixed
// absolute tolerance: the right tolerance is the standard error.
void check(const std::string& what, double mc, double se, double exact) {
    bool ok = std::abs(mc - exact) <= 4.0 * se + 1e-9;
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << std::left << std::setw(40) << what << std::right
              << std::setw(14) << mc << " +/- " << std::setw(9) << se << "  exact " << std::setw(14) << exact << "\n";
    if (!ok)
        ++failures;
}

// Boolean check with the same reporting.
void check(const std::string& what, bool ok) {
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << what << "\n";
    if (!ok)
        ++failures;
}

} // namespace

// Walk-through, in the order a user would write it:
//   1 number of paths and seed (in PricerConfig)   2 time (Date = int, today = 0)
//   3 events (past fixings)                        4 payoffs (the trades)
//   5 portfolio                                    6 PricerConfig (model choice per asset)
//   7 market data                                  8 the Market object
// followed by the two pricing modes and the validation against closed forms.
int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "====================================================================\n";
    std::cout << "  Multi-Model MC Engine: Market, per-asset config, two pricing modes\n";
    std::cout << "  Settlement after fixing, FX conversion in the USD numeraire (FE5226)\n";
    std::cout << "====================================================================\n\n";

    // Date is a plain int: days from today (t0 = 0). A real system would use a Date class.
    const Date today = 0;

    // ------------------------------------------------------------------------
    // 1. Market: curves, spots, vols, correlations. FX drifts are implied by the curves.
    // ------------------------------------------------------------------------
    // Flat annual rates and vols of the two quoted assets. Note that NO drift is specified: the
    // Market derives it from the curves (mu = r_USD - r_foreign), here mu_eur and mu_gbp, which are
    // also needed below only to compute the closed-form benchmarks.
    const double r_usd = 0.04, r_eur = 0.025, r_gbp = 0.03;
    const double vol_eur = 0.12, vol_gbp = 0.14, rho = 0.50;
    const double mu_eur = r_usd - r_eur, mu_gbp = r_usd - r_gbp;

    // Market quotes only CCY.USD assets; USD is the numeraire. Any other pair (EUR.GBP below) is
    // derived from them by the pricers.
    Market market(today);
    market.setFlatCurve("USD", r_usd);
    market.setFlatCurve("EUR", r_eur);
    market.setFlatCurve("GBP", r_gbp);
    market.setSpot("EUR.USD", 1.10);
    market.setSpot("GBP.USD", 1.28);
    market.setVol("EUR.USD", vol_eur);
    market.setVol("GBP.USD", vol_gbp);
    // Correlation of the two assets; it also determines the volatility of the EUR.GBP cross.
    market.setCorrelation("EUR.USD", "GBP.USD", rho);

    // ------------------------------------------------------------------------
    // 2. Historical fixings
    // ------------------------------------------------------------------------
    // Fixings already observed on the valuation date. Only GBP.USD has any: the Asian option below
    // is "seasoned", i.e. it started in the past, and everything else is fully in the future.
    Events events = {
        makeFixingEvent("GBP.USD", -90, 1.25),
        makeFixingEvent("GBP.USD",   0, 1.28)
    };

    // ------------------------------------------------------------------------
    // 3. Payoffs. Settlement dates are on or after the fixing dates.
    // ------------------------------------------------------------------------
    // Common notional of all trades, to make the numbers comparable.
    const double N = 1000000.0;
    // Trade 0: at-the-money EUR.USD call, 6 months, paid at expiry: the reference trade with a
    // Black-Scholes closed form.
    auto vanilla_call = std::make_shared<PayoffVanilla>(
        "EUR.USD", N, 1.10, +1, /*expiry*/ 180, /*settle*/ 180, "USD");
    auto asian_call = std::make_shared<PayoffAsian>(
        "GBP.USD", N, 1.26, +1, std::vector<Date>{ -90, 0, 90, 180 }, /*settle*/ 180, "USD");
    auto worst_of_call = std::make_shared<PayoffWorstOfCall>(
        "EUR.USD", "GBP.USD", 1.10, 1.28, N, 1.00, /*expiry*/ 180, /*settle*/ 180, "USD");
    // Trade 3: same option, but cash settles 30 days AFTER the fixing: tests the forward discount
    // P(Tf,Ts); closed form = Black-Scholes discounted to the settlement date.
    auto vanilla_late = std::make_shared<PayoffVanilla>(          // pays 30 days after expiry
        "EUR.USD", N, 1.10, +1, 180, /*settle*/ 210, "USD");
    // Trade 4: payoff N * X(T) paid in GBP, where X is GBP.USD (strike 0, so no optionality). It
    // exercises the foreign-currency path: conversion with the contingent FX spot at the fixing and the
    // GBP curve between fixing and payment. It has a closed form because E[X^2] is lognormal.
    auto gbp_paid = std::make_shared<PayoffVanilla>(              // GBP.USD (K=0) paid in GBP
        "GBP.USD", N, 0.0, +1, 180, 210, "GBP");
    // Trade 5: option on the CROSS EUR.GBP. No EUR.GBP model exists: its spot is the ratio of the
    // simulated EUR.USD and GBP.USD, hence it depends on their correlation.
    auto cross_call = std::make_shared<PayoffVanilla>(            // EUR.GBP, simulated via USD pairs
        "EUR.GBP", N, 0.86, +1, 180, 210, "USD");

    std::vector<ppayoff_t> portfolio = { vanilla_call, asian_call, worst_of_call, vanilla_late, gbp_paid, cross_call };
    const std::vector<std::string> names = {
        "EUR vanilla call, settle = expiry", "GBP Asian call (seasoned)", "EUR/GBP worst-of call",
        "EUR vanilla call, settle +30d",     "GBP.USD K=0 paid in GBP",   "EUR.GBP cross call, settle +30d" };

    std::cout << "[IPayoff query methods of the portfolio as a SuperPayoff]\n";
    // Config 1: closed-form Black-Scholes everywhere (the default). 256000 paths -> standard errors
    // of order 0.2% of the price for the vanilla.
    PricerConfig config_closed;                       // BS_CLOSED for all assets
    config_closed.n_paths = 256000;
    // Config 2: same settings, but Euler discretization for every asset.
    PricerConfig config_euler = config_closed;
    config_euler.default_model = "BS_EULER";
    auto super_closed = std::make_shared<SuperPayoff>(portfolio, config_closed);
    std::cout << "  settlementDates() : ";
    for (Date d : super_closed->settlementDates(events)) std::cout << d << " ";
    std::cout << "\n  maturityDate()    : " << super_closed->maturityDate(events) << "\n  cashflow ccys     : ";
    for (const auto& c : super_closed->getCashFlowCurrencies(events)) std::cout << c << " ";
    std::cout << "\n\n";

    // ------------------------------------------------------------------------
    // 4. Closed form benchmarks (all in USD, BS_CLOSED, as of today)
    // ------------------------------------------------------------------------
    const double T = 180.0 / 365.0;
    std::vector<double> exact(portfolio.size(), std::nan(""));
    // Closed-form benchmarks, only for trades that have one. All use the same forward drift and the
    // same discount factors as the engine, so a difference is Monte Carlo (or discretization) error.
    exact[0] = blackScholesPrice(1.10, 1.10, market.df("USD", 180), mu_eur, vol_eur, T, +1, N);
    exact[3] = blackScholesPrice(1.10, 1.10, market.df("USD", 210), mu_eur, vol_eur, T, +1, N);
    // Payoff N*X(T) in GBP, X = GBP.USD: PV = Pd(0,Tf) * Pg(Tf,Ts) * N * E[X(T)^2]
    exact[4] = N * market.df("USD", 180) * market.df("GBP", 210) / market.df("GBP", 180)
               * 1.28 * 1.28 * std::exp((2.0 * mu_gbp + vol_gbp * vol_gbp) * T);
    // Cross EUR.GBP = EUR.USD / GBP.USD is lognormal under the USD measure
    {
        double mu_c  = mu_eur - mu_gbp + vol_gbp * vol_gbp - rho * vol_eur * vol_gbp;
        double vol_c = std::sqrt(vol_eur * vol_eur + vol_gbp * vol_gbp - 2.0 * rho * vol_eur * vol_gbp);
        exact[5]     = blackScholesPrice(1.10 / 1.28, 0.86, market.df("USD", 210), mu_c, vol_c, T, +1, N);
    }

    // ------------------------------------------------------------------------
    // 5. Mode 1: one pricer and one MC per payoff
    // ------------------------------------------------------------------------
    // MODE 1: one simulation per trade. Each result has its own standard error.
    std::cout << "[Mode 1: one MC per payoff, BS_CLOSED, N = " << config_closed.n_paths << " paths]\n";
    auto res1 = priceEach(portfolio, events, market, config_closed);
    for (size_t i = 0; i < portfolio.size(); ++i) {
        if (std::isnan(exact[i]))
            std::cout << "  " << std::left << std::setw(45) << names[i] << std::right << std::setw(14)
                      << res1[i].pv << " +/- " << std::setw(9) << res1[i].std_err << "\n";
        else
            check(names[i], res1[i].pv, res1[i].std_err, exact[i]);
    }

    // ------------------------------------------------------------------------
    // 6. Mode 2: the portfolio is a SuperPayoff, one MC, total and breakdown
    // ------------------------------------------------------------------------
    // MODE 2: the same trades as ONE payoff, one simulation: trades see the same paths. The total and
    // the breakdown come from that single run; each component must agree with mode 1 within the
    // combined statistical error (the random numbers differ because the set of models differs).
    std::cout << "\n[Mode 2: SuperPayoff, one MC for the portfolio]\n";
    PriceResult res2 = priceSuper(*super_closed, events, market);
    double sum_parts = 0.0;
    for (size_t i = 0; i < portfolio.size(); ++i) {
        const PriceResult& c = res2.components[i];
        sum_parts += c.pv;
        std::cout << "  " << std::left << std::setw(45) << names[i] << std::right << std::setw(14) << c.pv
                  << " +/- " << std::setw(9) << c.std_err << "\n";
        if (!std::isnan(exact[i]))
            check(names[i] + " vs closed form", c.pv, c.std_err, exact[i]);
        check(names[i] + " mode 1 vs mode 2", std::abs(c.pv - res1[i].pv) <=
                  4.0 * std::hypot(c.std_err, res1[i].std_err));
    }
    std::cout << "  " << std::left << std::setw(45) << "TOTAL" << std::right << std::setw(14) << res2.pv
              << " +/- " << std::setw(9) << res2.std_err << "\n";
    check("total equals sum of components", std::abs(res2.pv - sum_parts) < 1e-6 * std::abs(sum_parts));

    // ------------------------------------------------------------------------
    // 7. Per-asset config. Mode 1: EUR.USD trades closed form, the EUR.GBP cross with Euler
    //    (its USD legs are simulated with Euler for that trade only).
    // ------------------------------------------------------------------------
    std::cout << "\n[Per-asset config, mode 1: EUR.USD closed form, EUR.GBP Euler (30d steps)]\n";
    // Per-asset configuration: closed form for EUR.USD, Euler for the cross EUR.GBP (and therefore for
    // the two USD legs of that trade). Fine in mode 1, where each trade has its own simulation.
    PricerConfig config_mixed = config_closed;
    config_mixed.setModel("EUR.USD", "BS_CLOSED");
    config_mixed.setModel("EUR.GBP", "BS_EULER");
    auto res_mixed = priceEach({ vanilla_late, cross_call }, events, market, config_mixed);
    std::cout << "  EUR vanilla, closed form : " << std::setw(14) << res_mixed[0].pv << " +/- " << res_mixed[0].std_err << "\n";
    std::cout << "  EUR.GBP cross, Euler     : " << std::setw(14) << res_mixed[1].pv << " +/- " << res_mixed[1].std_err
              << "  (closed form: " << exact[5] << ")\n";

    // The same config in a single simulation would simulate EUR.USD twice: rejected.
    // In mode 2 the same config would need two models for EUR.USD in ONE simulation (closed form for
    // the vanilla, Euler for the cross): the engine must refuse instead of silently picking one.
    bool rejected = false;
    try {
        priceSuper(SuperPayoff({ vanilla_late, cross_call }, config_mixed), events, market);
    } catch (const std::invalid_argument& e) {
        rejected = true;
        std::cout << "  Mode 2 with this config is rejected: " << e.what() << "\n";
    }
    check("conflicting per-asset models in one MC rejected", rejected);

    std::cout << "\n[Mode 2 with BS_EULER for all assets]\n";
    PriceResult res_euler = priceSuper(SuperPayoff(portfolio, config_euler), events, market);
    for (size_t i = 0; i < portfolio.size(); ++i)
        std::cout << "  " << std::left << std::setw(45) << names[i] << std::right << std::setw(14)
                  << res_euler.components[i].pv << " +/- " << std::setw(9) << res_euler.components[i].std_err
                  << "  (closed: " << res2.components[i].pv << ")\n";
    std::cout << "  " << std::left << std::setw(45) << "TOTAL" << std::right << std::setw(14) << res_euler.pv << "\n";

    // ------------------------------------------------------------------------
    // 8. Common Random Numbers (CRN) finite-difference delta on EUR.USD
    // ------------------------------------------------------------------------
    // Delta by central finite difference with Common Random Numbers: both bumped runs use the same
    // seed, hence the same random numbers, so the sampling noise cancels in the difference. With
    // independent noise the error would be of order (price noise) / (2 * bump) and unusable.
    const double S0   = 1.10;
    const double bump = 1e-4;

    market.setSpot("EUR.USD", S0 + bump);
    auto pv_up = priceEach({ vanilla_call }, events, market, config_closed);
    market.setSpot("EUR.USD", S0 - bump);
    auto pv_dn = priceEach({ vanilla_call }, events, market, config_closed);
    market.setSpot("EUR.USD", S0);

    double mc_delta = (pv_up[0].pv - pv_dn[0].pv) / (2.0 * bump);
    double d1       = (mu_eur + 0.5 * vol_eur * vol_eur) * T / (vol_eur * std::sqrt(T));
    double bs_del   = N * std::exp((mu_eur - r_usd) * T) * normCDF(d1);

    std::cout << "\n[CRN finite-difference delta (dS = " << bump << ")]\n";
    std::cout << "  MC delta " << mc_delta << ", analytic " << bs_del << "\n";
    check("CRN delta within 0.1% of analytic", std::abs(mc_delta - bs_del) < 1e-3 * bs_del);

    std::cout << "\n" << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << "\n";
    return failures == 0 ? 0 : 1;
}
