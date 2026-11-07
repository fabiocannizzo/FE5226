#pragma once

#include "Interfaces.h"
#include "Market.h"
#include "SuperPayoff.h"
#include <cstdint>
#include <string>
#include <vector>

namespace mc {

// Monte Carlo estimate of a PV in the numeraire currency (USD) with its standard error.
// For a composite (SuperPayoff) 'components' holds the weighted breakdown of the leaf trades.
// Result of a valuation: the Monte Carlo estimate (USD, as of today) and its standard error, so
// the user knows how much to trust the digits. For a composite (SuperPayoff) `components` holds
// the weighted contribution of each leaf trade, in the order of the portfolio.
struct PriceResult {
    double                   pv      = 0.0;
    double                   std_err = 0.0;
    std::vector<PriceResult> components;
};

// One simulation for a set of pricers. Composite pricers are flattened into their leaves, all
// leaves share the simulation, and the result has one entry per pricer passed in.
// MC settings (paths, seed, antithetic) and the per-asset model choice come from PricerConfig;
// curves, spots, vols, correlations and the models themselves come from the Market.
// The Monte Carlo engine. It knows nothing about trade types or model types: it only talks to
// IPricer, IModelMC and Market. What it does, in order (see MCEngine.cpp):
//   1. flatten composites into leaf pricers;
//   2. collect the dates and models the leaves need; build the models (Market) and let them
//      add their own dates; configure and bind everything;
//   3. compute the path-independent discount factor of each cash flow;
//   4. assemble the correlation matrix, Cholesky-factor it;
//   5. loop on paths: draw, correlate, advance every model ONCE (step by step), price every leaf on the same path;
//   6. average, with standard error; recombine leaves into the requested results.
// All leaves share the simulation, so pricing N trades costs one simulation, not N.
class MCEngine {
public:
    // `pricers`: one result will be returned per entry (an entry may itself be a composite).
    // `config`: n_paths, seed, antithetic and the per-asset model parameters used by Market.
    MCEngine(std::vector<pmcpricer_t> pricers, PricerConfig config = PricerConfig{});

    // Run the simulation against `mkt`. Const and re-entrant: all work state is local, and the RNG is
    // re-seeded on every call, so the same inputs give the same output (and bumping one market input
    // re-uses the same random numbers: Common Random Numbers).
    std::vector<PriceResult> price(const Market& mkt) const;

private:
    // A leaf pricer (index into m_leaf_pricers) and the weight with which it enters a requested result.
    struct Term {
        size_t leaf;
        double weight;
    };
    // Recursively expand a (possibly composite) pricer into weighted leaves.
    void flatten(const pmcpricer_t& p, double weight, std::vector<Term>& terms);

    PricerConfig                  m_config;
    std::vector<pmcpricer_t>      m_leaf_pricers;
    // m_terms[r] = leaves (with weights) whose sum is the r-th requested result.
    std::vector<std::vector<Term>> m_terms; // per input pricer: weighted leaves
    std::vector<bool>             m_composite;
};

// ----------------------------------------------------------------------------
// Pricing modes
// ----------------------------------------------------------------------------

// Mode 1: one pricer, hence one MC simulation, per payoff
// PRICING MODE 1: loop over the payoffs, get a pricer from each, and run ONE Monte Carlo per trade.
// Pro: each trade can have its own model choice, and trades are fully independent.
// Con: n trades cost n simulations, and trades do not share paths.
inline std::vector<PriceResult> priceEach(const std::vector<ppayoff_t>& payoffs,
                                          const Events&                 events,
                                          const Market&                 mkt,
                                          const PricerConfig&           config) {
    std::vector<PriceResult> out;
    out.reserve(payoffs.size());
    for (const auto& payoff : payoffs) {
        MCEngine engine({ payoff->getPricer(events, config) }, config);
        out.push_back(engine.price(mkt).front());
    }
    return out;
}

// Mode 2: the whole portfolio is one payoff, simulated once with its own config;
// the result has the total and the breakdown by component
// PRICING MODE 2: the whole portfolio is a SuperPayoff (itself a payoff) priced in one Monte Carlo
// with its own, consistent config. Returns the total and the breakdown by component.
inline PriceResult priceSuper(const SuperPayoff& super, const Events& events, const Market& mkt) {
    MCEngine engine({ super.getPricer(events) }, super.config());
    return engine.price(mkt).front();
}

} // namespace mc
