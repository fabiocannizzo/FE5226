#pragma once

#include "Interfaces.h"
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace mc {

class MCEngine {
public:
    // Construct from already-generated pricers
    MCEngine(std::vector<pmcpricer_t> pricers,
             Date                     today       = 0,
             size_t                   n_paths     = 100000,
             uint64_t                 seed        = 42,
             bool                     antithetic  = true);

    // Construct directly from a vector of payoffs (with events and pricer config)
    MCEngine(const std::vector<ppayoff_t>& payoffs,
             const Events&                 events,
             const PricerConfig&           config      = PricerConfig{},
             Date                          today       = 0,
             size_t                        n_paths     = 100000,
             uint64_t                      seed        = 42,
             bool                          antithetic  = true);

    // Price all registered pricers against prepopulated MarketDataServer
    std::vector<double> price(const MarketDataServer& mkt) const;

    // Price all registered pricers against the polymorphic market dictionary
    std::vector<double> price(const std::map<std::string, IMarketObject*>& mkt) const;

    // Pricing function taking a vector of payoffs directly against MarketDataServer
    std::vector<double> price(const std::vector<ppayoff_t>& payoffs,
                              const Events&                 events,
                              const PricerConfig&           config,
                              const MarketDataServer&       mkt) const;

    const std::vector<Date>&        masterDates() const  { return m_dates; }
    const std::vector<std::string>& masterModels() const { return m_models; }

private:
    void buildMasterTimelineAndModels();

    std::vector<pmcpricer_t>         m_pricers;
    Date                             m_today;
    size_t                           m_n_paths;
    uint64_t                         m_seed;
    bool                             m_antithetic;
    mutable std::vector<Date>        m_dates;  // Merged & sorted master timeline
    mutable std::vector<std::string> m_models; // Unique & sorted required models
};

// Standalone pricing function taking a vector of payoffs
inline std::vector<double> price(const std::vector<ppayoff_t>& payoffs,
                                 const Events&                 events,
                                 const PricerConfig&           config,
                                 const MarketDataServer&       mkt,
                                 Date                          today      = 0,
                                 size_t                        n_paths    = 100000,
                                 uint64_t                      seed       = 42,
                                 bool                          antithetic = true) {
    MCEngine engine(payoffs, events, config, today, n_paths, seed, antithetic);
    return engine.price(mkt);
}

} // namespace mc
