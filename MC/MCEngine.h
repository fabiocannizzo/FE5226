#pragma once

#include "Interfaces.h"
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace mc {

class MCEngine {
public:
    MCEngine(std::vector<pmcpricer_t> pricers,
             Date                     today,
             size_t                   n_paths,
             uint64_t                 seed = 42,
             bool                     antithetic = true);

    // Price all registered pricers against prepopulated MarketDataServer
    std::vector<double> price(const MarketDataServer& mkt) const;

    // Price all registered pricers against the polymorphic market dictionary
    std::vector<double> price(const std::map<std::string, IMarketObject*>& mkt) const;

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

} // namespace mc
