#pragma once

#include "Interfaces.h"
#include <algorithm>
#include <memory>
#include <vector>

namespace mc {

// ----------------------------------------------------------------------------
// Composite pricer: a weighted list of component pricers. It owns no cash flow; the engine
// flattens it into its leaf pricers, simulates once, and reports total and breakdown.
// ----------------------------------------------------------------------------
// Composite pattern: the weighted list of the component pricers. It has no cash flow of its own
// (fixingDate / settlementDate / pricePath are the throwing defaults of IPricer); the engine
// asks components() and componentWeights() and flattens the tree down to its leaves.
// The date / model / binding members just forward to the components.
class SuperPricer : public IPricer {
public:
    SuperPricer(std::vector<pmcpricer_t> components, std::vector<double> weights)
        : m_components(std::move(components)), m_weights(std::move(weights)) {}

    std::vector<Date> requiredDates() const override {
        std::vector<Date> out;
        for (const auto& c : m_components) {
            auto d = c->requiredDates();
            out.insert(out.end(), d.begin(), d.end());
        }
        return out;
    }
    std::vector<std::string> modelNeeded() const override {
        std::vector<std::string> out;
        for (const auto& c : m_components) {
            auto m = c->modelNeeded();
            out.insert(out.end(), m.begin(), m.end());
        }
        return out;
    }
    void bindIndices(const std::vector<Date>& dates, const std::vector<std::string>& models_names,
                     const std::vector<const IModelMC*>& models) override {
        for (auto& c : m_components)
            c->bindIndices(dates, models_names, models);
    }
    // This is what makes the engine treat a portfolio like a trade.
    std::vector<pmcpricer_t> components() const override { return m_components; }
    std::vector<double>      componentWeights() const override { return m_weights; }

private:
    std::vector<pmcpricer_t> m_components;
    std::vector<double>      m_weights;
};

// ----------------------------------------------------------------------------
// SuperPayoff: a portfolio that is itself a payoff. It carries its own PricerConfig, used
// consistently for all components (the config passed to getPricer is ignored), so the whole
// portfolio is priced in a single MC with one model per asset.
// ----------------------------------------------------------------------------
// A portfolio expressed as ONE payoff. Differences from pricing the trades one by one:
//   - it owns a PricerConfig, applied consistently to every component (the config passed to
//     getPricer is ignored), so all components use the same model for a given asset;
//   - the whole book is simulated once, on the same paths;
//   - the result carries the total and the breakdown by component, with standard errors.
class SuperPayoff : public IPayoff {
public:
    // weights: position size of each component (default 1 = long one unit)
    SuperPayoff(std::vector<ppayoff_t> components, PricerConfig config, std::vector<double> weights = {})
        : m_components(std::move(components)), m_config(std::move(config)), m_weights(std::move(weights)) {
        MYASSERT(!m_components.empty(), "SuperPayoff requires at least one component");
        if (m_weights.empty())
            m_weights.assign(m_components.size(), 1.0);
        MYASSERT(m_weights.size() == m_components.size(), "SuperPayoff: one weight per component required");
    }

    const PricerConfig& config() const { return m_config; }
    size_t              size() const { return m_components.size(); }

    // Builds the component pricers with the super payoff's own config.
    std::shared_ptr<IPricer> getPricer(const Events& events, const PricerConfig& = PricerConfig{}) const override {
        std::vector<pmcpricer_t> pricers;
        for (const auto& p : m_components)
            pricers.push_back(p->getPricer(events, m_config));
        return std::make_shared<SuperPricer>(std::move(pricers), m_weights);
    }

    std::vector<std::string> getCashFlowCurrencies(const Events& events) const override {
        std::vector<std::string> out;
        for (const auto& p : m_components) {
            auto v = p->getCashFlowCurrencies(events);
            out.insert(out.end(), v.begin(), v.end());
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }
    std::vector<Date> observationDates(const Events& events) const override {
        return unionOf([&](const IPayoff& p) { return p.observationDates(events); });
    }
    std::vector<Date> settlementDates(const Events& events) const override {
        return unionOf([&](const IPayoff& p) { return p.settlementDates(events); });
    }
    Date maturityDate(const Events& events) const override {
        Date m = m_components.front()->maturityDate(events);
        for (const auto& p : m_components)
            m = std::max(m, p->maturityDate(events));
        return m;
    }

private:
    // Sorted union of a per-component list of dates.
    template <typename F>
    std::vector<Date> unionOf(F f) const {
        std::vector<Date> out;
        for (const auto& p : m_components) {
            auto v = f(*p);
            out.insert(out.end(), v.begin(), v.end());
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    std::vector<ppayoff_t> m_components;
    PricerConfig           m_config;
    std::vector<double>    m_weights;
};

} // namespace mc
