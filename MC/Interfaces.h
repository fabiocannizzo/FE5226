#pragma once

#include "Matrix.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace mc {

// Base polymorphic interface
struct IObject {
    virtual ~IObject() = default;
};

// Date represented as integer day offset from valuation date (T_0 = 0)
using Date = int;

// ----------------------------------------------------------------------------
// Event hierarchy: polymorphic lifecycle events (market fixings, payments, etc.)
// ----------------------------------------------------------------------------
enum class EventTypeEnum {
    Fixing,
    Cashflow,
    Decision,
    Touch
};

struct IEvent : IObject {
    virtual ~IEvent() = default;
    virtual EventTypeEnum EventType() const = 0;
    virtual Date          date() const = 0;
};

class BaseEvent : public IEvent {
public:
    explicit BaseEvent(Date date) : m_date(date) {}

    Date date() const override { return m_date; }

protected:
    Date m_date = 0;
};

class FixingEvent : public BaseEvent {
public:
    FixingEvent(std::string asset, Date date, double value)
        : BaseEvent(date), m_asset(std::move(asset)), m_value(value) {}

    EventTypeEnum      EventType() const override { return EventTypeEnum::Fixing; }
    const std::string& asset() const { return m_asset; }
    double             value() const { return m_value; }

protected:
    std::string m_asset;
    double      m_value = 0.0;
};

template <typename T>
using EventPtr = std::shared_ptr<T>;

using Event  = std::shared_ptr<IEvent>;
using Events = std::vector<Event>;

inline bool lookupEvent(const Events&      events,
                        const std::string& asset,
                        Date               d,
                        double&            out) {
    for (const auto& ev : events) {
        if (ev && ev->date() == d && ev->EventType() == EventTypeEnum::Fixing) {
            const auto* fe = static_cast<const FixingEvent*>(ev.get());
            if (fe->asset() == asset) {
                out = fe->value();
                return true;
            }
        }
    }
    return false;
}

inline double getEvent(const Events&      events,
                       const std::string& asset,
                       Date               d) {
    double out = 0.0;
    MYASSERT(lookupEvent(events, asset, d, out),
             "Fixing event not found: " << asset << " at date " << d);
    return out;
}

inline std::shared_ptr<FixingEvent> makeFixingEvent(std::string asset, Date date, double value) {
    return std::make_shared<FixingEvent>(std::move(asset), date, value);
}

// Forward declarations
struct IModelMC;
using IModel = IModelMC;

// ----------------------------------------------------------------------------
// 1. Pricer Interface (IPricer / IMCPricer) & Configuration
// ----------------------------------------------------------------------------
struct PricerConfig {
    std::string                        default_model = "BS_CLOSED"; // e.g. "BS_CLOSED", "BS_EULER"
    std::map<std::string, std::string> model_overrides;             // asset -> model family
    std::string                        base_ccy      = "USD";       // reporting currency
};

struct IPricer : IObject {
    // 1. Future simulation dates required by this pricer (dates not yet observed)
    virtual std::vector<Date> requiredDates() const = 0;

    // 2. Concrete diffusion models needed (e.g. {"BS_CLOSED.EUR.USD"})
    virtual std::vector<std::string> modelNeeded() const = 0;

    // 3. Settlement date & currency for discounting
    virtual Date        settlementDate() const = 0;
    virtual std::string baseCurrency() const = 0;

    // 4. Pre-bind integer indices into the engine's master timeline, model list & model objects
    virtual void bindIndices(const std::vector<Date>&            master_dates,
                             const std::vector<std::string>&     master_models,
                             const std::vector<const IModelMC*>& models) = 0;

    // 5. Evaluate discounted PV on one simulated trajectory
    //    Stitches past events and reconstructed spot prices
    virtual double pricePath(const std::vector<Matrix>& model_paths,
                             double                     df_settle,
                             double                     fx_to_base) const = 0;
};

using IMCPricer = IPricer;
using pmcpricer_t = std::shared_ptr<IPricer>;
using ppricer_t = std::shared_ptr<IPricer>;

// ----------------------------------------------------------------------------
// 2. Payoff Interface (Derivative Contract Specification & Pricer Factory)
// ----------------------------------------------------------------------------
struct IPayoff : IObject {
    // Factory: generates dedicated pricer lifecycled against historical events and configured with model choices
    virtual std::shared_ptr<IPricer> getPricer(const Events&       events,
                                               const PricerConfig& config = PricerConfig{}) const = 0;

    // Currencies of cashflows generated by this trade
    virtual std::vector<std::string> getCashFlowCurrencies(const Events& events) const = 0;

    // Observation dates (market fixing / evaluation dates)
    virtual std::vector<Date> observationDates(const Events& events) const = 0;

    // Trade maturity date (final expiry / observation date)
    virtual Date maturityDate(const Events& events) const = 0;

    // Cashflow payment / settlement dates
    virtual std::vector<Date> settlementDates(const Events& events) const = 0;
};

using ppayoff_t = std::shared_ptr<const IPayoff>;

// ----------------------------------------------------------------------------
// 3. Market Objects, Models, Discount Curves & Correlations
// ----------------------------------------------------------------------------
struct IMarketObject : IObject {
    virtual std::string name() const = 0;
};

struct IModelMC : IMarketObject {
    // Number of driving standard Brownian motions (NFactors)
    virtual size_t numberOfBrownians() const = 0;

    // Number of diffused state variables per time step (NStates)
    virtual size_t numberOfStateVariables() const = 0;

    // Optional hook to insert intermediate discretization dates (e.g. monthly for Euler)
    virtual std::vector<Date> mergeDates(const std::vector<Date>& dates) const {
        return dates;
    }

    // Pre-compute time-step constants on the merged timeline
    virtual void configureTimeline(const std::vector<Date>& dates) = 0;

    // Diffuse state variables across all time steps:
    // rndNumbers: correlated standard normals [numberOfBrownians() x n_steps]
    // stateVariables: diffused state variables [numberOfStateVariables() x n_dates]
    virtual void diffuse(const ConstMatrixView& rndNumbers,
                         MatrixView&            stateVariables) const = 0;

    // Observable function: queries model object passing stateVariables to reconstruct observable (e.g. Spot)
    virtual double observable(size_t date_idx, const Matrix& stateVariables) const = 0;

    // Convenience alias for Spot observable
    virtual double spot(size_t date_idx, const Matrix& stateVariables) const {
        return observable(date_idx, stateVariables);
    }
};

// ----------------------------------------------------------------------------
// 4. Prepopulated Market Data Server
// ----------------------------------------------------------------------------
struct ICurve : IMarketObject {
    virtual double df(Date d) const = 0;
};

struct ICorrelation : IMarketObject {
    virtual double value() const = 0;
};

struct MarketDataServer {
    std::map<std::string, double>  spots;
    std::map<std::string, double>  vols;
    std::map<std::string, double>  drifts;
    std::map<std::string, double>  correlations;
    std::map<std::string, ICurve*> discount_curves;

    double getSpot(const std::string& asset) const {
        auto it = spots.find(asset);
        MYASSERT(it != spots.end(), "Spot not found in MarketDataServer: " << asset);
        return it->second;
    }
    double getVol(const std::string& asset) const {
        auto it = vols.find(asset);
        MYASSERT(it != vols.end(), "Vol not found in MarketDataServer: " << asset);
        return it->second;
    }
    double getDrift(const std::string& asset) const {
        auto it = drifts.find(asset);
        MYASSERT(it != drifts.end(), "Drift not found in MarketDataServer: " << asset);
        return it->second;
    }
    double getCorrelation(const std::string& rfA, const std::string& rfB) const {
        if (rfA == rfB)
            return 1.0;
        std::string key = (rfA < rfB) ? "CORR." + rfA + ":" + rfB : "CORR." + rfB + ":" + rfA;
        auto        it  = correlations.find(key);
        if (it != correlations.end())
            return it->second;
        return 0.0;
    }
    ICurve* getCurve(const std::string& ccy) const {
        auto it = discount_curves.find(ccy);
        MYASSERT(it != discount_curves.end(), "Discount curve not found: " << ccy);
        return it->second;
    }
};

// Canonical dictionary key for pairwise correlation between risk factors rfA and rfB
// e.g. corrKey("BS.EUR.USD.0", "BS.GBP.USD.0") -> "CORR.BS.EUR.USD.0:BS.GBP.USD.0"
inline std::string corrKey(const std::string& rfA, const std::string& rfB) {
    if (rfA < rfB)
        return "CORR." + rfA + ":" + rfB;
    return "CORR." + rfB + ":" + rfA;
}

inline double getCorrelation(const std::map<std::string, IMarketObject*>& mkt,
                             const std::string&                           rfA,
                             const std::string&                           rfB) {
    if (rfA == rfB)
        return 1.0; // Unit diagonal
    auto it = mkt.find(corrKey(rfA, rfB));
    if (it == mkt.end())
        return 0.0; // Uncorrelated by default if omitted
    auto* c = dynamic_cast<ICorrelation*>(it->second);
    MYASSERT(c != nullptr, "Market object is not an ICorrelation: " << corrKey(rfA, rfB));
    return c->value();
}

} // namespace mc
