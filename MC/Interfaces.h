#pragma once

#include "Matrix.h"
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>
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
// Kinds of events. An event is a FACT that has already happened by the valuation date
// (a market fixing, a payment made, an exercise decision, a barrier touch). Events are what
// distinguishes a seasoned trade from a new one: a pricer replaces the random, simulated value
// of a past observation by the recorded one.
// Only Fixing is implemented in this example; the others are placeholders showing how the
// hierarchy would be extended.
enum class EventTypeEnum {
    Fixing,    // e.g. stock closeing price
    Cashflow,  // e.g. amount of USD paid
    Decision,  // e.g. for a vanilla option
    Touch      // e.g. if a certain price has been touched
};

// Interface of every event. EventType() allows safe down-casting (see lookupEvent) without
// RTTI; date() is the day the event happened, as an offset from the valuation date (< 0 or = 0).
struct IEvent : IObject {
    virtual ~IEvent() = default;
    virtual EventTypeEnum EventType() const = 0;
    virtual Date          date() const = 0;
};

// Convenience base class: stores the date once, so concrete events only add their own payload.
class BaseEvent : public IEvent {
public:
    explicit BaseEvent(Date date) : m_date(date) {}

    Date date() const override { return m_date; }

protected:
    Date m_date = 0;
};

// A recorded market observation: "asset X was worth `value` on `date`".
// The asset name follows the library-wide convention BASE.QUOTE (e.g. "GBP.USD" = USD per 1 GBP).
// The very same string is the key of spots, vols, correlations and model tags, so a fixing is
// matched to its asset by plain string equality.
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

// Search the events for a Fixing of `asset` on date `d`.
// Returns true and sets `out` if found, false otherwise (the caller then knows the observation
// lies in the future and must be simulated).
// Implementation: linear scan on (date, type, asset). Adequate for the handful of fixings of a
// demo; a production system would index events by (asset, date) in a map.
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

// Like lookupEvent, but a missing fixing is an error (throws via MYASSERT).
inline double getEvent(const Events&      events,
                       const std::string& asset,
                       Date               d) {
    double out = 0.0;
    MYASSERT(lookupEvent(events, asset, d, out),
             "Fixing event not found: " << asset << " at date " << d);
    return out;
}

// Factory for fixings; returns shared_ptr<FixingEvent>, which converts implicitly to Event
// (shared_ptr<IEvent>) so it can be stored directly in an Events vector.
inline std::shared_ptr<FixingEvent> makeFixingEvent(std::string asset, Date date, double value) {
    return std::make_shared<FixingEvent>(std::move(asset), date, value);
}

// Forward declarations
struct IModelMC;
using IModel = IModelMC;

// ----------------------------------------------------------------------------
// 1. Configuration & Pricer Interface
// ----------------------------------------------------------------------------

// Numeraire currency of every simulation. For simplicity it is always USD, and every
// FX asset is simulated as a "<CCY>.USD" pair (crosses are derived from USD pairs).
// A more sophisticated algorithm could choose the most natural numeraire per payoff
// (particularly for multi-asset payoffs); this is deliberately not done here.
inline const std::string NUMERAIRE_CCY = "USD";

// Per-asset simulation settings. Which stochastic model simulates an asset, and with which
// numerical parameters, is NOT a property of a trade or of the market: it is a numerical choice,
// hence it lives in PricerConfig. Assets without an entry use PricerConfig::default_model.
// Model choice and numerical parameters for one asset
struct AssetConfig {
    std::string family;                  // "BS_CLOSED" or "BS_EULER"; empty = use default_model
    int         euler_max_step = 30;     // BS_EULER: max step size in days
};

// Everything that controls HOW the simulation is run (as opposed to WHAT is priced = payoffs,
// or the state of the world = Market):
//   - which model simulates each asset (default + per-asset overrides),
//   - the Monte Carlo settings: number of paths, random seed, antithetic variates.
// The seed makes a run reproducible; re-using it after bumping a market input gives Common Random
// Numbers, which is what makes finite-difference Greeks stable.
struct PricerConfig {
    std::string                        default_model = "BS_CLOSED";
    std::map<std::string, AssetConfig> assets;      // per asset or per pair, e.g. "EUR.USD", "EUR.GBP"
    size_t                             n_paths    = 100000;
    uint64_t                           seed       = 42;
    bool                               antithetic = true;

    void setModel(const std::string& asset, const std::string& family) { assets[asset].family = family; }

    // Model family used to simulate the USD pair 'leg' when requested on behalf of 'pair'.
    // The requested pair's own entry (e.g. cross EUR.GBP) takes precedence over the leg's
    // entry (EUR.USD), which takes precedence over the default.
    std::string family(const std::string& leg, const std::string& pair = "") const {
        for (const std::string* key : { &pair, &leg }) {
            auto it = assets.find(*key);
            if (it != assets.end() && !it->second.family.empty())
                return it->second.family;
        }
        return default_model;
    }

    // Numerical parameters (e.g. the Euler step size) of a simulated USD asset
    AssetConfig assetConfig(const std::string& asset) const {
        auto it = assets.find(asset);
        return it != assets.end() ? it->second : AssetConfig{};
    }
};

struct IPricer;
using pmcpricer_t = std::shared_ptr<IPricer>;

// A pricer is either a leaf (one cash flow) or a composite (a weighted list of components).
// Pricer interface. A pricer is created by IPayoff::getPricer() and has this life cycle:
//   1. CONSTRUCTION  - lifecycles the contract against the Events: past fixings are frozen,
//                      only future observations are kept as "required dates".
//   2. DECLARATION   - the engine asks requiredDates() and modelNeeded() of every pricer and
//                      merges the answers into one master timeline and one list of models.
//   3. BINDING       - bindIndices() turns dates / model tags into integer positions once, so the
//                      path loop never searches or compares strings.
//   4. EVALUATION    - pricePath() values the contract on one simulated path, many times.
// Two kinds of pricer implement it:
//   - LEAF: one cash flow (fixing date, settlement date, payment currency, pricePath).
//   - COMPOSITE: a weighted list of pricers (components / componentWeights), e.g. a portfolio.
// The members of the other kind throw std::logic_error if called.
struct IPricer : IObject {
    // 1. Future simulation dates required by this pricer (dates not yet observed)
    virtual std::vector<Date> requiredDates() const = 0;

    // 2. Concrete diffusion models needed (e.g. {"BS_CLOSED.EUR.USD"}), including the FX model
    //    needed to convert a foreign-currency cash flow into the numeraire
    virtual std::vector<std::string> modelNeeded() const = 0;

    // 3. Pre-bind integer indices into the engine's master timeline, model list & model objects
    virtual void bindIndices(const std::vector<Date>&            master_dates,
                             const std::vector<std::string>&     master_models,
                             const std::vector<const IModelMC*>& models) = 0;

    // 4. Cash flow description (leaf pricers only). The payoff is determined on the fixing date
    //    and paid on the settlement date (>= fixing date) in the payment currency.
    virtual Date        fixingDate() const      { notACashflow(); }
    virtual Date        settlementDate() const  { notACashflow(); }
    virtual std::string paymentCurrency() const { notACashflow(); }

    // 5. Value on one simulated trajectory (leaf pricers only), in numeraire currency, as of today.
    //    'scale' is the path-independent discount factor computed by the engine from the market:
    //       scale = P_num(0,Tf) * P_ccy(Tf,Ts)
    //    i.e. discounting to today from the fixing date, and the conditional expectation on the
    //    fixing date of the (deterministic-rate) discount to the settlement date. The pricer
    //    multiplies by the contingent FX spot ccy.USD(Tf) if the payment currency is not USD.
    virtual double pricePath(const std::vector<Matrix>& /*model_paths*/, double /*scale*/) const {
        notACashflow();
    }

    // 6. Composite pricers (e.g. a portfolio) return their components and weights
    virtual std::vector<pmcpricer_t> components() const { return {}; }
    virtual std::vector<double>      componentWeights() const { return {}; }

private:
    // Called by the default implementations of the leaf-only members when invoked on a composite.
    [[noreturn]] static void notACashflow() {
        throw std::logic_error("IPricer: composite pricer has no single cash flow");
    }
};

using IMCPricer = IPricer;
using ppricer_t = std::shared_ptr<IPricer>;

// ----------------------------------------------------------------------------
// 2. Payoff Interface (Derivative Contract Specification & Pricer Factory)
// ----------------------------------------------------------------------------
// Payoff interface = a contract SPECIFICATION plus a pricer factory. It holds only the terms of
// the trade (underlying, strike, dates, currency) and is immutable, so one payoff object can be
// shared and priced under many markets, configs and event sets. The query methods describe the
// trade without pricing it (useful for reporting and for scheduling market data).
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

// Payoffs are shared as pointers-to-const: pricing never modifies a contract.
using ppayoff_t = std::shared_ptr<const IPayoff>;

// ----------------------------------------------------------------------------
// 3. Market Objects, Models, Discount Curves & Correlations
// ----------------------------------------------------------------------------
// Anything the market can hold and name: models, curves, correlations.
struct IMarketObject : IObject {
    virtual std::string name() const = 0;
};

// A stochastic model of ONE market asset (or of one family of related random factors), able to
// simulate its own state along the master timeline. It is the only place where the dynamics live.
// Division of labour:
//   - the ENGINE draws correlated normals and, step after step, gives each model the normals of
//     that step (one value per Brownian of the model);
//   - the MODEL turns them, together with its previous states, into the new states
//     (generateStates) and reconstructs observables from the stored states (observable);
//   - the PRICER reads observables, e.g. the spot, and never sees the state layout.
// Dimensions (one path; a "step" goes from date k to date k+1):
//   - rndNumbers of a step      : numberOfBrownians() values (correlated standard normals)
//   - previousStates / states   : numberOfStateVariables() values each
//   - state history kept by the engine for the pricers: a Matrix
//     [numberOfStateVariables() x n_dates], column 0 = today, n_dates = n_steps + 1.
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

    // States at the first date (today), e.g. W = 0 or S = S0.
    // states: numberOfStateVariables() values, filled by the model.
    virtual void initializeStates(std::span<double> states) const = 0;

    // One diffusion step, from date `timeIndex` to date `timeIndex + 1` of the master timeline.
    //   timeIndex     : index of the step (0 .. n_steps-1); selects the step length dt
    //   rndNumbers    : the correlated standard normals of this step and this path,
    //                   numberOfBrownians() values
    //   previousStates: states at date `timeIndex`, numberOfStateVariables() values
    //   states        : OUTPUT, states at date `timeIndex + 1`, numberOfStateVariables() values
    // The model only contains its one-step recursion; the engine owns the loop over steps
    // and over paths, which is what lets it advance several models together step by step.
    virtual void generateStates(size_t                    timeIndex,
                                std::span<const double>   rndNumbers,
                                std::span<const double>   previousStates,
                                std::span<double>         states) const = 0;

    // Observable function: queries model object passing stateVariables to reconstruct observable (e.g. Spot)
    virtual double observable(size_t date_idx, const Matrix& stateVariables) const = 0;

    // Convenience alias for Spot observable
    virtual double spot(size_t date_idx, const Matrix& stateVariables) const {
        return observable(date_idx, stateVariables);
    }
};

// ----------------------------------------------------------------------------
// 4. Curves & Correlations
// ----------------------------------------------------------------------------
// Discount curve: df(d) = P(0, d) is the price today of 1 unit of the curve's currency paid on
// date d. Deterministic: interest rates are not simulated in this example.
struct ICurve : IMarketObject {
    virtual double df(Date d) const = 0;
};

// Correlation between two risk factors (see Market::setCorrelation for the naming convention).
struct ICorrelation : IMarketObject {
    virtual double value() const = 0;
};

} // namespace mc
