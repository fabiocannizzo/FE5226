#pragma once

#include "Interfaces.h"
#include <algorithm>
#include <array>
#include <vector>

namespace mc {

// ----------------------------------------------------------------------------
// Building blocks shared by all pricers
// ----------------------------------------------------------------------------

// Reference to a model of the engine's master list, resolved at bind time
// Handle to one model of the engine's master list. Before binding it only knows the model TAG;
// bind() resolves it once into (column, pointer), so evaluating a path is a plain array access.
// An empty tag means "not used" (e.g. the numerator of USD.JPY, or the FX model of a USD flow).
struct ModelRef {
    std::string    tag;
    size_t         col = 0;
    const IModelMC* ptr = nullptr;

    bool active() const { return !tag.empty(); }

    // Position of the tag in the engine's sorted model list = index of the model's state matrix in
    // the `paths` vector passed to pricePath().
    void bind(const std::vector<std::string>& master_models, const std::vector<const IModelMC*>& models) {
        if (!active())
            return;
        auto it = std::find(master_models.begin(), master_models.end(), tag);
        MYASSERT(it != master_models.end(), "Model not found in master models: " << tag);
        col = static_cast<size_t>(std::distance(master_models.begin(), it));
        ptr = models[col];
    }
    // Value of the model's observable (the spot) at timeline column `date_col` on the current path.
    double observe(size_t date_col, const std::vector<Matrix>& paths) const {
        return ptr->observable(date_col, paths[col]);
    }
};

// Model tag = "<FAMILY>.<ASSET>" for the simulated USD asset `leg`, requested on behalf of
// `pair` (the family comes from the config, see PricerConfig::family).
inline std::string modelTag(const PricerConfig& config, const std::string& leg, const std::string& pair) {
    return config.family(leg, pair) + "." + leg;
}

// Spot of an FX pair BASE.QUOTE (price of one BASE in QUOTE) built from the USD-simulated
// assets: S(BASE.QUOTE) = X(BASE.USD) / X(QUOTE.USD), with X(USD.USD) = 1.
// Direct USD pairs use one model, crosses (e.g. EUR.GBP) and inverse pairs (USD.JPY) are derived.
// Spot of an FX pair BASE.QUOTE built from the simulated "<CCY>.USD" assets:
//     S(BASE.QUOTE) = X(BASE.USD) / X(QUOTE.USD),   with X(USD.USD) = 1.
//   - EUR.USD: one model (numerator only)
//   - EUR.GBP: two models, EUR.USD / GBP.USD (a cross), inherits their correlation
//   - USD.JPY: one model, 1 / JPY.USD (inverse)
// Every pricer reads its underlyings through this class, so direct, cross and inverse pairs are
// handled uniformly and the engine only ever simulates USD pairs.
class PairObservable {
public:
    // Split "BASE.QUOTE" and request a model for each non-USD side.
    PairObservable(const std::string& pair, const PricerConfig& config) {
        auto dot = pair.find('.');
        MYASSERT(dot != std::string::npos, "Invalid FX pair: " << pair);
        const std::string base = pair.substr(0, dot), quote = pair.substr(dot + 1);
        MYASSERT(base != quote, "Invalid FX pair: " << pair);
        if (base != NUMERAIRE_CCY)
            m_num.tag = modelTag(config, base + "." + NUMERAIRE_CCY, pair);
        if (quote != NUMERAIRE_CCY)
            m_den.tag = modelTag(config, quote + "." + NUMERAIRE_CCY, pair);
    }

    // Model tags this observable needs (0, 1 or 2).
    std::vector<std::string> models() const {
        std::vector<std::string> out;
        if (m_num.active()) out.push_back(m_num.tag);
        if (m_den.active()) out.push_back(m_den.tag);
        return out;
    }
    // Bind the (up to two) models of the pair.
    void bind(const std::vector<std::string>& master_models, const std::vector<const IModelMC*>& models) {
        m_num.bind(master_models, models);
        m_den.bind(master_models, models);
    }
    // Spot of the pair at timeline column `date_col` on the current path.
    double value(size_t date_col, const std::vector<Matrix>& paths) const {
        double v = m_num.active() ? m_num.observe(date_col, paths) : 1.0;
        if (m_den.active())
            v /= m_den.observe(date_col, paths);
        return v;
    }

private:
    ModelRef m_num, m_den;
};

// Base of all leaf pricers: one cash flow, determined at 'fixing' and paid at 'settle' in 'ccy'.
// Valuation on the fixing date Tf is the conditional expectation of the cash flow:
//    V(Tf) = payoff * P_ccy(Tf,Ts) * FX(ccy.USD)(Tf)            [FX = 1 for USD]
// and the numeraire (USD) discounts it to today with P_USD(0,Tf). The deterministic part
// P_USD(0,Tf) * P_ccy(Tf,Ts) is the 'scale' computed by the engine; the contingent FX spot
// is read from the simulated path here. With the FX drift implied by the curves this makes
// the foreign cash flow arbitrage-consistent with the foreign curve.
// Common base of all single-cash-flow pricers. The contract fixes its amount on `fixing` and
// pays it on `settle` >= `fixing` in currency `ccy`. On the fixing date the amount is KNOWN, so
// the value at that date is its conditional expectation; with deterministic rates:
//     V(Tf) = payoff * P_ccy(Tf, Ts)   [* X_ccy.USD(Tf) to express it in USD]
// Discounting to today with the USD numeraire gives, per path,
//     PV = P_USD(0,Tf) * P_ccy(Tf,Ts) * X_ccy.USD(Tf) * payoff
// The first two factors are path independent and arrive in `scale` (computed by the engine from
// the market); the contingent FX spot X(Tf) is read from the simulated path here (fx()).
// For a USD payment X = 1 and the product reduces to P_USD(0, Ts).
class LeafPricer : public IPricer {
public:
    LeafPricer(std::string ccy, Date fixing, Date settle, const PricerConfig& config)
        : m_ccy(std::move(ccy)), m_fixing(fixing), m_settle(settle) {
        MYASSERT(settle >= fixing, "Settlement date " << settle << " must not precede fixing date " << fixing);
        // A foreign-currency payment needs the FX rate "<ccy>.USD" at the fixing date: ask the engine for
        // that model, even if the underlying of the trade is a different pair.
        if (m_ccy != NUMERAIRE_CCY)
            m_fx.tag = modelTag(config, m_ccy + "." + NUMERAIRE_CCY, m_ccy + "." + NUMERAIRE_CCY);
    }

    Date        fixingDate() const override { return m_fixing; }
    Date        settlementDate() const override { return m_settle; }
    std::string paymentCurrency() const override { return m_ccy; }

protected:
    // Resolve the FX model and the timeline column at which the FX rate is read.
    void bindFx(const std::vector<Date>& master_dates, const std::vector<std::string>& master_models,
                const std::vector<const IModelMC*>& models) {
        m_fx.bind(master_models, models);
        if (m_fx.active()) {
            // A fixing in the past is converted with today's spot (first column of the timeline)
            // A fixing already in the past is converted with today's spot (the first column of the timeline).
            const Date d  = std::max(m_fixing, master_dates.front());
            auto       it = std::lower_bound(master_dates.begin(), master_dates.end(), d);
            MYASSERT(it != master_dates.end() && *it == d, "Fixing date not found in master timeline: " << d);
            m_fx_date_col = static_cast<size_t>(std::distance(master_dates.begin(), it));
        }
    }
    // Contingent FX spot ccy.USD at the fixing date on the current path; 1 for USD payments.
    double fx(const std::vector<Matrix>& paths) const {
        return m_fx.active() ? m_fx.observe(m_fx_date_col, paths) : 1.0;
    }
    void addFxModel(std::vector<std::string>& out) const {
        if (m_fx.active()) out.push_back(m_fx.tag);
    }
    // Index of date `d` in the sorted master timeline (binary search); `d` must be present.
    static size_t columnOf(const std::vector<Date>& master_dates, Date d) {
        auto it = std::lower_bound(master_dates.begin(), master_dates.end(), d);
        MYASSERT(it != master_dates.end() && *it == d, "Date not found in master timeline: " << d);
        return static_cast<size_t>(std::distance(master_dates.begin(), it));
    }

private:
    std::string m_ccy;
    Date        m_fixing;
    Date        m_settle;
    ModelRef    m_fx;
    size_t      m_fx_date_col = 0;
};

// ----------------------------------------------------------------------------
// 1. Vanilla (European) Option Pricer
// ----------------------------------------------------------------------------
// Pricer of PayoffVanilla. Life cycle in the constructor: if the expiry fixing is already in the
// events the spot is frozen and no date is requested; otherwise the expiry is requested.
class VanillaMCPricer : public LeafPricer {
public:
    VanillaMCPricer(const std::string& pair, double notional, double strike, int call_put, Date expiry,
                    Date settle, const std::string& ccy, const PricerConfig& config, const Events& events)
        : LeafPricer(ccy, expiry, settle, config)
        , m_obs(pair, config)
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_expiry(expiry) {
        // Lifecycle against historical events
        // Lifecycle: has the observation already happened?
        if (lookupEvent(events, pair, m_expiry, m_past_val))
            m_is_historical = true;
        else
            m_future_dates.push_back(m_expiry);
    }

    // Only FUTURE dates: past observations never enter the simulation.
    std::vector<Date> requiredDates() const override { return m_future_dates; }
    // Models of the underlying, plus the FX model when the payment currency is not USD.
    std::vector<std::string> modelNeeded() const override {
        auto out = m_obs.models();
        addFxModel(out);
        return out;
    }

    // Called once after the engine has fixed the timeline and the models: resolve dates and tags
    // into integer positions, so that pricePath() does no searching.
    void bindIndices(const std::vector<Date>& master_dates, const std::vector<std::string>& master_models,
                     const std::vector<const IModelMC*>& models) override {
        if (!m_is_historical)
            m_timeline_col = columnOf(master_dates, m_expiry);
        m_obs.bind(master_models, models);
        bindFx(master_dates, master_models, models);
    }

    // Value of the contract on one path, in USD as of today:
    //     scale * X_ccy.USD(Tf) * notional * max(omega * (S - K), 0)
    // S is the recorded fixing if the expiry is in the past, the simulated spot otherwise.
    double pricePath(const std::vector<Matrix>& model_paths, double scale) const override {
        double S = m_is_historical ? m_past_val : m_obs.value(m_timeline_col, model_paths);
        double payoff = m_notional * std::max(m_omega * (S - m_strike), 0.0);
        return scale * fx(model_paths) * payoff;
    }

private:
    PairObservable    m_obs;
    double            m_notional;
    double            m_strike;
    int               m_omega;
    Date              m_expiry;
    bool              m_is_historical = false;
    double            m_past_val      = 0.0;
    std::vector<Date> m_future_dates;
    size_t            m_timeline_col = 0;
};

// ----------------------------------------------------------------------------
// 2. Arithmetic Average Asian Option Pricer
// ----------------------------------------------------------------------------
// Pricer of PayoffAsian. Some fixings may be past (recorded in the events) and others future
// (simulated). m_future_dates holds the sorted unique future dates; m_future_slot[j] maps fixing j
// of the schedule to its position in that list; m_bound_timeline_cols maps the list to timeline
// columns. The cash flow's fixing date is the LAST fixing of the schedule.
class AsianMCPricer : public LeafPricer {
public:
    AsianMCPricer(const std::string& pair, double notional, double strike, int call_put,
                  const std::vector<Date>& fix_dates, Date settle, const std::string& ccy,
                  const PricerConfig& config, const Events& events)
        : LeafPricer(ccy, fix_dates.empty() ? 0 : *std::max_element(fix_dates.begin(), fix_dates.end()),
                     settle, config)
        , m_obs(pair, config)
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_fix_dates(fix_dates) {
        MYASSERT(!m_fix_dates.empty(), "Asian option requires at least one fixing date");

        const size_t m = m_fix_dates.size();
        m_is_historical.assign(m, false);
        m_past_val.assign(m, 0.0);
        m_future_slot.assign(m, 0);

        for (size_t j = 0; j < m; ++j) {
            double val = 0.0;
            if (lookupEvent(events, pair, m_fix_dates[j], val)) {
                m_is_historical[j] = true;
                m_past_val[j]      = val;
            } else {
                m_future_dates.push_back(m_fix_dates[j]);
            }
        }
        // Sorted unique future dates first, then map each future fixing to its slot
        // First sort and de-duplicate the future dates, THEN compute each fixing's slot: computing the
        // slot before sorting would give wrong columns for an unsorted schedule.
        std::sort(m_future_dates.begin(), m_future_dates.end());
        m_future_dates.erase(std::unique(m_future_dates.begin(), m_future_dates.end()), m_future_dates.end());
        for (size_t j = 0; j < m; ++j) {
            if (!m_is_historical[j])
                m_future_slot[j] = static_cast<size_t>(
                    std::lower_bound(m_future_dates.begin(), m_future_dates.end(), m_fix_dates[j]) -
                    m_future_dates.begin());
        }
    }

    // Only the future fixings.
    std::vector<Date> requiredDates() const override { return m_future_dates; }
    std::vector<std::string> modelNeeded() const override {
        auto out = m_obs.models();
        addFxModel(out);
        return out;
    }

    void bindIndices(const std::vector<Date>& master_dates, const std::vector<std::string>& master_models,
                     const std::vector<const IModelMC*>& models) override {
        m_bound_timeline_cols.clear();
        for (Date d : m_future_dates)
            m_bound_timeline_cols.push_back(columnOf(master_dates, d));
        m_obs.bind(master_models, models);
        bindFx(master_dates, master_models, models);
    }

    // Average of the fixings (recorded or simulated) -> call / put on the average -> scale and FX.
    double pricePath(const std::vector<Matrix>& model_paths, double scale) const override {
        double       sum = 0.0;
        const size_t m   = m_fix_dates.size();
        for (size_t j = 0; j < m; ++j) {
            sum += m_is_historical[j] ? m_past_val[j]
                                      : m_obs.value(m_bound_timeline_cols[m_future_slot[j]], model_paths);
        }
        double avg    = sum / static_cast<double>(m);
        double payoff = m_notional * std::max(m_omega * (avg - m_strike), 0.0);
        return scale * fx(model_paths) * payoff;
    }

private:
    PairObservable      m_obs;
    double              m_notional;
    double              m_strike;
    int                 m_omega;
    std::vector<Date>   m_fix_dates;
    std::vector<bool>   m_is_historical;
    std::vector<double> m_past_val;
    std::vector<size_t> m_future_slot;
    std::vector<Date>   m_future_dates;
    std::vector<size_t> m_bound_timeline_cols;
};

// ----------------------------------------------------------------------------
// 3. Multi-Asset Worst-Of Basket Call Pricer
// ----------------------------------------------------------------------------
// Pricer of PayoffWorstOfCall: two underlyings, each observed on the expiry date, either from the
// events (past) or from the simulation. If both are past no date is requested at all.
class WorstOfMCPricer : public LeafPricer {
public:
    WorstOfMCPricer(const std::string& pair0, const std::string& pair1, double spot0_0, double spot0_1,
                    double notional, double strike, Date expiry, Date settle, const std::string& ccy,
                    const PricerConfig& config, const Events& events)
        : LeafPricer(ccy, expiry, settle, config)
        , m_obs{ PairObservable(pair0, config), PairObservable(pair1, config) }
        , m_spot0{ spot0_0, spot0_1 }
        , m_notional(notional)
        , m_strike(strike)
        , m_expiry(expiry) {
        const std::array<const std::string*, 2> pairs{ &pair0, &pair1 };
        for (size_t a = 0; a < 2; ++a) {
            double val = 0.0;
            if (lookupEvent(events, *pairs[a], m_expiry, val)) {
                m_is_historical[a] = true;
                m_past_val[a]      = val;
            }
        }
        // Ask for the expiry date only if at least one leg still has to be simulated.
        if (!m_is_historical[0] || !m_is_historical[1])
            m_future_dates.push_back(m_expiry);
    }

    std::vector<Date> requiredDates() const override { return m_future_dates; }
    std::vector<std::string> modelNeeded() const override {
        std::vector<std::string> out;
        for (const auto& o : m_obs) {
            auto m = o.models();
            out.insert(out.end(), m.begin(), m.end());
        }
        addFxModel(out);
        return out;
    }

    void bindIndices(const std::vector<Date>& master_dates, const std::vector<std::string>& master_models,
                     const std::vector<const IModelMC*>& models) override {
        if (!m_future_dates.empty())
            m_timeline_col = columnOf(master_dates, m_expiry);
        for (auto& o : m_obs)
            o.bind(master_models, models);
        bindFx(master_dates, master_models, models);
    }

    // Worst performance of the two legs -> call on it -> scale and FX.
    double pricePath(const std::vector<Matrix>& model_paths, double scale) const override {
        double worst = 0.0;
        for (size_t a = 0; a < 2; ++a) {
            double S    = m_is_historical[a] ? m_past_val[a] : m_obs[a].value(m_timeline_col, model_paths);
            double perf = S / m_spot0[a];
            worst       = (a == 0) ? perf : std::min(worst, perf);
        }
        double payoff = m_notional * std::max(worst - m_strike, 0.0);
        return scale * fx(model_paths) * payoff;
    }

private:
    std::array<PairObservable, 2> m_obs;
    std::array<double, 2>         m_spot0;
    double                        m_notional;
    double                        m_strike;
    Date                          m_expiry;
    std::array<bool, 2>           m_is_historical{ false, false };
    std::array<double, 2>         m_past_val{ 0.0, 0.0 };
    std::vector<Date>             m_future_dates;
    size_t                        m_timeline_col = 0;
};

} // namespace mc
