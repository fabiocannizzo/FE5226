#pragma once

#include "Interfaces.h"
#include <algorithm>
#include <array>
#include <vector>

namespace mc {

// ----------------------------------------------------------------------------
// 1. Vanilla (European) Option Pricer
// ----------------------------------------------------------------------------
class VanillaMCPricer : public IPricer {
public:
    VanillaMCPricer(std::string                pair,
                    double                     notional,
                    double                     strike,
                    int                        call_put,
                    Date                       expiry,
                    Date                       settle,
                    std::string                ccy,
                    std::string                model_name,
                    const std::vector<Events>& events)
        : m_pair(std::move(pair))
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_expiry(expiry)
        , m_settle(settle)
        , m_ccy(std::move(ccy))
        , m_model_name(std::move(model_name)) {
        // Lifecycle against historical events
        double past_val = 0.0;
        if (lookupEvent(events, m_pair, m_expiry, past_val)) {
            m_is_historical = true;
            m_past_val      = past_val;
        } else {
            m_is_historical = false;
            m_future_dates.push_back(m_expiry);
        }
    }

    std::vector<Date> requiredDates() const override { return m_future_dates; }
    std::vector<std::string> modelNeeded() const override { return { m_model_name }; }
    Date settlementDate() const override { return m_settle; }
    std::string baseCurrency() const override { return m_ccy; }

    void bindIndices(const std::vector<Date>&          master_dates,
                     const std::vector<std::string>&   master_models,
                     const std::vector<const IModel*>& models) override {
        if (!m_is_historical) {
            auto it = std::lower_bound(master_dates.begin(), master_dates.end(), m_expiry);
            MYASSERT(it != master_dates.end() && *it == m_expiry,
                     "Expiry date not found in master timeline: " << m_expiry);
            m_timeline_col = static_cast<size_t>(std::distance(master_dates.begin(), it));
        }

        auto mit = std::find(master_models.begin(), master_models.end(), m_model_name);
        MYASSERT(mit != master_models.end(), "Model not found in master models: " << m_model_name);
        m_model_col = static_cast<size_t>(std::distance(master_models.begin(), mit));
        m_model     = models[m_model_col];
    }

    double pricePath(const std::vector<Matrix>& model_paths,
                     double                     df_settle,
                     double                     fx_to_base) const override {
        double S = m_is_historical
                       ? m_past_val
                       : m_model->observable(m_timeline_col, model_paths[m_model_col]);
        double payoff = m_notional * std::max(m_omega * (S - m_strike), 0.0);
        return df_settle * fx_to_base * payoff;
    }

private:
    std::string       m_pair;
    double            m_notional;
    double            m_strike;
    int               m_omega;
    Date              m_expiry;
    Date              m_settle;
    std::string       m_ccy;
    std::string       m_model_name;
    bool              m_is_historical = false;
    double            m_past_val      = 0.0;
    std::vector<Date> m_future_dates;
    size_t            m_timeline_col = 0;
    size_t            m_model_col    = 0;
    const IModel*     m_model        = nullptr;
};

// ----------------------------------------------------------------------------
// 2. Arithmetic Average Asian Option Pricer
// ----------------------------------------------------------------------------
class AsianMCPricer : public IPricer {
public:
    AsianMCPricer(std::string                pair,
                  double                     notional,
                  double                     strike,
                  int                        call_put,
                  std::vector<Date>          fix_dates,
                  Date                       settle,
                  std::string                ccy,
                  std::string                model_name,
                  const std::vector<Events>& events)
        : m_pair(std::move(pair))
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_fix_dates(std::move(fix_dates))
        , m_settle(settle)
        , m_ccy(std::move(ccy))
        , m_model_name(std::move(model_name)) {
        MYASSERT(!m_fix_dates.empty(), "Asian option requires at least one fixing date");

        const size_t m = m_fix_dates.size();
        m_is_historical.assign(m, false);
        m_past_val.assign(m, 0.0);
        m_future_slot.assign(m, 0);

        for (size_t j = 0; j < m; ++j) {
            Date   d   = m_fix_dates[j];
            double val = 0.0;
            if (lookupEvent(events, m_pair, d, val)) {
                m_is_historical[j] = true;
                m_past_val[j]      = val;
            } else {
                m_is_historical[j] = false;
                auto it            = std::find(m_future_dates.begin(), m_future_dates.end(), d);
                if (it == m_future_dates.end()) {
                    m_future_slot[j] = m_future_dates.size();
                    m_future_dates.push_back(d);
                } else {
                    m_future_slot[j] = std::distance(m_future_dates.begin(), it);
                }
            }
        }
        std::sort(m_future_dates.begin(), m_future_dates.end());
    }

    std::vector<Date> requiredDates() const override { return m_future_dates; }
    std::vector<std::string> modelNeeded() const override { return { m_model_name }; }
    Date settlementDate() const override { return m_settle; }
    std::string baseCurrency() const override { return m_ccy; }

    void bindIndices(const std::vector<Date>&          master_dates,
                     const std::vector<std::string>&   master_models,
                     const std::vector<const IModel*>& models) override {
        m_bound_timeline_cols.clear();
        for (Date d : m_future_dates) {
            auto it = std::lower_bound(master_dates.begin(), master_dates.end(), d);
            MYASSERT(it != master_dates.end() && *it == d,
                     "Asian fixing date not found in master timeline: " << d);
            m_bound_timeline_cols.push_back(
                static_cast<size_t>(std::distance(master_dates.begin(), it)));
        }

        auto mit = std::find(master_models.begin(), master_models.end(), m_model_name);
        MYASSERT(mit != master_models.end(), "Model not found in master models: " << m_model_name);
        m_model_col = static_cast<size_t>(std::distance(master_models.begin(), mit));
        m_model     = models[m_model_col];
    }

    double pricePath(const std::vector<Matrix>& model_paths,
                     double                     df_settle,
                     double                     fx_to_base) const override {
        double        sum  = 0.0;
        const size_t  m    = m_fix_dates.size();
        const Matrix& path = model_paths[m_model_col];

        for (size_t j = 0; j < m; ++j) {
            if (m_is_historical[j]) {
                sum += m_past_val[j];
            } else {
                size_t slot = m_future_slot[j];
                size_t col  = m_bound_timeline_cols[slot];
                sum += m_model->observable(col, path);
            }
        }

        double avg    = sum / static_cast<double>(m);
        double payoff = m_notional * std::max(m_omega * (avg - m_strike), 0.0);
        return df_settle * fx_to_base * payoff;
    }

private:
    std::string         m_pair;
    double              m_notional;
    double              m_strike;
    int                 m_omega;
    std::vector<Date>   m_fix_dates;
    Date                m_settle;
    std::string         m_ccy;
    std::string         m_model_name;
    std::vector<bool>   m_is_historical;
    std::vector<double> m_past_val;
    std::vector<size_t> m_future_slot;
    std::vector<Date>   m_future_dates;
    std::vector<size_t> m_bound_timeline_cols;
    size_t              m_model_col = 0;
    const IModel*       m_model     = nullptr;
};

// ----------------------------------------------------------------------------
// 3. Multi-Asset Worst-Of Basket Call Pricer
// ----------------------------------------------------------------------------
class WorstOfMCPricer : public IPricer {
public:
    WorstOfMCPricer(std::string                pair0,
                    std::string                pair1,
                    double                     spot0_0,
                    double                     spot0_1,
                    double                     notional,
                    double                     strike,
                    Date                       expiry,
                    Date                       settle,
                    std::string                ccy,
                    std::string                model0,
                    std::string                model1,
                    const std::vector<Events>& events)
        : m_pairs{ std::move(pair0), std::move(pair1) }
        , m_spot0{ spot0_0, spot0_1 }
        , m_notional(notional)
        , m_strike(strike)
        , m_expiry(expiry)
        , m_settle(settle)
        , m_ccy(std::move(ccy))
        , m_model_names{ std::move(model0), std::move(model1) } {
        // Lifecycle each asset
        for (size_t a = 0; a < 2; ++a) {
            double val = 0.0;
            if (lookupEvent(events, m_pairs[a], m_expiry, val)) {
                m_is_historical[a] = true;
                m_past_val[a]      = val;
            } else {
                m_is_historical[a] = false;
            }
        }
        if (!m_is_historical[0] || !m_is_historical[1]) {
            m_future_dates.push_back(m_expiry);
        }
    }

    std::vector<Date> requiredDates() const override { return m_future_dates; }
    std::vector<std::string> modelNeeded() const override {
        return { m_model_names[0], m_model_names[1] };
    }
    Date settlementDate() const override { return m_settle; }
    std::string baseCurrency() const override { return m_ccy; }

    void bindIndices(const std::vector<Date>&          master_dates,
                     const std::vector<std::string>&   master_models,
                     const std::vector<const IModel*>& models) override {
        if (!m_future_dates.empty()) {
            auto it = std::lower_bound(master_dates.begin(), master_dates.end(), m_expiry);
            MYASSERT(it != master_dates.end() && *it == m_expiry, "Expiry date not in master timeline");
            m_timeline_col = static_cast<size_t>(std::distance(master_dates.begin(), it));
        }
        for (size_t a = 0; a < 2; ++a) {
            auto mit = std::find(master_models.begin(), master_models.end(), m_model_names[a]);
            MYASSERT(mit != master_models.end(), "Model not found: " << m_model_names[a]);
            m_model_cols[a] = static_cast<size_t>(std::distance(master_models.begin(), mit));
            m_models[a]     = models[m_model_cols[a]];
        }
    }

    double pricePath(const std::vector<Matrix>& model_paths,
                     double                     df_settle,
                     double                     fx_to_base) const override {
        double S0 = m_is_historical[0]
                        ? m_past_val[0]
                        : m_models[0]->observable(m_timeline_col, model_paths[m_model_cols[0]]);
        double S1 = m_is_historical[1]
                        ? m_past_val[1]
                        : m_models[1]->observable(m_timeline_col, model_paths[m_model_cols[1]]);

        double perf0  = S0 / m_spot0[0];
        double perf1  = S1 / m_spot0[1];
        double worst  = std::min(perf0, perf1);
        double payoff = m_notional * std::max(worst - m_strike, 0.0);
        return df_settle * fx_to_base * payoff;
    }

private:
    std::array<std::string, 2>   m_pairs;
    std::array<double, 2>        m_spot0;
    double                       m_notional;
    double                       m_strike;
    Date                         m_expiry;
    Date                         m_settle;
    std::string                  m_ccy;
    std::array<std::string, 2>   m_model_names;
    std::array<bool, 2>          m_is_historical{ false, false };
    std::array<double, 2>        m_past_val{ 0.0, 0.0 };
    std::vector<Date>            m_future_dates;
    size_t                       m_timeline_col = 0;
    std::array<size_t, 2>        m_model_cols{ 0, 0 };
    std::array<const IModel*, 2> m_models{ nullptr, nullptr };
};

} // namespace mc
