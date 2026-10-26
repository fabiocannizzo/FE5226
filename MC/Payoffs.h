#pragma once

#include "Interfaces.h"
#include "Pricers.h"
#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace mc {

// ----------------------------------------------------------------------------
// 1. Vanilla (European) Option Payoff: serves as contract spec & pricer factory
// ----------------------------------------------------------------------------
class PayoffVanilla : public IPayoff {
public:
    PayoffVanilla(std::string pair,
                  double      notional,
                  double      strike,
                  int         call_put,
                  Date        expiry,
                  Date        settle,
                  std::string ccy          = "USD",
                  std::string model_family = "BS_CLOSED")
        : m_pair(std::move(pair))
        , m_pay_ccy(std::move(ccy))
        , m_model_family(std::move(model_family))
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_expiry(expiry)
        , m_settle(settle) {}

    std::vector<std::string> underlyingNames() const override { return { m_pair }; }
    Date                     settlementDate() const override { return m_settle; }
    std::string              payoffCurrency() const override { return m_pay_ccy; }

    std::shared_ptr<IPricer> getPricer(const std::vector<Events>& events) const override {
        return std::make_shared<VanillaMCPricer>(
            m_pair, m_notional, m_strike, m_omega, m_expiry, m_settle,
            m_pay_ccy, m_model_family + "." + m_pair, events);
    }

private:
    std::string m_pair;
    std::string m_pay_ccy;
    std::string m_model_family;
    double      m_notional;
    double      m_strike;
    int         m_omega; // +1 for Call, -1 for Put
    Date        m_expiry;
    Date        m_settle;
};

using PayoffEuropean = PayoffVanilla;

// ----------------------------------------------------------------------------
// 2. Arithmetic Average Asian Option Payoff
// ----------------------------------------------------------------------------
class PayoffAsian : public IPayoff {
public:
    PayoffAsian(std::string       pair,
                double            notional,
                double            strike,
                int               call_put,
                std::vector<Date> fix_dates,
                Date              settle,
                std::string       ccy          = "USD",
                std::string       model_family = "BS_CLOSED")
        : m_pair(std::move(pair))
        , m_pay_ccy(std::move(ccy))
        , m_model_family(std::move(model_family))
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_fix_dates(std::move(fix_dates))
        , m_settle(settle) {
        MYASSERT(!m_fix_dates.empty(), "Asian option requires at least one fixing date");
    }

    std::vector<std::string> underlyingNames() const override { return { m_pair }; }
    Date                     settlementDate() const override { return m_settle; }
    std::string              payoffCurrency() const override { return m_pay_ccy; }

    std::shared_ptr<IPricer> getPricer(const std::vector<Events>& events) const override {
        return std::make_shared<AsianMCPricer>(
            m_pair, m_notional, m_strike, m_omega, m_fix_dates, m_settle,
            m_pay_ccy, m_model_family + "." + m_pair, events);
    }

private:
    std::string       m_pair;
    std::string       m_pay_ccy;
    std::string       m_model_family;
    double            m_notional;
    double            m_strike;
    int               m_omega;
    std::vector<Date> m_fix_dates;
    Date              m_settle;
};

// ----------------------------------------------------------------------------
// 3. Multi-Asset Worst-Of Basket Call Option Payoff
// ----------------------------------------------------------------------------
class PayoffWorstOfCall : public IPayoff {
public:
    PayoffWorstOfCall(std::string pair0,
                      std::string pair1,
                      double      spot0_0,
                      double      spot0_1,
                      double      notional,
                      double      strike,
                      Date        expiry,
                      Date        settle,
                      std::string ccy    = "USD",
                      std::string model0 = "BS_CLOSED",
                      std::string model1 = "BS_CLOSED")
        : m_pairs{ std::move(pair0), std::move(pair1) }
        , m_spot0{ spot0_0, spot0_1 }
        , m_pay_ccy(std::move(ccy))
        , m_models{ std::move(model0), std::move(model1) }
        , m_notional(notional)
        , m_strike(strike)
        , m_expiry(expiry)
        , m_settle(settle) {}

    std::vector<std::string> underlyingNames() const override {
        return { m_pairs[0], m_pairs[1] };
    }
    Date        settlementDate() const override { return m_settle; }
    std::string payoffCurrency() const override { return m_pay_ccy; }

    std::shared_ptr<IPricer> getPricer(const std::vector<Events>& events) const override {
        return std::make_shared<WorstOfMCPricer>(
            m_pairs[0], m_pairs[1], m_spot0[0], m_spot0[1], m_notional, m_strike,
            m_expiry, m_settle, m_pay_ccy,
            m_models[0] + "." + m_pairs[0],
            m_models[1] + "." + m_pairs[1],
            events);
    }

private:
    std::array<std::string, 2> m_pairs;
    std::array<double, 2>      m_spot0;
    std::string                m_pay_ccy;
    std::array<std::string, 2> m_models;
    double                     m_notional;
    double                     m_strike;
    Date                       m_expiry;
    Date                       m_settle;
};

} // namespace mc
