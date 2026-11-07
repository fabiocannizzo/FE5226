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
// European option on an FX pair BASE.QUOTE:  notional * max(omega * (S(expiry) - K), 0)
// omega = +1 (call) / -1 (put). The amount is determined at `expiry` (fixing date) and paid at
// `settle` (>= expiry) in currency `ccy`. The payment currency may differ from the QUOTE currency
// of the pair: the pricer converts it.
class PayoffVanilla : public IPayoff {
public:
    PayoffVanilla(std::string pair,
                  double      notional,
                  double      strike,
                  int         call_put,
                  Date        expiry,
                  Date        settle,
                  std::string ccy = "USD")
        : m_pair(std::move(pair))
        , m_pay_ccy(std::move(ccy))
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_expiry(expiry)
        , m_settle(settle) {
        // A contract cannot pay before its amount is known.
        MYASSERT(settle >= expiry, "Settlement date must not precede expiry");
    }

    // Factory: hands the contract terms, the config and the known fixings to the pricer, which then
    // decides what remains to be simulated. The payoff itself keeps no state.
    std::shared_ptr<IPricer> getPricer(const Events&       events,
                                       const PricerConfig& config = PricerConfig{}) const override {
        return std::make_shared<VanillaMCPricer>(
            m_pair, m_notional, m_strike, m_omega, m_expiry, m_settle, m_pay_ccy, config, events);
    }

    std::vector<std::string> getCashFlowCurrencies(const Events& /*events*/) const override {
        return { m_pay_ccy };
    }

    // Contractual dates at which market data is read (independent of events: a past date is still an
    // observation date of the contract).
    std::vector<Date> observationDates(const Events& /*events*/) const override {
        return { m_expiry };
    }

    Date maturityDate(const Events& /*events*/) const override {
        return m_expiry;
    }

    std::vector<Date> settlementDates(const Events& /*events*/) const override {
        return { m_settle };
    }

private:
    std::string m_pair;
    std::string m_pay_ccy;
    double      m_notional;
    double      m_strike;
    // Contract terms; all immutable after construction.
    int         m_omega; // +1 for Call, -1 for Put
    Date        m_expiry;
    Date        m_settle;
};

// Alternative name for the same class.
using PayoffEuropean = PayoffVanilla;

// ----------------------------------------------------------------------------
// 2. Arithmetic Average Asian Option Payoff
// ----------------------------------------------------------------------------
// Arithmetic average (Asian) option:  notional * max(omega * (mean_j S(t_j) - K), 0)
// The schedule may contain PAST dates: the average then mixes recorded fixings (from Events)
// with simulated ones. Payment on `settle`, which must not precede the last fixing.
class PayoffAsian : public IPayoff {
public:
    PayoffAsian(std::string       pair,
                double            notional,
                double            strike,
                int               call_put,
                std::vector<Date> fix_dates,
                Date              settle,
                std::string       ccy = "USD")
        : m_pair(std::move(pair))
        , m_pay_ccy(std::move(ccy))
        , m_notional(notional)
        , m_strike(strike)
        , m_omega(call_put)
        , m_fix_dates(std::move(fix_dates))
        , m_settle(settle) {
        MYASSERT(!m_fix_dates.empty(), "Asian option requires at least one fixing date");
        MYASSERT(settle >= *std::max_element(m_fix_dates.begin(), m_fix_dates.end()),
                 "Settlement date must not precede the last fixing");
    }

    std::shared_ptr<IPricer> getPricer(const Events&       events,
                                       const PricerConfig& config = PricerConfig{}) const override {
        return std::make_shared<AsianMCPricer>(
            m_pair, m_notional, m_strike, m_omega, m_fix_dates, m_settle, m_pay_ccy, config, events);
    }

    std::vector<std::string> getCashFlowCurrencies(const Events& /*events*/) const override {
        return { m_pay_ccy };
    }

    std::vector<Date> observationDates(const Events& /*events*/) const override {
        return m_fix_dates;
    }

    Date maturityDate(const Events& /*events*/) const override {
        return m_fix_dates.empty() ? 0 : m_fix_dates.back();
    }

    std::vector<Date> settlementDates(const Events& /*events*/) const override {
        return { m_settle };
    }

private:
    std::string       m_pair;
    std::string       m_pay_ccy;
    double            m_notional;
    double            m_strike;
    int               m_omega;
    std::vector<Date> m_fix_dates;
    Date              m_settle;
};

// ----------------------------------------------------------------------------
// 3. Multi-Asset Worst-Of Basket Call Option Payoff
// ----------------------------------------------------------------------------
// Worst-of call on two FX pairs, in performance terms:
//     notional * max( min(S1(T)/S1_0, S2(T)/S2_0) - K, 0 )
// `spot0_*` are the CONTRACTUAL initial levels (fixed at trade inception), not today's market
// spot; the strike K is a performance level (1.00 = at the money). Needs two correlated models.
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
                      std::string ccy = "USD")
        : m_pairs{ std::move(pair0), std::move(pair1) }
        , m_spot0{ spot0_0, spot0_1 }
        , m_pay_ccy(std::move(ccy))
        , m_notional(notional)
        , m_strike(strike)
        , m_expiry(expiry)
        , m_settle(settle) {
        MYASSERT(settle >= expiry, "Settlement date must not precede expiry");
    }

    std::shared_ptr<IPricer> getPricer(const Events&       events,
                                       const PricerConfig& config = PricerConfig{}) const override {
        return std::make_shared<WorstOfMCPricer>(
            m_pairs[0], m_pairs[1], m_spot0[0], m_spot0[1], m_notional, m_strike,
            m_expiry, m_settle, m_pay_ccy, config, events);
    }

    std::vector<std::string> getCashFlowCurrencies(const Events& /*events*/) const override {
        return { m_pay_ccy };
    }

    std::vector<Date> observationDates(const Events& /*events*/) const override {
        return { m_expiry };
    }

    Date maturityDate(const Events& /*events*/) const override {
        return m_expiry;
    }

    std::vector<Date> settlementDates(const Events& /*events*/) const override {
        return { m_settle };
    }

private:
    std::array<std::string, 2> m_pairs;
    std::array<double, 2>      m_spot0;
    std::string                m_pay_ccy;
    double                     m_notional;
    double                     m_strike;
    Date                       m_expiry;
    Date                       m_settle;
};

} // namespace mc
