#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

using namespace std;

// cumulative standard normal distribution function
double cnorm(double x)
{
    return 0.5 * erfc(x / (-sqrt(2.0)));
}

// this is an abstract class, because it contains a pure virtual method
class EuropeanOption
{
    double m_notional;
    double m_tExpiry;
    double m_tSettle;
    double m_strike;
    bool m_isCall;

public:
    EuropeanOption(double notional, double tExpiry, double tSettle, double strike, bool isCall)
        : m_notional(notional), m_tExpiry(tExpiry), m_tSettle(tSettle), m_strike(strike), m_isCall(isCall) {}

    double strike() const { return m_strike; }
    bool isCall() const { return m_isCall; }

    virtual double undiscounted(double fwd, double stdev) = 0; // pure virtual

    double price(double today, double fwd, double rate, double vol)
    {
        double timeToExp = m_tExpiry - today;
        double timeToSettle = m_tSettle - today;
        double discfact = exp(-rate * timeToSettle);
        double stDev = sqrt(timeToExp) * vol;
        double undisc = undiscounted(fwd, stDev);
        return m_notional * discfact * undisc;
    }
};

class EuropeanVanilla : public EuropeanOption
{
public:
    using EuropeanOption::EuropeanOption; // reuse constructor of base class

    // implemente virtual method invoiled from the base class
    virtual double undiscounted(double fwd, double stdev)
    {
        double d1 = log(fwd / strike()) / stdev + 0.5 * stdev;
        double d2 = d1 - stdev;
        double nd1 = cnorm(d1);
        double nd2 = cnorm(d2);
        double u = fwd * nd1 - strike() * nd2;
        // if it is a put we use put/call parity:  F-K = Call-Put => Put = Call - (F-K)
        return isCall() ? u : u - (fwd - strike());
    }
};

class EuropeanDigital : public EuropeanOption
{
public:
    using EuropeanOption::EuropeanOption; // resuse constructor of new class

    // implemente virtual method invoiled from the base class
    virtual double undiscounted(double fwd, double stdev)
    {
        double d1 = log(fwd / strike()) / stdev + 0.5 * stdev;
        double d2 = d1 - stdev;
        double nd2 = cnorm(d2);

        // if it is a put we use digital put-call parity: P(S<K)+P(S>K)=1
        return isCall() ? nd2 : 1 - nd2;
    }
};


int main()
{
    // alternative syntax for a defining a new type (typedef)
    using europtr = shared_ptr<EuropeanOption>;

    vector<europtr> trades;
    // add various types of trades to the vector
    trades.push_back(europtr(new EuropeanVanilla(100, 2, 2.001, 10.1, true)));
    trades.push_back(europtr(new EuropeanVanilla(100, 2, 2.001, 9.9, false)));
    trades.push_back(europtr(new EuropeanDigital(100, 2, 2.001, 9.9, true)));
    trades.push_back(europtr(new EuropeanDigital(100, 2, 2.001, 9.9, false)));

    // compute the price for all options
    for (auto t : trades)
        cout << t->price(0, 10, 0.05, 0.30) << endl;

    return 0;
}

