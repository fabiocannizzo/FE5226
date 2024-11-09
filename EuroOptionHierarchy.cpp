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
    double m_tExpiry;   // absolute time specified in years
    double m_tSettle;   // absolute time specified in years
    double m_strike;
    bool   m_isCall;  // call or put option

public:
    EuropeanOption(double notional, double tExpiry, double tSettle, double strike, bool isCall)
        : m_notional(notional), m_tExpiry(tExpiry), m_tSettle(tSettle), m_strike(strike), m_isCall(isCall) {}

    double strike() const { return m_strike; }
    bool isCall() const { return m_isCall; }

    virtual double undiscounted(double fwd, double stdev) const = 0; // pure virtual

    double price(double today, double fwd, double rate, double vol) const
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
    /*
       We want to have a constructor with the same signature of the base class one and simlpy
       forward the arguments to the base class constructor. I.e. we would want to define a constructor
       like this:

       EuropeanVanilla(double notional, double tExpiry, double tSettle, double strike, bool isCall)
            : EuropeanOption(notional, tExpiry, tSettle, strike, isCall) {}

       This is tedious. What we can do instead is to simply inherits the constructors of the base class
    */
    using EuropeanOption::EuropeanOption; // inherit constructor of base class

    // implements virtual method declared in the base class
    double undiscounted(double fwd, double stdev) const override
    {
        double d1 = log(fwd / strike()) / stdev + 0.5 * stdev;
        double d2 = d1 - stdev;
        double nd1 = cnorm(d1);
        double nd2 = cnorm(d2);
        double u = fwd * nd1 - strike() * nd2;
        // if it is a put we use put-call parity:  F-K = Call-Put => Put = Call - (F-K)
        return isCall() ? u : u - (fwd - strike());
    }
};

class EuropeanDigital : public EuropeanOption
{
public:
    using EuropeanOption::EuropeanOption; // resuse constructor of new class

    // implements virtual method declared in the base class
    double undiscounted(double fwd, double stdev) const override
    {
        double d2 = log(fwd / strike()) / stdev - 0.5 * stdev;
        double nd2 = cnorm(d2);

        // if it is a put we use digital put-call parity: P(S<K)+P(S>K)=1
        return isCall() ? nd2 : 1 - nd2;
    }
};


int main()
{
    // alternative syntax for a defining a new type (typedef)
    // shared_ptr is a particular type of smart pointer defined in the STL
    using europtr = shared_ptr<EuropeanOption>;

    vector<europtr> trades;
    // add various types of trades to the vector
    trades.push_back(europtr(new EuropeanVanilla(100, 2, 2.001, 10.1, true)));
    trades.push_back(europtr(new EuropeanVanilla(100, 2, 2.001, 9.9, false)));
    trades.push_back(europtr(new EuropeanDigital(100, 2, 2.001, 9.9, true)));
    trades.push_back(europtr(new EuropeanDigital(100, 2, 2.001, 9.9, false)));

    // compute the price for all options
    for (const auto& t : trades)
        cout << t->price(0, 10, 0.05, 0.30) << endl;

    return 0;
}

