#include <iostream>
#include <utility>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

struct VolMark
{
    double strike;
    double volatility;
};

// cumulative standard normal distribution function
double cnorm(double x)
{
    return 0.5 * erfc(x / (-sqrt(2.0)));
}

// compute the undiscounted price of a call option using Black formula
double ucall(double fwd, double strike, double vol, double expiry)
{
    double stdev = vol * sqrt(expiry);
    double d1 = log(fwd / strike) / stdev + 0.5 * stdev;
    double d2 = d1 - stdev;
    double nd1 = cnorm(d1);
    double nd2 = cnorm(d2);
    return fwd * nd1 - strike * nd2;
}

// return false if it find an arbitrage, true otherwise
// Efficient implementation: it does not allocate temporary memory using pOld and pNew
bool check_for_arbitrages(const vector<VolMark>& marks, double fwd, double expiry)
{
    struct StrikeAndPrice
    {
        double strike;
        double price;
        double slope;
    };

    StrikeAndPrice pOld = {0, fwd, -1.0};
    for (const auto& [strike, volatility] : marks) {
        StrikeAndPrice pNew;
        pNew.strike = strike;
        pNew.price = ucall(fwd, strike, volatility, expiry);
        pNew.slope = (pNew.price - pOld.price) / (pNew.strike - pOld.strike);
        if (pOld.strike >= pNew.strike) // strikes must be in increasing order
            return false;
        if (pOld.slope >= pNew.slope)   // slopes must be in increasing order
            return false;

        pOld = pNew; // corkscrew
    }
    if (pOld.slope >= 0) // last slope must be < 0
        return false;

    return true;
}

int main()
{
    double fwd = 1.0;
    double expiry = 1.0;
    vector<VolMark> marks =
        { { 0.412605041576915,	0.35 }
        , { 0.627556851754149,	0.32 }
        , { 1.0,                0.30 }
        , { 1.59348113045318,	0.31 } // try with vol=0.61
        , { 2.4236252567766,	0.32 }
        };

    bool ok = check_for_arbitrages(marks, fwd, expiry);
    cout << "the marks are " << (ok ? "" : "not ") << "valid" << endl;

    return 0;
}

