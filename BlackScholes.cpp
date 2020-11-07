#include <cmath>
#include <iostream>

using namespace std;

// cumulative standard normal distribution function
double cnorm(double x)
{
    return 0.5 * erfc(x / (-sqrt(2.0)));
}

int main()
{
    double vol, T, K, F, r;

    cout << "Enter vol: ";
    cin >> vol;

    cout << "Enter time to maturity in years: ";
    cin >> T;

    cout << "Enter forward price: ";
    cin >> F;

    cout << "Enter continuously compounded interest rate: ";
    cin >> r;

    cout << "Enter strike: ";
    cin >> K;

    double stdev = vol * sqrt(T);
    double d1 = log( F / K) / stdev + 0.5 * stdev;
    double d2 = d1 - stdev;

    // cumnorm(x) = 0.5 + 0.5 * erf(x/sqrt(2))
    double nd1 = cnorm(d1);
    double nd2 = cnorm(d2);

    double disc_fact = exp(-r * T);

    double price = disc_fact * (F * nd1 - K * nd2);

    cout << "The price of the European Call ooption is: " << price << endl;
}
