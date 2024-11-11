#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;


auto evalCallOption(double ST, double k)
{ 
    return max(ST-k, 0.0);
}

struct Mkt
{
    double S0;
    double vol;
    double mu;
    double ir;
};

template <typename F>
double mcEval(const Mkt& mkt, double T, F&& f)
{
    const size_t nSim = 100000;

    // Create a Mersenne Twister random number generator with a seed
    std::mt19937 rng(std::random_device{}());

    // Define a normal distribution with a mean of 0.0 and standard deviation of 1.0
    std::normal_distribution<double> distribution(0.0, 1.0);

    const double df = exp(-mkt.ir*T);
    const double drift = mkt.mu * T;
    const double stdev = mkt.vol * sqrt(T);

    double sum = 0;
    for (size_t i = 0; i < nSim; ++i) {
        double eps = distribution(rng);  // generate gaussian in N(0,1)
        double ST = mkt.S0*exp(drift + stdev*eps);
        double po = f(ST);
        sum += po;
    }

    return df * (sum / nSim);
}

int main()
{
    auto callWithStrike12 = [](double ST) {return evalCallOption(ST, 1.02);};

    std::cout << mcEval(Mkt{1.0, 0.3, 0.05, 0.1}, 2.0, callWithStrike12) << "\n";

    return 0;
}
