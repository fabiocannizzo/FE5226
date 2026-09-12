#include <iostream>
#include <cmath>

using namespace std;

// Generic root search using bisection and a callback
double rootSearch(double a, double b, double (*fun)(double)) {
    double fa = (*fun)(a);
    double c = a;
    while ((b - a) > 1e-6) {
        c = 0.5 * (a + b);
        double fc = (*fun)(c);
        if (fc == 0.0)
            return c;
        if (fa * fc < 0.0)
            b = c;
        else {
            a = c;
            fa = fc;
        }
    }
    return 0.5 * (a + b);
}

// Target function: y = x^2 - 1
double myFunc(double x) {
    return x * x - 1.0;
}

int main() {
    // Find root of x^2 - 1 = 0 in [0, 2]
    double root = rootSearch(0.0, 2.0, &myFunc);
    cout << "Root of x^2 - 1 in [0, 2]: " << root << endl;
    return 0;
}
