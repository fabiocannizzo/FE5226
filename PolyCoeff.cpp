#include <iostream>
#include <cmath>

// coefficients of the polynomial: a^2 + bx + c
struct PolyCoef {
    double a, b, c;
};

// roots of the polynomial: a^2 + bx + c
struct Roots {
    double r1, r2;
};

// compute the roots for a polynomial: a^2 + bx + c
Roots roots(const PolyCoef& coef)
{
    auto &[a, b, c] = coef;  // create aliases for convenience to the data members in PolyCoef

    auto delta = b * b - 4.0 * a * c;
    auto sdelta = std::sqrt(delta);

    return { (-b + sdelta) /(2.0*a), (-b - sdelta) /(2.0*a)}; // initializer list for return type Roots
}

int main()
{
    // we use initializer list for the argument of typePolyCoef
    // auto r = roots(PolyCoef{1.0, -5.0, 6.0});
    // we can also omit the type name
    auto r = roots({1.0, -5.0, 6.0});
    std::cout << "roots: " << r.r1 << ", " << r.r2 << '\n';
    return 0;
}

