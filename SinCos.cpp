/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program


 Output the values of sin(x)*cos(x) in the range 0-2pi at intervals equally spaced of size pi/10
 Justify the output to the right and print exactly 4 decimals
*/

// define is a compiler directive
// NOTE: if we were working with C++20 standard (or higher)
// we would not this. We would instead import the header <numbers>,
// which makes available various numerical constant
#define _USE_MATH_DEFINES  // cause import of some macros like M_PI

#include <iostream>
#include <iomanip>  // needed to format the output
#include <cmath>    // needed for sin, cos and M_PI

// import only the symbols needed
using namespace std;

int main()
{
    const unsigned nPoints = 10;
    const double step = M_PI / nPoints;  // here nPoints is lifted to double

    cout << std::fixed << std::setprecision(3)
    for (unsigned i = 0; i < 2 * nPoints; ++i) {
        double x = i*step;
        cout << std::setw(8) << x << ", "
             << std::setw(8) << sin(x)*cos(x)
             << "\n";
    }

    cout << "\n";
}

