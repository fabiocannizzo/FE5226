/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program


 Output the values of sin(x)*cos(x) in the range 0-2pi at intervals equally spaced of size pi/10
 Justify the output to the right and print exactly 4 decimals
*/

#include <iostream>
#include <numbers>  // we use this for PI (requires c++20)
#include <iomanip>  // needed to format the output
#include <cmath>    // needed for sin, cos

// import only the symbols needed
using namespace std;

int main()
{
    const unsigned nPoints = 10;
    const double step = std::numbers::pi / nPoints;  // here nPoints is lifted to double

    cout << std::fixed << std::setprecision(3);
    for (unsigned i = 0; i < 2 * nPoints; ++i) {
        double x = i*step;
        cout << std::setw(8) << x << ", "
             << std::setw(8) << sin(x)*cos(x)
             << "\n";
    }
    cout << "\n";
}

