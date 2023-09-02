#include <iostream>
#include <cmath>

int main()
{
    double series = 1.0;
    const double onethird = 1.0 / 3.0;
    double factor = onethird;
    for (size_t i = 1; i <= 100; ++i) {
        series += factor;
        factor *= onethird;
    }

    std::cout << "The sum of the first 101 terms of the series is: " << series << "\n";

    return 0;
}
