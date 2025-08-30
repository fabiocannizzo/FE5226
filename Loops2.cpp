#include <iostream>

int main()
{
    // compute the product of (n^3-1) / (n^3+1) for n=2 to 1000
    double product = 1.0;
    for (int n = 2; n <= 1000; ++n) {
        // Note that if n*n*n > 2^31-1, then n3 overflows!
        // To be safe we could write: auto n3 = static_cast<double>(n) * n * n;
        // Converting the first n to double, then we have 2 multiplications of type double
        // With n=1000 we are still ok, and we can use int
        auto n3 = n * n * n;
        auto numerator = n3 - 1;
        auto denominator = n3 + 1;
        if (n < 10) {
            std::cout << numerator << "/" << denominator << "\n";
        }
        product *= static_cast<double>(numerator) / denominator; // remember to cast to double here, we want a division of type double
    }

    std::cout << "The product is: " << product << "\n";

    return 0;
}
