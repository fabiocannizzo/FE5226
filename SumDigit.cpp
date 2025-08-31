#include <iostream>

int main()
{
    int n = 3481;
    int nsave = n;

    int sumDigits = 0;

    while(n) {
        auto digit = n % 10;
        std::cout << "lowest digit: " << digit << "\n";
        std::cout << "residual:  (" << n << "-" << digit << ")/10 : " << (n /= 10) << "\n";

        sumDigits += digit;
    }

    std::cout << "sum of digits for " << nsave << " is " << sumDigits << "\n";

    return 0;
}

