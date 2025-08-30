#include <iostream>
#include <cmath>

int main()
{
    int n = 9;

    // check if n is prime
    int divisor = 0;
    if (n > 2 && n % 2 == 0) {
        divisor = 2;
    }
    else {
        for (int i = 3, m = (int)std::sqrt(n); i <= m; i += 2) {
            if (n % i == 0) { // is n divisible by i ?
                divisor = i;
                break;
            }
        }
    }

    std::cout << "the number " << n << " is ";
    if (divisor==0)
        std::cout << "prime\n";
    else {
        std::cout << "divisible by " << divisor << "\n";
    }
    return 0;
}
