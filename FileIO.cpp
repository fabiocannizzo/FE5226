#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

int main()
{
    // write
    {
        std::ofstream os("sqrt.txt");
        os << std::setprecision(3) << std::fixed;
        for (size_t i = 1; i < 10; ++i)
            os << std::setw(2) << i << ": "
            << std::setw(10) << std::sqrt(i) << '\n';
    }

    // read
    {
        std::ifstream is("sqrt.txt");

        std::cout << std::setprecision(3) << std::fixed;

        while (!is.eof()) {
            int i; double y; char c;

            is >> i >> c >> y;

            if (is.fail())
                break;

            std::cout << std::setw(2) << i << ": "
                << std::setw(10) << y << '\n';
        }
    }

    return 0;
}
