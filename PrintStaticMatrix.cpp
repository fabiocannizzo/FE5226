#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main()
{

    const int nr = 4;
    const int nc = 3;

    int M[nr *nc] =
        { 1, 2, 3
        , 4, 5 ,6
        , 7, 8, 9
        , 9, 10, 11
        };

    // row-major interpretation of the storage M
    for (int r = 0; r < nr; ++r) {
        for (int c = 0; c < nc; ++c)
            std::cout << (c>0 ? ", ": "") << std::setw(2) << M[r*nc + c];
        std::cout << '\n';
    }
    std::cout << '\n';

    return 0;
}

