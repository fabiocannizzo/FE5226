#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdint>
using namespace std;

int main()
{
    int A[20];

    // fill up with ranfom values in [0..9]
    for (auto& a : A)
        a = rand() % 10;

    // print values
    for (auto a : A)
        std::cout << a << ',';
    cout << '\n';

    // compute the sum using a while-loop based on pointer arithmetic
    {
        const auto *begin = A;
        const auto *end = begin + sizeof(A)/sizeof(*A);
        double s = 0;
        while(begin != end)
            s += *begin++; // get the value pointed by begin, then increments the pointer
        cout << "sum: " << s << '\n';
    }

    // compute the sum using a for-loop based on pointer arithmetic
    {
        double s = 0;
        for(const auto *begin = A, *end = begin + sizeof(A)/sizeof(*A); begin != end; ++begin)
            s += *begin; // get the value pointed by begin
        cout << "sum: " << s << '\n';
    }


    cout << endl;
}

