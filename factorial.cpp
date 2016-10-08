#include <iostream>
#include <cmath>
using namespace std;

unsigned factorial1(unsigned n)
{
    if (n > 1)
        return n*factorial1(n - 1);
    return 1;
}

unsigned factorial2(unsigned n)
{
    unsigned result = 1;
    for (; n>1; --n)
        result *= n;
    return result;
}


int main()
{
    cout << factorial1(4) << "\n";
    cout << factorial2(4) << "\n";
}

