#include <iostream>
#include <vector>

#include "GetInteger.h"

using namespace std;

class FactorialCalculator
{
    vector<unsigned long> m_cache;
public:
    FactorialCalculator()
        : m_cache(1,1)  // init with 0!
    {
    }

    unsigned long operator()(unsigned long n)
    {
        if (n + 1 > m_cache.size()) {
            for (size_t i = m_cache.size(); i <= n; ++i)
                m_cache.push_back(i * m_cache.back());
        }
        return m_cache[n];
    }
};

int main()
{
    FactorialCalculator f;
    while(true)
    {
         int n = getInteger("Enter a positive integer (-1 to terminate): ");
         if (n < 0)
             break;
         cout << "Factorial of " << n << " is " << f(n) << endl;
    }
    return 0;
}
