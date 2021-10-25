#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <set>

using namespace std;

// This is a 'stateful' generator
// It generate integer random number in the range [1-3]
// and adds it to the last date
struct IncreasingDateGenerator
{
    int m_t;
    IncreasingDateGenerator() : m_t(0) {}
    int operator()()
    {
        // update the generator state
        m_t += 1 + (rand() % 3);
        return m_t;
    }
};

void print(const vector<int>& v)
{
    for (size_t i = 0; i < v.size(); ++i)
        cout << v[i] << ", ";
    cout << endl;
}

// This preprocessor macro is used in the compiler directives below to
// decide which of the 4 alternative implementation we want to use
#define VERSION_TO_USE 1

int main()
{
    vector<int> v1(10);
    vector<int> v2(15);
    print(v1);
    print(v2);

    // we start with begin()+1, because we are happy for the first element to remain zero
    generate(v1.begin() + 1, v1.end(), IncreasingDateGenerator());
    generate(v2.begin() + 1, v2.end(), IncreasingDateGenerator());
    print(v1);
    print(v2);

#if VERSION_TO_USE==1
    // we do not know in advance the size of the results,
    // but we know it cannot be big the the sum of the two sizes
    // This method requires that v1 and v2 are sorted
    vector<int> res(v1.size() + v2.size());
    // compute the union
    auto it = set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin());  // O(n)
    it = unique(res.begin(), it);
    res.erase(it, res.end());
#elif VERSION_TO_USE==2
    // 2 operations with cost O(n logn) + memory allocation
    // This method does not require that v1 and v2 are sorted
    vector<int> res;
    set<int> temp(v1.begin(), v1.end());
    temp.insert(v2.begin(), v2.end());
    res.assign(temp.begin(), temp.end());
#else
    // we do not know in advance the size of the results,
    // but we know it cannot be big the the sum of the two sizes
    // This method does not require that v1 and v2 are sorted
    vector<int> res;
    res.reserve(v1.size() + v2.size());
    res = v1;
    res.insert(res.end(), v2.begin(), v2.end());
    sort(res.begin(), res.end());  // full sort operation O(n logn)
    auto it = unique(res.begin(), res.end());
    res.resize(std::distance(res.begin(), it));
#endif

    print(res);
}
