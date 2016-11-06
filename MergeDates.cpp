#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>

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

void print( const vector<int>& v )
{
    for (size_t i = 0; i < v.size(); ++i)
        cout << v[i] << ", ";
    cout << endl;
}

int main()
{
    vector<int> v1(10);
    vector<int> v2(15);
    print(v1);
    print(v2);

    // we do not know in advance the size of the results,
    // but we know it cannot be big the the sum of the two sizes 
    vector<int> res(v1.size()+v2.size());

    // we start with begin()+1, because we are happy for the first element to remain zero
    generate(v1.begin()+1, v1.end(), IncreasingDateGenerator());
    generate(v2.begin()+1, v2.end(), IncreasingDateGenerator());
    print(v1);
    print(v2);

    // compute the union
    auto it = set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin());
    res.resize( distance(res.begin(), it) );
    print(res);
}
