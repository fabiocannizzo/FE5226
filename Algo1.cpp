#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>

using namespace std;

// generate integer random number in the range [a,b)
// assumes that b>a
int myrand(int a, int b)
{
    int div = b-a;
    return a + (rand() % (b-a));
}

void print( const vector<int>& v )
{
    for (size_t i = 0; i < v.size(); ++i)
        cout << v[i] << ", ";
    cout << endl;
}

int main()
{
    vector<int> v(20);

    // fill vector with random numbers in the range [1,16), i.e. [1,15]
    generate( v.begin(), v.end(), bind(myrand, 1, 16) );
    print(v);

    sort(v.begin(), v.end());
    print(v);

    auto it = unique(v.begin(), v.end());
    print(v);

    size_t newsize = distance(v.begin(), it);
    v.resize( newsize );
    print(v);
}
