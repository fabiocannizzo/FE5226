#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>

using namespace std;

// generate integer random number in the range [a,b)
// assumes that b>a
int myrand(int a, int b)
{
    return a + (rand() % (b-a));
}

void print( const vector<int>& v )
{
    for (auto x :  v)
        cout << x << ", ";
    cout << endl;
}

int main()
{
    vector<int> v(20);

    // fill vector with random numbers in the range [1,16), i.e. [1,15]
    generate( v.begin(), v.end(), bind(myrand, 1, 16) );
    //generate( v.begin(), v.end(), [](){return 1 + (rand() % 15)} );
    print(v);

    sort(v.begin(), v.end());
    print(v);

    // remove duplicate elements (the function unique operates on sorted vectors)
    auto it = unique(v.begin(), v.end());
    print(v);
    v.erase( it, v.end() );
    print(v);
}
