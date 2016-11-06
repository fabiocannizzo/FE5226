#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <functional>

using namespace std;

// generate integer random number in the range [a,b)
// assumes that b>a
int myrand(int a, int b)
{
    int div = b-a;
    return a + (rand() % (b-a));
}

// accepts any container, not only vectors
template <typename Container>
void print( const Container& v )
{
    for ( auto i = v.begin(); i != v.end(); ++i)
        cout << *i << ", ";
    cout << endl;
}

int main()
{
    vector<int> v(20);

    // fill vector with random numbers in the range [1,16), i.e. [1,15]
    generate( v.begin(), v.end(), bind(myrand, 1, 16) );
    print(v);

    // copy all elements into as set, which will sort and remove duplicates
    set<int> tmp(v.begin(), v.end());
    print(v);  // note we use the same print function as per a vector, thanks to templates
    
    // copy back to the vector 
    v.assign(tmp.begin(), tmp.end());
    print(v);
}