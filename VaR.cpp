#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

using namespace std;

int main()
{
    vector<double> pl(1000);

    // generate portfolio P&L with random numbers
    // let's pretend these are portfolio P&L in different possible future scenarios
    generate(pl.begin(), pl.end(), []() { return (rand() % 201) - 100;} );


    // We want to find what is our maximum loss with 95% confidence

    // We could sort the vector, then read the 50-th element,
    // but that would be an overkill, beacuse we do not need to sort all vector
    // We just need a partial sorting so that all elements before element 50
    // are smaller than element 50 and all elements after element 50 are larger
    // than elemenet 50

    // iterator pointing to the 50-th element in the vector
    auto element50 = pl.begin() + 50;
    // partial sort
    nth_element(pl.begin(), element50, pl.end());
    cout << "The portfolio VaR is: "
         << *element50
         << endl;
}
