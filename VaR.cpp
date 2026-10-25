#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;

int main()
{
    vector<double> pl(100000);

    // generate portfolio P&L with random numbers
    // let's pretend these are portfolio P&L in different possible future scenarios
    generate(pl.begin(), pl.end(), []() { return (rand() % 201) - 100;} );

    // We want to find what is our maximum loss with 95% confidence

    // We could sort the vector, then read the 500-th element,
    // but that would be an overkill, beacuse we do not need to sort all vector
    // We just need a partial sorting so that all elements before element 500
    // are smaller than element 500 and all elements after element 501 are larger
    // than elemenet 500
    // 95% left-quantile of the vector is element number 500

    const size_t quantileIndex = pl.size() * 5 /100;

    //  
    // solve the problem using nth_element
    //

    {
        // create a copy, so that we do not modify the orignal vector
        vector<double> plcopy(pl);
        
        // iterator pointing to the 500-th element in the vector
        auto quantileIter = plcopy.begin() + quantileIndex;
    	
        // Record the starting time point
        auto start = std::chrono::high_resolution_clock::now();

        // partial sort
        nth_element(plcopy.begin(), quantileIter, plcopy.end());

        // Record the ending time point
        auto end = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        cout << "The portfolio VaR is: "
         << *quantileIter
         << ", computed with nth_element in " << duration.count() << " microseconds"
         << endl;
    }

    //
    // solve the problem using sort
    //

    {
    	// create a copy, so that we do not modify the orignal vector
        vector<double> plcopy(pl);
        
        // iterator pointing to the 500-th element in the vector
        auto quantileIter = plcopy.begin() + quantileIndex;
    	
        // Record the starting time point
        auto start = std::chrono::high_resolution_clock::now();

        // partial sort
        sort(plcopy.begin(), plcopy.end());

        // Record the ending time point
        auto end = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        cout << "The portfolio VaR is: "
            << *quantileIter
            << ", computed with std::sort in " << duration.count() << " microseconds"
            << endl;
    }
}