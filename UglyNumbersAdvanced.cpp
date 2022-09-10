#include <iostream>
#include <algorithm>
using namespace std;

/*
    Ugly numbers are numbers whose only prime factors are 2, 3 or 5.
    The sequence:
        1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, ...
    shows the first 11 ugly numbers.
    By convention, 1 is included.
    Given a number n, the task is to find the n-th Ugly number.
*/

// see explanation for this algorithm here:
// https://www.geeksforgeeks.org/ugly-numbers/  (method 2)
size_t getNthUglyNo(size_t n)
{
    // allocate memory
    size_t* ugly = new size_t[n];

    ugly[0] = 1; // the first ugly number is 1

    // ugly numbers are only dividible by these factors
    size_t factors[3] = { 2, 3, 5 };

    // these are the next ugly numbers candidate, one in correspondence of each factor
    size_t candidates[3] = { 2, 3, 5 };

    // indices of the next ugly number which will be used as a base for the calculation of the
    // next ugly number in correspondence of each factor
    size_t indices[3] = { 0, 0, 0 };

    // find and store the next n-1 ugly numbers
    for (size_t i = 1, j2 = 0, j3 = 0, j4 = 0; i < n; ) {

        // find the pointer to the smallest candidate
        size_t* smallest_candidate_ptr = min_element(candidates, candidates +3);

        // add the next ugly number to the array
        if (ugly[i - 1] < *smallest_candidate_ptr) { // there may be duplicates
            ugly[i] = *smallest_candidate_ptr;
            //std::cout << ugly[i] << ", ";
            ++i;
        }

        // use pointer arithmetic to determine the index of the selected candidate
        size_t indexOfUsedCandidate = smallest_candidate_ptr - candidates; // pointer arithmetic

        // increase the ugly number index associated with the selected factor
        auto index = ++indices[indexOfUsedCandidate];

        // update the selected candidate with the next ugly number obtainable with the corresponding factor
        candidates[indexOfUsedCandidate] = factors[indexOfUsedCandidate] * ugly[index];
    }

    size_t result = ugly[n - 1];

    // release memory
    delete[] ugly;

    return result;
}

int main()
{
    size_t n;
    cout << "Which is the index (starting from 1) of the ugly number you want to find? ";
    cin >> n;
    if (cin.fail() || n == 0) {
        cout << "invalid index entered\n";
        return -1; // terminate the program with an error code
    }

    if (n == 0)
        return 0;

    cout << "The " << n << "-th ugly number is " << getNthUglyNo(n) << "\n";

    return 0;
}

