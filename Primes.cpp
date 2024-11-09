#include <iostream>
#include <cmath>
using namespace std;

int main()
{
    // we want to find all primes smaller than this upperBound
    unsigned upperBound = 20;

    // how many primes we found so far
    unsigned nPrimesFound = 0;

    // as we find new primes we store them here (if we find more than 80000, we tell the user)
    unsigned primes[80000]; // over-dimensioned

    // we check all prime candidates in [2,upperBound]
    for (unsigned i = 2; i < upperBound; ++i) {
        // there is no need to check divisors greater than sqrt(i)
        unsigned largest = static_cast<unsigned>(sqrt(static_cast<double>(i)));

        bool isPrime = true; // postulate the number is prime
        // if it can be divided by any of the prime number found so far, than it is not prime
        for (unsigned j = 0; j < nPrimesFound; ++j) {
            if (primes[j] > largest) // all remaining primes are >sqrt(n), so there is no need to check
                break; // exit from the loop in j
            if (i % primes[j] == 0) { // divisible by primes[j]?
                isPrime = false;
                break; // exit from the loop in j
            }
        }

        if (isPrime) {
            if (nPrimesFound == sizeof(primes) / sizeof(*primes)) {
                cout << "Excedeed maximum storage capacity.\n"
                     << "The largest prime numbers will be missing.\n"
                     << "Use a larger array.\n";
                break;
            }
            // add i to the array and increment nPrimesFound
            primes[nPrimesFound++] = i;
        }
    }

    // print prime numbers found
    for (unsigned i = 0; i < nPrimesFound; ++i)
        cout << primes[i] << ", ";
    cout << endl;

    return 0;
}
