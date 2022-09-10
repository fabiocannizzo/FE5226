#include <iostream>
using namespace std;

/*
    Ugly numbers are numbers whose only prime factors are 2, 3 or 5.
    The sequence:
        1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, ...
    shows the first 11 ugly numbers.
    By convention, 1 is included.
    Given a number n, the task is to find the n-th Ugly number.
*/

// If n has m has a factor, even multiple times, remove the factor
// Note that n is passed as a reference
size_t removeFactor(size_t& n, size_t factor)
{
    while (n % factor == 0)  // is n divisible by factor?
        n /= factor;  // divide by factor
    // we also return the value, which is convenient to concatenate checks in the caller function
    return n;
}

bool isUgly(size_t n)
{
    // because of short boolean evaluation, so if removeFactor(n, 2) == 1, then the rest of
    // the expression is not executed
    return removeFactor(n, 2) == 1 || removeFactor(n, 3) == 1 || removeFactor(n, 5) == 1;
}

size_t getNthUglyNo(size_t n)
{
    size_t candidate = 0; // it is incremented at the begin of the loop, so the first candidate checked will be 1

    do {
        ++candidate;
        if (isUgly(candidate))
            --n;
    } while (n);

    return candidate;
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

