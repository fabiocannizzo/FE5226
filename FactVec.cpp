#include <iostream>
#include <exception>
//#include <vector>

//#include "GetInteger.h"

using namespace std;

struct A {};
struct B {};
struct C {};
void foo(int x)
{
    if (x > 10)
        throw A{};
    else if (x < -10)
        throw B{};
    else
        throw C{};
}
void bar(int x)
{
    try {
        foo(x);
    }
    catch (A&) { cout << "A exception"; } // catch A exceptions, others pass through
}
int main()
{
    std::cout << "running\n";
    try { bar(20); }
    catch (B&) {
        cout << "B exception";
    } // catch B exceptions, others pass through
    catch (...) {
        cout << "Some other exception";
    } // catch everything else. No special handler for C
    return 0;
}