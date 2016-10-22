#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

int main()
{
    const char *filename = "input.txt";
    const unsigned n = 10;

    ofstream of(filename);
    of << setprecision(3) << fixed;
    for (unsigned i = 0; i < n; ++i)
        of << i << " " << sqrt(i) << endl;
    of.close();

    ifstream ff(filename);
    while (!ff.eof()) {
        unsigned i;
        double s;
        ff >> i >> s;  // note that what we read has only 3 digits precision
        cout << i << " " << s << endl;
    }
    ff.close();

    return 0;
}
