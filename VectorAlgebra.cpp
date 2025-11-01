#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

void print(const char *msg, const vector<double>& v)
{
    cout << msg << ": ";
    for (size_t i = 0; i < v.size(); ++i)
        cout << v[i] << ", ";
    cout << endl;
}


namespace Hidden {

//
// Template functions with generic function F
//

// generic pointer based function on two vectors
template <typename F>
void vecVecOp(double *dest, const double *src1, const double *src2, size_t n, F&& f)
{
    std::transform(src1, src1+n, src2, dest, f);
}

// generic pointer based function transformation on vector
template <typename F>
void vecOp(double *dest, const double *src, size_t n, F&& f)
{
    std::transform(src, src+n, dest, f);
}

} // namespace Hidden

//
// Actual toolkit functions
//

// vector based add function
void vectorAdd(std::vector<double>& dest, const std::vector<double>& src1, const std::vector<double>& src2)
{
    // FIXME: here we should assert that the three vectros have the same size
    Hidden::vecVecOp(dest.data(), src1.data(), src2.data(), src1.size(), std::plus<double>{});
}

// vector based add function
void vectorAdd(std::vector<double>& dest, const std::vector<double>& src1, double src2)
{
    // FIXME: here we should assert that the three vectros have the same size
    Hidden::vecOp(dest.data(), src1.data(), src1.size(), [src2](auto s2){return src2 + s2; } );
}

void vectorAdd(std::vector<double>& dest, double src1, const std::vector<double>& src2)
{
    // use the commutative property of addition
    vectorAdd(dest, src2, src1);
}

// ... other functions

int main()
{
    std::vector<double> v1 = {1,2,4,2};
    std::vector<double> v2 = {4,1,0,1};

    std::vector<double> v3(v1.size());
    vectorAdd(v3,v1,v2);

    print("v1", v1);
    print("v2", v2);
    print("v3", v3);
}
