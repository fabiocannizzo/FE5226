#include <vector>

// generic pointer based function
template <typename F>
void binaryVectorOp(double *dest, const double *src1, const double *src2, size_t n, F&& f)
{
    std::tranform(src1, src1+n, src2, dest, f);
}


// pointer based add function
template <typename F>
void vectorAdd(double *dest, const double *src1, const double *src2, size_t n)
{
    binaryVectorOp(dest, src1, src2, n, std::plus<double>{});
}

// vector based add function
template <typename F>
void vectorAdd(std::vector<double>& dest, const std::vector<double>& src1, const std::vector<double>& src2)
{
    // FIXME: here we should assert that the three vectros have the same size
    vectorAdd(dest.data(), src1.data(), src2.data(), src1.size());
}

// ... other functions



