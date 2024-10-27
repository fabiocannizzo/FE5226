#include <vector>

// generic pointer based function
template <typename F>
void binaryVectorOp(double *res, const double *src1, const double *src2, size_t n, F&& f)
{
    std::tranform(src1, src2+n, src2, res, f);
}


// pointer based add function
template <typename F>
void vectorAdd(double *res, const double *src1, const double *src2, size_t n)
{
    binaryVectorOp(res, src1, src2, n, std::multiplies<double>{});
}

// vector based add function
template <typename F>
void vectorAdd(std::vector<double>& res, const std::vector<double>& src1, const std::vector<double>& src2)
{
    // FIXME: here we should assert that the three vectros have the same size
    vectorAdd(res.data(), src1.data(), src2.data(), src1.size());
}

// ... other functions



