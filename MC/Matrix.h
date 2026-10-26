#pragma once

#include <cmath>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define MYASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::ostringstream os; \
            os << msg; \
            throw std::invalid_argument(os.str()); \
        } \
    } while (0)

namespace mc {

class Matrix;
class MatrixView;

// Read-only view into a contiguous block of rows of a Matrix
class ConstMatrixView {
public:
    ConstMatrixView(const double* data, size_t rows, size_t cols)
        : m_data(data), m_rows(rows), m_cols(cols) {}

    ConstMatrixView(const Matrix& m);
    ConstMatrixView(const MatrixView& mv);

    size_t        rows() const { return m_rows; }
    size_t        cols() const { return m_cols; }
    const double* data() const { return m_data; }

    double operator()(size_t r, size_t c) const {
        return m_data[r * m_cols + c];
    }

    ConstMatrixView rowSlice(size_t r_start, size_t r_end) const {
        MYASSERT(r_start <= r_end && r_end <= m_rows, "Invalid rowSlice range");
        return ConstMatrixView(m_data + r_start * m_cols, r_end - r_start, m_cols);
    }

private:
    const double* m_data;
    size_t        m_rows;
    size_t        m_cols;
};

// Read-write view into a contiguous block of rows of a Matrix
class MatrixView {
public:
    MatrixView(double* data, size_t rows, size_t cols)
        : m_data(data), m_rows(rows), m_cols(cols) {}

    MatrixView(Matrix& m);

    size_t        rows() const { return m_rows; }
    size_t        cols() const { return m_cols; }
    double*       data()       { return m_data; }
    const double* data() const { return m_data; }

    double& operator()(size_t r, size_t c) {
        return m_data[r * m_cols + c];
    }

    double operator()(size_t r, size_t c) const {
        return m_data[r * m_cols + c];
    }

    MatrixView rowSlice(size_t r_start, size_t r_end) {
        MYASSERT(r_start <= r_end && r_end <= m_rows, "Invalid rowSlice range");
        return MatrixView(m_data + r_start * m_cols, r_end - r_start, m_cols);
    }

private:
    double* m_data;
    size_t  m_rows;
    size_t  m_cols;
};

// Contiguous row-major 2D matrix
class Matrix {
public:
    Matrix() : m_rows(0), m_cols(0) {}

    Matrix(size_t rows, size_t cols, double val = 0.0)
        : m_rows(rows), m_cols(cols), m_data(rows * cols, val) {}

    size_t        rows() const { return m_rows; }
    size_t        cols() const { return m_cols; }
    double*       data()       { return m_data.data(); }
    const double* data() const { return m_data.data(); }

    double& operator()(size_t r, size_t c) {
        return m_data[r * m_cols + c];
    }

    double operator()(size_t r, size_t c) const {
        return m_data[r * m_cols + c];
    }

    ConstMatrixView rowSlice(size_t r_start, size_t r_end) const {
        MYASSERT(r_start <= r_end && r_end <= m_rows, "Invalid rowSlice range");
        return ConstMatrixView(m_data.data() + r_start * m_cols, r_end - r_start, m_cols);
    }

    MatrixView rowSlice(size_t r_start, size_t r_end) {
        MYASSERT(r_start <= r_end && r_end <= m_rows, "Invalid rowSlice range");
        return MatrixView(m_data.data() + r_start * m_cols, r_end - r_start, m_cols);
    }

private:
    size_t              m_rows;
    size_t              m_cols;
    std::vector<double> m_data;
};

inline ConstMatrixView::ConstMatrixView(const Matrix& m)
    : m_data(m.data()), m_rows(m.rows()), m_cols(m.cols()) {}

inline ConstMatrixView::ConstMatrixView(const MatrixView& mv)
    : m_data(mv.data()), m_rows(mv.rows()), m_cols(mv.cols()) {}

inline MatrixView::MatrixView(Matrix& m)
    : m_data(m.data()), m_rows(m.rows()), m_cols(m.cols()) {}

// Cholesky decomposition: Sigma = L * L^T for symmetric positive-definite Sigma
inline Matrix cholesky(const Matrix& Sigma) {
    MYASSERT(Sigma.rows() == Sigma.cols(), "Correlation matrix must be square");
    size_t M = Sigma.rows();
    Matrix L(M, M, 0.0);
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double sum = Sigma(i, j);
            for (size_t k = 0; k < j; ++k)
                sum -= L(i, k) * L(j, k);
            if (i == j) {
                MYASSERT(sum > 0.0, "Correlation matrix is not symmetric positive-definite!");
                L(i, i) = std::sqrt(sum);
            } else {
                L(i, j) = sum / L(j, j);
            }
        }
    }
    return L;
}

} // namespace mc
