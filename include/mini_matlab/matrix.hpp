#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace mini_matlab {

// One base type so the REPL can catch every engine error with a single clause.
class MatrixError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class DimensionMismatch : public MatrixError {
public:
    using MatrixError::MatrixError;
};

class OutOfRange : public MatrixError {
public:
    using MatrixError::MatrixError;
};

// A dense, row-major matrix of doubles.
//
// Invariant: data_.size() == rows_ * cols_.
//
// The invariant does NOT hold for a moved-from Matrix: data_ is emptied but
// rows_/cols_ are scalars and keep their old values. As with the standard
// library containers, the only operations supported on a moved-from Matrix are
// destruction and assignment.
class Matrix {
public:
    Matrix() = default;

    // Zero-filled.
    Matrix(std::size_t rows, std::size_t cols);

    // Row-major; values.size() must equal rows * cols.
    Matrix(std::size_t rows, std::size_t cols, std::vector<double> values);

#ifdef MINI_MATLAB_TRACE_SPECIAL
    // Present only under -DMINI_MATLAB_TRACE_SPECIAL, to make copy-vs-move
    // visible. These are exactly what the compiler would generate. Shipping
    // them unconditionally would be a liability: the user-declared destructor
    // alone suppresses the implicit move operations, silently downgrading every
    // move to a copy.
    Matrix(const Matrix& other);
    Matrix(Matrix&& other) noexcept;
    Matrix& operator=(const Matrix& other);
    Matrix& operator=(Matrix&& other) noexcept;
    ~Matrix();
#endif

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    std::size_t size() const noexcept { return data_.size(); }
    bool empty() const noexcept { return data_.empty(); }

    // Unchecked. Out-of-range indices are undefined behaviour, not an exception.
    double& operator()(std::size_t r, std::size_t c) noexcept;
    double operator()(std::size_t r, std::size_t c) const noexcept;

    // Bounds-checked; throws OutOfRange.
    double& at(std::size_t r, std::size_t c);
    double at(std::size_t r, std::size_t c) const;

    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator-=(const Matrix& rhs);
    Matrix& operator*=(double s) noexcept;

    Matrix transpose() const;

    const double* data() const noexcept { return data_.data(); }

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<double> data_;
};

// Defined inline: element access is the innermost operation in every loop in
// the engine, and a non-inlined call here would dominate the arithmetic.
inline double& Matrix::operator()(std::size_t r, std::size_t c) noexcept {
    return data_[r * cols_ + c];
}

inline double Matrix::operator()(std::size_t r, std::size_t c) const noexcept {
    return data_[r * cols_ + c];
}

// lhs by value: `a + b` copies once (same as any other form), but `f() + b`
// moves, and `a + b + c` reuses the first temporary instead of allocating again.
Matrix operator+(Matrix lhs, const Matrix& rhs);
Matrix operator-(Matrix lhs, const Matrix& rhs);

Matrix operator*(const Matrix& lhs, const Matrix& rhs);
Matrix operator*(Matrix m, double s);
Matrix operator*(double s, Matrix m);

// Exact comparison. NaN != NaN, and -0.0 == 0.0. Suitable for structural tests,
// not for comparing results of different floating-point computations.
bool operator==(const Matrix& a, const Matrix& b) noexcept;
bool operator!=(const Matrix& a, const Matrix& b) noexcept;

} // namespace mini_matlab
