#include "mini_matlab/matrix.hpp"

#include <limits>
#include <string>
#include <utility>

#ifdef MINI_MATLAB_TRACE_SPECIAL
#include <iostream>
#endif

namespace mini_matlab {
namespace {

std::string shape(std::size_t r, std::size_t c) {
    return std::to_string(r) + "x" + std::to_string(c);
}

std::string shape(const Matrix& m) {
    return shape(m.rows(), m.cols());
}

// rows * cols can wrap around size_t. Wrapping would produce a small,
// successful allocation and an object whose invariant is a lie, so check
// before allocating anything.
std::size_t checked_area(std::size_t rows, std::size_t cols) {
    if (rows != 0 && cols > std::numeric_limits<std::size_t>::max() / rows) {
        throw DimensionMismatch("dimensions " + shape(rows, cols) + " overflow size_t");
    }
    return rows * cols;
}

} // namespace

Matrix::Matrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), data_(checked_area(rows, cols)) {}

Matrix::Matrix(std::size_t rows, std::size_t cols, std::vector<double> values)
    : rows_(rows), cols_(cols), data_(std::move(values)) {
    if (data_.size() != checked_area(rows, cols)) {
        throw DimensionMismatch("expected " + std::to_string(checked_area(rows, cols)) +
                                " values for a " + shape(rows, cols) + " matrix, got " +
                                std::to_string(data_.size()));
    }
}

#ifdef MINI_MATLAB_TRACE_SPECIAL

Matrix::Matrix(const Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {
    std::cerr << "[Matrix] copy ctor    " << shape(*this)
              << "  (allocates " << data_.size() << " doubles)\n";
}

Matrix::Matrix(Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), data_(std::move(other.data_)) {
    std::cerr << "[Matrix] move ctor    " << shape(*this)
              << "  (steals buffer; source left claiming " << shape(other)
              << " with " << other.data_.size() << " elements)\n";
}

Matrix& Matrix::operator=(const Matrix& other) {
    std::cerr << "[Matrix] copy assign  " << shape(other) << "\n";
    rows_ = other.rows_;
    cols_ = other.cols_;
    data_ = other.data_;
    return *this;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept {
    std::cerr << "[Matrix] move assign  " << shape(other) << "\n";
    rows_ = other.rows_;
    cols_ = other.cols_;
    data_ = std::move(other.data_);
    return *this;
}

Matrix::~Matrix() {
    std::cerr << "[Matrix] dtor         " << shape(*this) << "\n";
}

#endif // MINI_MATLAB_TRACE_SPECIAL

double& Matrix::at(std::size_t r, std::size_t c) {
    if (r >= rows_ || c >= cols_) {
        throw OutOfRange("index (" + std::to_string(r) + "," + std::to_string(c) +
                         ") out of range for " + shape(*this) + " matrix");
    }
    return data_[r * cols_ + c];
}

double Matrix::at(std::size_t r, std::size_t c) const {
    if (r >= rows_ || c >= cols_) {
        throw OutOfRange("index (" + std::to_string(r) + "," + std::to_string(c) +
                         ") out of range for " + shape(*this) + " matrix");
    }
    return data_[r * cols_ + c];
}

Matrix& Matrix::operator+=(const Matrix& rhs) {
    if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
        throw DimensionMismatch("cannot add " + shape(*this) + " and " + shape(rhs));
    }
    for (std::size_t i = 0; i < data_.size(); ++i) {
        data_[i] += rhs.data_[i];
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
    if (rows_ != rhs.rows_ || cols_ != rhs.cols_) {
        throw DimensionMismatch("cannot subtract " + shape(rhs) + " from " + shape(*this));
    }
    for (std::size_t i = 0; i < data_.size(); ++i) {
        data_[i] -= rhs.data_[i];
    }
    return *this;
}

Matrix& Matrix::operator*=(double s) noexcept {
    for (double& x : data_) {
        x *= s;
    }
    return *this;
}

Matrix& Matrix::operator/=(double s) {
    if (s==0.0) {
        throw DivisionByZero("division by zero");
    }
    for (double& x : data_) {
        x /= s;
    }
    return *this;
}

Matrix Matrix::transpose() const {
    Matrix out(cols_, rows_);
    for (std::size_t r = 0; r < rows_; ++r) {
        for (std::size_t c = 0; c < cols_; ++c) {
            out(c, r) = (*this)(r, c);
        }
    }
    return out;
}

Matrix operator+(Matrix lhs, const Matrix& rhs) {
    lhs += rhs;
    return lhs;
}

Matrix operator-(Matrix lhs, const Matrix& rhs) {
    lhs -= rhs;
    return lhs;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.cols() != rhs.rows()) {
        throw DimensionMismatch("cannot multiply " + shape(lhs) + " by " + shape(rhs));
    }

    Matrix out(lhs.rows(), rhs.cols());

    // i-k-j rather than the textbook i-j-k. Both are O(m*n*p) in arithmetic,
    // but this order walks rhs and out along rows in the innermost loop, so
    // every access is a sequential stride through a contiguous buffer. i-j-k
    // walks rhs down a column, striding by rhs.cols() doubles and touching a
    // new cache line on nearly every step.
    for (std::size_t i = 0; i < lhs.rows(); ++i) {
        for (std::size_t k = 0; k < lhs.cols(); ++k) {
            const double a = lhs(i, k);
            for (std::size_t j = 0; j < rhs.cols(); ++j) {
                out(i, j) += a * rhs(k, j);
            }
        }
    }
    return out;
}

Matrix operator*(Matrix m, double s) {
    m *= s;
    return m;
}

Matrix operator*(double s, Matrix m) {
    m *= s;
    return m;
}

Matrix operator/(Matrix m, double s) {
    m /= s;
    return m;
}
bool operator==(const Matrix& a, const Matrix& b) noexcept {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a.data()[i] != b.data()[i]) {
            return false;
        }
    }
    return true;
}

bool operator!=(const Matrix& a, const Matrix& b) noexcept {
    return !(a == b);
}

} // namespace mini_matlab
