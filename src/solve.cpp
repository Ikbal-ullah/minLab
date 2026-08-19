#include "mini_matlab/solve.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace mini_matlab {
namespace {

std::string shape(const Matrix& m) {
    return std::to_string(m.rows()) + "x" + std::to_string(m.cols());
}

void swap_rows(Matrix& m, std::size_t r1, std::size_t r2) {
    if (r1 == r2) {
        return;
    }
    for (std::size_t j = 0; j < m.cols(); ++j) {
        std::swap(m(r1, j), m(r2, j));
    }
}

// A NaN or Inf anywhere in the system makes every downstream test meaningless:
// comparisons against NaN are all false, so a NaN would silently lose the pivot
// search and then be reported as "singular to working precision" -- a diagnosis
// that is simply untrue. Rejecting it here costs O(n^2) in front of an O(n^3)
// algorithm and lets the caller see what actually went wrong.
void require_finite(const Matrix& m, const char* what) {
    for (std::size_t i = 0; i < m.size(); ++i) {
        if (!std::isfinite(m.data()[i])) {
            throw MatrixError(std::string(what) +
                              " contains a non-finite value (NaN or Inf); the system has no "
                              "meaningful solution");
        }
    }
}

// The largest magnitude in the original matrix, used to make the singularity
// test scale-invariant. An absolute threshold would call a well-conditioned
// system singular merely because it was expressed in different units.
double max_magnitude(const Matrix& m) {
    double best = 0.0;
    for (std::size_t i = 0; i < m.size(); ++i) {
        const double v = std::fabs(m.data()[i]);
        if (v > best) {
            best = v;
        }
    }
    return best;
}

} // namespace

Matrix solve(Matrix a, Matrix b) {
    if (a.rows() != a.cols()) {
        throw DimensionMismatch("solve requires a square matrix, got " + shape(a));
    }
    if (a.rows() != b.rows()) {
        throw DimensionMismatch("cannot solve " + shape(a) + " against a right-hand side of " +
                                shape(b));
    }

    const std::size_t n = a.rows();
    const std::size_t m = b.cols();

    // A 0x0 system is vacuously solved by a 0-by-m matrix. Returning early also
    // keeps `n - 1` out of the code below, where it would wrap to SIZE_MAX.
    if (n == 0) {
        return b;
    }

    require_finite(a, "coefficient matrix");
    require_finite(b, "right-hand side");

    // Computed before elimination begins: the threshold must be relative to the
    // problem as posed, not to the partially eliminated matrix, whose entries
    // grow by the growth factor as elimination proceeds.
    const double tol = max_magnitude(a) * static_cast<double>(n) *
                       std::numeric_limits<double>::epsilon();

    // ---- forward elimination -------------------------------------------
    // k runs to n-1, not n-2. The final iteration performs no elimination
    // (the inner loops are empty) but does test the last pivot, which
    // back-substitution is about to divide by.
    for (std::size_t k = 0; k < n; ++k) {
        // Partial pivoting: the largest magnitude at or below the diagonal in
        // column k. This forces |multiplier| <= 1, which bounds the growth
        // factor and with it the backward error.
        // The scan runs from i = k seeded with 0.0, rather than from i = k+1
        // seeded with |a(k,k)|. Identical results, but this form never depends
        // on a comparison against the incumbent being meaningful -- which is
        // exactly the property the seeded form lacks, and why require_finite
        // above is not merely defensive tidiness. This is LAPACK's idamax.
        std::size_t pivot_row = k;
        double pivot_mag = 0.0;
        for (std::size_t i = k; i < n; ++i) {
            const double v = std::fabs(a(i, k));
            if (v > pivot_mag) {
                pivot_mag = v;
                pivot_row = i;
            }
        }

        if (pivot_mag <= tol) {
            throw SingularMatrix("matrix is singular to working precision (pivot " +
                                 std::to_string(pivot_mag) + " at column " + std::to_string(k) +
                                 " is below tolerance " + std::to_string(tol) + ")");
        }

        swap_rows(a, k, pivot_row);
        swap_rows(b, k, pivot_row);

        const double pivot = a(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = a(i, k) / pivot;

            // Assigned rather than computed: a(i,k) - factor*pivot is zero only
            // up to rounding. Writing the exact zero keeps the matrix genuinely
            // upper triangular, so the full-row swaps above stay correct and the
            // lower triangle cannot feed stale values back in.
            a(i, k) = 0.0;

            // j starts at k+1: column k was just set, and columns below k are
            // already zero from earlier steps.
            for (std::size_t j = k + 1; j < n; ++j) {
                a(i, j) -= factor * a(k, j);
            }
            for (std::size_t j = 0; j < m; ++j) {
                b(i, j) -= factor * b(k, j);
            }
        }
    }

    // ---- back substitution ---------------------------------------------
    Matrix x(n, m);
    for (std::size_t col = 0; col < m; ++col) {
        // `i-- > 0` rather than `i >= 0`: i is unsigned, so `i >= 0` is always
        // true and the loop never terminates. This form tests before
        // decrementing, so the body runs with i = n-1 down to 0 and exits when
        // the test sees 0.
        for (std::size_t i = n; i-- > 0;) {
            double sum = b(i, col);
            for (std::size_t j = i + 1; j < n; ++j) {
                sum -= a(i, j) * x(j, col);
            }
            x(i, col) = sum / a(i, i);
        }
    }

    return x;
}

} // namespace mini_matlab
