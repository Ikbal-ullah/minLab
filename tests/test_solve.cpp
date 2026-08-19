#include "test_harness.hpp"

#include "mini_matlab/matrix.hpp"
#include "mini_matlab/solve.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using mini_matlab::DimensionMismatch;
using mini_matlab::Matrix;
using mini_matlab::MatrixError;
using mini_matlab::SingularMatrix;
using mini_matlab::solve;

namespace {

// Infinity norm of a column vector / matrix: the largest magnitude anywhere.
double norm_inf(const Matrix& m) {
    double best = 0.0;
    for (std::size_t i = 0; i < m.size(); ++i) {
        best = std::fmax(best, std::fabs(m.data()[i]));
    }
    return best;
}

double residual(const Matrix& a, const Matrix& x, const Matrix& b) {
    return norm_inf(a * x - b);
}

void test_2x2_known_answer() {
    testharness::section("2x2 with a known answer");

    // 2x + 1y =  5
    // 1x + 3y = 10   ->  x = 1, y = 3
    Matrix a(2, 2, {2.0, 1.0, 1.0, 3.0});
    Matrix b(2, 1, {5.0, 10.0});

    Matrix x = solve(a, b);
    CHECK(x.rows() == 2 && x.cols() == 1);
    CHECK_NEAR(x(0, 0), 1.0, 1e-12);
    CHECK_NEAR(x(1, 0), 3.0, 1e-12);
    CHECK(residual(a, x, b) < 1e-12);

    // solve takes its arguments by value; the caller's matrices are untouched.
    CHECK(a == Matrix(2, 2, {2.0, 1.0, 1.0, 3.0}));
    CHECK(b == Matrix(2, 1, {5.0, 10.0}));

    // Solving twice must give the same answer -- the property an in-place
    // implementation would silently lose.
    CHECK(solve(a, b) == x);
}

void test_3x3_known_answer() {
    testharness::section("3x3 with a known answer");

    //  2x +  y -  z =   8
    // -3x -  y + 2z = -11   ->  x = 2, y = 3, z = -1
    // -2x +  y + 2z =  -3
    Matrix a(3, 3, {2.0, 1.0, -1.0, -3.0, -1.0, 2.0, -2.0, 1.0, 2.0});
    Matrix b(3, 1, {8.0, -11.0, -3.0});

    Matrix x = solve(a, b);
    CHECK(x.rows() == 3 && x.cols() == 1);
    CHECK_NEAR(x(0, 0), 2.0, 1e-12);
    CHECK_NEAR(x(1, 0), 3.0, 1e-12);
    CHECK_NEAR(x(2, 0), -1.0, 1e-12);
    CHECK(residual(a, x, b) < 1e-12);
}

void test_row_swap_on_first_pivot() {
    testharness::section("row swap required on the first pivot");

    // a(0,0) is exactly zero. Without pivoting the first multiplier is a
    // division by zero and the whole solve produces Inf/NaN.
    //   0x + 1y = 1
    //   1x + 0y = 2   ->  x = 2, y = 1
    Matrix a(2, 2, {0.0, 1.0, 1.0, 0.0});
    Matrix b(2, 1, {1.0, 2.0});

    Matrix x = solve(a, b);
    CHECK_NEAR(x(0, 0), 2.0, 1e-12);
    CHECK_NEAR(x(1, 0), 1.0, 1e-12);
    CHECK(residual(a, x, b) < 1e-12);

    // Same idea one size up, with the zero pivot in the middle column so the
    // swap happens at k = 1 rather than k = 0.
    Matrix a3(3, 3, {1.0, 2.0, 3.0, 2.0, 4.0, 7.0, 3.0, 5.0, 3.0});
    Matrix b3(3, 1, {14.0, 30.0, 22.0});
    Matrix x3 = solve(a3, b3);
    CHECK(residual(a3, x3, b3) < 1e-12);

    // The pathological case from the design notes. This matrix is
    // well-conditioned (kappa ~ 2.6), so any error here is the algorithm's
    // fault, not the problem's. Without pivoting the multiplier is 1e20, the
    // "1" in a(1,1) is absorbed below the last bit of it, and x(0,0) comes out
    // as 0 instead of 1.
    Matrix tiny(2, 2, {1e-20, 1.0, 1.0, 1.0});
    Matrix tb(2, 1, {1.0, 2.0});
    Matrix tx = solve(tiny, tb);
    CHECK_NEAR(tx(0, 0), 1.0, 1e-9);
    CHECK_NEAR(tx(1, 0), 1.0, 1e-9);
}

void test_exactly_singular() {
    testharness::section("exactly singular systems throw");

    // Row 1 is exactly 2x row 0. The second pivot eliminates to exactly 0.0.
    Matrix a(2, 2, {1.0, 2.0, 2.0, 4.0});
    CHECK_THROWS(SingularMatrix, solve(a, Matrix(2, 1, {3.0, 6.0})));  // consistent
    CHECK_THROWS(SingularMatrix, solve(a, Matrix(2, 1, {1.0, 1.0})));  // inconsistent
    CHECK_THROWS(MatrixError, solve(a, Matrix(2, 1, {3.0, 6.0})));     // REPL boundary

    // A zero column means no pivot exists anywhere in it.
    Matrix zero_col(2, 2, {0.0, 1.0, 0.0, 3.0});
    CHECK_THROWS(SingularMatrix, solve(zero_col, Matrix(2, 1, {1.0, 2.0})));

    // The all-zero matrix: max_magnitude is 0, so tol is 0, and the pivot 0.0
    // satisfies `0 <= 0`. Singular, as it must be.
    Matrix zeros(3, 3);
    CHECK_THROWS(SingularMatrix, solve(zeros, Matrix(3, 1, {1.0, 2.0, 3.0})));

    // Rank deficiency that only surfaces at the third pivot: row2 = row0 + row1.
    Matrix rank2(3, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 5.0, 7.0, 9.0});
    CHECK_THROWS(SingularMatrix, solve(rank2, Matrix(3, 1, {1.0, 2.0, 3.0})));

    // Scale invariance: the same singular matrix scaled by 1e12 is still
    // singular, and scaled by 1e-12 is still singular. An absolute pivot
    // tolerance would get one of these two wrong.
    CHECK_THROWS(SingularMatrix, solve(a * 1e12, Matrix(2, 1, {3.0, 6.0})));
    CHECK_THROWS(SingularMatrix, solve(a * 1e-12, Matrix(2, 1, {3.0, 6.0})));
}

void test_well_conditioned_but_tiny_scale() {
    testharness::section("scale invariance does not create false singulars");

    // Every entry is 1e-13 -- far below any plausible absolute tolerance --
    // but the system is perfectly well-conditioned and must solve cleanly.
    // This is the failure an absolute pivot threshold would produce.
    Matrix a(2, 2, {2e-13, 1e-13, 1e-13, 3e-13});
    Matrix b(2, 1, {5e-13, 1e-12});

    Matrix x = solve(a, b);
    CHECK_NEAR(x(0, 0), 1.0, 1e-6);
    CHECK_NEAR(x(1, 0), 3.0, 1e-6);
}

void test_ill_conditioned_residual() {
    testharness::section("ill-conditioned: assert on the residual, not on x");

    // Hilbert matrix, H(i,j) = 1/(i+j+1). kappa(H8) is roughly 1.5e10, so
    // double precision loses about ten of its sixteen digits on the *solution*.
    const std::size_t n = 8;
    Matrix h(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            h(i, j) = 1.0 / static_cast<double>(i + j + 1);
        }
    }

    // Manufacture a right-hand side from a known solution of all ones.
    Matrix x_true(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        x_true(i, 0) = 1.0;
    }
    Matrix b = h * x_true;

    Matrix x = solve(h, b);

    // Assert on the residual ||Hx - b||, NOT on ||x - x_true||.
    //
    // Partial pivoting makes the algorithm backward stable: the computed x is
    // the exact solution of a nearby system (H + dH)x = b with dH small. That
    // guarantee is about the residual, and it holds no matter how badly
    // conditioned H is.
    //
    // The forward error is a different quantity entirely and is bounded by
    // kappa(H) * (backward error). With kappa ~ 1.5e10 that permits x to be
    // wrong in the sixth decimal place while the residual stays at machine
    // precision. Asserting a tight bound on x would therefore be asserting
    // something the algorithm never promised -- it would be testing the
    // conditioning of Hilbert matrices, not the correctness of this code.
    const double res = residual(h, x, b);
    CHECK(res < 1e-12);

    // Printed rather than asserted: the exact figure depends on the compiler's
    // floating-point contraction and instruction selection, so pinning it would
    // make the test fragile. It is here to be looked at.
    double fwd = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        fwd = std::fmax(fwd, std::fabs(x(i, 0) - 1.0));
    }
    std::printf("   residual %.3e   forward error %.3e   (ratio %.1e)\n", res, fwd,
                fwd / (res > 0.0 ? res : 1.0));

    // The one thing worth asserting about x: it is not garbage. A loose bound
    // that documents "the answer is in the right neighbourhood" without
    // pretending to a precision the conditioning does not allow.
    CHECK(fwd < 1e-3);
}

void test_non_finite_input_throws() {
    testharness::section("non-finite input is rejected, not misreported");

    const double nan_v = std::nan("");
    const double inf_v = std::numeric_limits<double>::infinity();

    Matrix good_b(2, 1, {1.0, 2.0});

    // Every comparison against NaN is false, so without the up-front check a
    // NaN loses the pivot search, survives into the multiplier, and the solve
    // eventually reports "singular to working precision" -- which is not what
    // is wrong with this matrix. The error must name the actual problem.
    CHECK_THROWS(MatrixError, solve(Matrix(2, 2, {nan_v, 1.0, 1.0, 0.0}), good_b));
    CHECK_THROWS(MatrixError, solve(Matrix(2, 2, {1.0, nan_v, 1.0, 0.0}), good_b));
    CHECK_THROWS(MatrixError, solve(Matrix(2, 2, {inf_v, 1.0, 1.0, 0.0}), good_b));

    // Specifically NOT SingularMatrix -- the diagnosis must be distinguishable.
    bool reported_as_singular = false;
    try {
        solve(Matrix(2, 2, {nan_v, 1.0, 1.0, 0.0}), good_b);
    } catch (const SingularMatrix&) {
        reported_as_singular = true;
    } catch (const MatrixError&) {
        reported_as_singular = false;
    }
    CHECK(!reported_as_singular);

    // The right-hand side is checked too: a NaN there poisons back-substitution
    // just as thoroughly.
    Matrix good_a(2, 2, {2.0, 1.0, 1.0, 3.0});
    CHECK_THROWS(MatrixError, solve(good_a, Matrix(2, 1, {nan_v, 1.0})));
    CHECK_THROWS(MatrixError, solve(good_a, Matrix(2, 1, {1.0, inf_v})));

    // A finite system next to these must still solve.
    CHECK(residual(good_a, solve(good_a, good_b), good_b) < 1e-12);
}

void test_non_square_throws() {
    testharness::section("non-square and mismatched shapes throw");

    CHECK_THROWS(DimensionMismatch, solve(Matrix(2, 3), Matrix(2, 1, {1.0, 2.0})));
    CHECK_THROWS(DimensionMismatch, solve(Matrix(3, 2), Matrix(3, 1, {1.0, 2.0, 3.0})));

    // Square, but the right-hand side has the wrong number of rows.
    Matrix a(2, 2, {1.0, 0.0, 0.0, 1.0});
    CHECK_THROWS(DimensionMismatch, solve(a, Matrix(3, 1, {1.0, 2.0, 3.0})));
    CHECK_THROWS(DimensionMismatch, solve(a, Matrix(1, 1, {1.0})));

    // Shape is checked before singularity: a non-square zero matrix reports
    // the shape problem, which is the more actionable error.
    CHECK_THROWS(DimensionMismatch, solve(Matrix(2, 3), Matrix(2, 1)));

    CHECK_THROWS(MatrixError, solve(Matrix(2, 3), Matrix(2, 1)));  // REPL boundary
}

void test_identity() {
    testharness::section("identity");

    Matrix i3(3, 3, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
    Matrix b(3, 1, {7.0, -2.0, 0.5});

    // I \ b == b, exactly: every multiplier is 0 and every pivot is 1, so no
    // rounding occurs anywhere. Exact equality is the right assertion here.
    CHECK(solve(i3, b) == b);

    Matrix i1(1, 1, {1.0});
    CHECK(solve(i1, Matrix(1, 1, {42.0})) == Matrix(1, 1, {42.0}));

    // A diagonal matrix is the same story with divisions that do round.
    Matrix d(3, 3, {2.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 8.0});
    Matrix dx = solve(d, Matrix(3, 1, {2.0, 4.0, 8.0}));
    CHECK(dx == Matrix(3, 1, {1.0, 1.0, 1.0}));
}

void test_1x1_and_empty() {
    testharness::section("1x1 and empty systems");

    CHECK(solve(Matrix(1, 1, {4.0}), Matrix(1, 1, {12.0})) == Matrix(1, 1, {3.0}));
    CHECK_THROWS(SingularMatrix, solve(Matrix(1, 1, {0.0}), Matrix(1, 1, {1.0})));

    // A 0x0 system against a 0x1 right-hand side is vacuously solved.
    Matrix e = solve(Matrix(0, 0), Matrix(0, 1));
    CHECK(e.rows() == 0 && e.cols() == 1);
}

void test_multiple_right_hand_sides() {
    testharness::section("multiple right-hand sides");

    Matrix a(2, 2, {2.0, 1.0, 1.0, 3.0});

    // Two independent systems solved in one factorisation.
    Matrix b(2, 2, {5.0, 3.0, 10.0, 5.0});
    Matrix x = solve(a, b);
    CHECK(x.rows() == 2 && x.cols() == 2);
    CHECK(residual(a, x, b) < 1e-12);

    // Each column must match the single-RHS answer for that column.
    Matrix x0 = solve(a, Matrix(2, 1, {5.0, 10.0}));
    Matrix x1 = solve(a, Matrix(2, 1, {3.0, 5.0}));
    CHECK_NEAR(x(0, 0), x0(0, 0), 1e-15);
    CHECK_NEAR(x(1, 0), x0(1, 0), 1e-15);
    CHECK_NEAR(x(0, 1), x1(0, 0), 1e-15);
    CHECK_NEAR(x(1, 1), x1(1, 0), 1e-15);

    // Solving against the identity produces the inverse -- the expensive route
    // discussed in the design notes, here only as a correctness check.
    Matrix inv = solve(a, Matrix(2, 2, {1.0, 0.0, 0.0, 1.0}));
    CHECK(residual(a, inv, Matrix(2, 2, {1.0, 0.0, 0.0, 1.0})) < 1e-12);
}

void test_larger_system() {
    testharness::section("larger system, round trip");

    // A diagonally dominant 12x12 built deterministically: well-conditioned,
    // and large enough that every pivot search, swap and update path runs many
    // times. Round-trip through a known solution.
    const std::size_t n = 12;
    Matrix a(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            a(i, j) = 1.0 / static_cast<double>(1 + ((i * 7 + j * 3) % 11));
        }
        a(i, i) += 20.0;
    }

    Matrix x_true(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        x_true(i, 0) = static_cast<double>(i) - 5.5;
    }

    Matrix b = a * x_true;
    Matrix x = solve(a, b);

    CHECK(residual(a, x, b) < 1e-10);
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_NEAR(x(i, 0), x_true(i, 0), 1e-10);
    }
}

} // namespace

int main() {
    test_2x2_known_answer();
    test_3x3_known_answer();
    test_row_swap_on_first_pivot();
    test_exactly_singular();
    test_well_conditioned_but_tiny_scale();
    test_ill_conditioned_residual();
    test_non_finite_input_throws();
    test_non_square_throws();
    test_identity();
    test_1x1_and_empty();
    test_multiple_right_hand_sides();
    test_larger_system();

    return testharness::summary();
}
