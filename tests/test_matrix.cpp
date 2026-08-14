#include "test_harness.hpp"

#include "mini_matlab/matrix.hpp"

#include <exception>
#include <utility>
#include <vector>

using mini_matlab::DimensionMismatch;
using mini_matlab::Matrix;
using mini_matlab::MatrixError;
using mini_matlab::OutOfRange;

namespace {

// Assigning through references defeats -Wself-assign-overloaded and -Wself-move,
// which only fire on the syntactically obvious `a = a` form. Real self-assignment
// bugs arrive through aliased references exactly like these, so this is both a
// workaround for -Werror and the more realistic test.
void copy_assign(Matrix& dst, const Matrix& src) { dst = src; }
void move_assign(Matrix& dst, Matrix& src) { dst = std::move(src); }

void test_construction() {
    testharness::section("construction and element access");

    Matrix z(2, 3);
    CHECK(z.rows() == 2);
    CHECK(z.cols() == 3);
    CHECK(z.size() == 6);
    CHECK(!z.empty());
    CHECK_NEAR(z(1, 2), 0.0, 0.0);  // vector value-initialises: elements are zero

    Matrix m(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    CHECK_NEAR(m(0, 0), 1.0, 1e-12);
    CHECK_NEAR(m(0, 2), 3.0, 1e-12);
    CHECK_NEAR(m(1, 0), 4.0, 1e-12);
    CHECK_NEAR(m(1, 2), 6.0, 1e-12);

    // Row-major: (r,c) lives at r*cols + c.
    CHECK_NEAR(m.data()[1 * 3 + 2], 6.0, 1e-12);

    m(0, 0) = 99.0;
    CHECK_NEAR(m(0, 0), 99.0, 1e-12);

    // Wrong number of values for the stated shape.
    CHECK_THROWS(DimensionMismatch, Matrix(2, 3, std::vector<double>{1.0, 2.0}));
    CHECK_THROWS(DimensionMismatch, Matrix(0, 0, std::vector<double>{1.0}));
}

void test_add_sub() {
    testharness::section("addition and subtraction");

    Matrix a(2, 2, {1.0, 2.0, 3.0, 4.0});
    Matrix b(2, 2, {10.0, 20.0, 30.0, 40.0});

    CHECK(a + b == Matrix(2, 2, {11.0, 22.0, 33.0, 44.0}));
    CHECK(b - a == Matrix(2, 2, {9.0, 18.0, 27.0, 36.0}));
    CHECK(a - a == Matrix(2, 2));

    // operator+ takes lhs by value; the caller's object must be untouched.
    Matrix sum = a + b;
    CHECK(sum == Matrix(2, 2, {11.0, 22.0, 33.0, 44.0}));
    CHECK(a == Matrix(2, 2, {1.0, 2.0, 3.0, 4.0}));

    CHECK(a + b + a == Matrix(2, 2, {12.0, 24.0, 36.0, 48.0}));
}

void test_multiply() {
    testharness::section("multiplication");

    Matrix a(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    Matrix b(3, 2, {7.0, 8.0, 9.0, 10.0, 11.0, 12.0});

    // [1 2 3] [ 7  8]   [ 58  64]
    // [4 5 6] [ 9 10] = [139 154]
    //         [11 12]
    Matrix p = a * b;
    CHECK(p.rows() == 2 && p.cols() == 2);
    CHECK_NEAR(p(0, 0), 58.0, 1e-12);
    CHECK_NEAR(p(0, 1), 64.0, 1e-12);
    CHECK_NEAR(p(1, 0), 139.0, 1e-12);
    CHECK_NEAR(p(1, 1), 154.0, 1e-12);

    // Not commutative: b*a is 3x3, a*b is 2x2.
    CHECK((b * a).rows() == 3);
    CHECK((b * a).cols() == 3);

    Matrix i3(3, 3, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
    CHECK(a * i3 == a);

    CHECK(a * 2.0 == Matrix(2, 3, {2.0, 4.0, 6.0, 8.0, 10.0, 12.0}));
    CHECK(2.0 * a == a * 2.0);
}

void test_dimension_mismatch() {
    testharness::section("dimension mismatch throws");

    Matrix a(2, 3);
    Matrix b(3, 2);
    Matrix c(2, 3);

    CHECK_THROWS(DimensionMismatch, a + b);
    CHECK_THROWS(DimensionMismatch, b + a);
    CHECK_THROWS(DimensionMismatch, a - b);

    // a is 2x3 and c is 2x3: inner dimensions 3 and 2 do not meet.
    CHECK_THROWS(DimensionMismatch, a * c);

    // Must be catchable at the REPL boundary as the base type.
    CHECK_THROWS(MatrixError, a + b);
    CHECK_THROWS(std::exception, a + b);

    // operator+= validates before touching anything, so a failed += leaves the
    // target exactly as it was (strong exception guarantee).
    Matrix keep(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    const Matrix before = keep;
    CHECK_THROWS(DimensionMismatch, keep += b);
    CHECK(keep == before);
}

void test_transpose() {
    testharness::section("transpose");

    Matrix m(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    Matrix t = m.transpose();

    // Non-square: the shape must actually flip.
    CHECK(t.rows() == 3);
    CHECK(t.cols() == 2);
    CHECK(t == Matrix(3, 2, {1.0, 4.0, 2.0, 5.0, 3.0, 6.0}));

    CHECK(t.transpose() == m);  // (M')' == M

    Matrix b(3, 2, {7.0, 8.0, 9.0, 10.0, 11.0, 12.0});
    CHECK((m * b).transpose() == b.transpose() * m.transpose());  // (AB)' == B'A'

    // The result is a fresh buffer, not a view onto the source.
    t(0, 0) = 42.0;
    CHECK_NEAR(m(0, 0), 1.0, 1e-12);
}

void test_value_semantics() {
    testharness::section("value semantics");

    Matrix a(2, 2, {1.0, 2.0, 3.0, 4.0});

    Matrix b = a;  // copy ctor
    b(0, 0) = 99.0;
    CHECK_NEAR(a(0, 0), 1.0, 1e-12);  // deep copy, not aliased
    CHECK_NEAR(b(0, 0), 99.0, 1e-12);

    Matrix c(1, 1);
    c = a;  // copy assignment must reshape, not just overwrite elements
    CHECK(c.rows() == 2 && c.cols() == 2 && c.size() == 4);
    CHECK(c == a);
}

void test_self_assignment() {
    testharness::section("self-assignment");

    Matrix a(2, 2, {1.0, 2.0, 3.0, 4.0});

    copy_assign(a, a);
    CHECK(a.rows() == 2 && a.cols() == 2);
    CHECK(a == Matrix(2, 2, {1.0, 2.0, 3.0, 4.0}));

    // Self-move-assignment leaves std::vector "valid but unspecified", so the
    // only thing guaranteed is that `a` remains destructible and assignable.
    // Assert that, and nothing about its contents.
    move_assign(a, a);
    a = Matrix(1, 1, {7.0});
    CHECK(a == Matrix(1, 1, {7.0}));
}

void test_moved_from_state() {
    testharness::section("moved-from state");

    Matrix src(2, 2, {1.0, 2.0, 3.0, 4.0});
    Matrix dst = std::move(src);

    // The part that actually matters: the destination holds the whole value.
    CHECK(dst == Matrix(2, 2, {1.0, 2.0, 3.0, 4.0}));

    // The sharp edge of Rule of Zero. vector's move ctor is required to be
    // O(1), which forces it to steal the buffer, so src's data is empty. But
    // rows_/cols_ are scalars -- "moving" a scalar copies it -- so src still
    // claims to be 2x2 while holding zero elements. The invariant
    // data_.size() == rows*cols is broken. This documents observed behaviour;
    // it is not a guarantee the class makes.
    CHECK(src.size() == 0);
    CHECK(src.rows() == 2);  // stale
    CHECK(src.cols() == 2);  // stale

    // Consequence: src(0,0) would index a zero-length buffer. Undefined
    // behaviour, so it is not called here. Assignment is the supported way back
    // to a consistent object.
    src = Matrix(3, 1, {1.0, 2.0, 3.0});
    CHECK(src.rows() == 3 && src.cols() == 1 && src.size() == 3);
    CHECK(src == Matrix(3, 1, {1.0, 2.0, 3.0}));
}

void test_1x1() {
    testharness::section("1x1 matrices");

    Matrix a(1, 1, {3.0});
    Matrix b(1, 1, {4.0});

    CHECK(a.rows() == 1 && a.cols() == 1 && a.size() == 1);
    CHECK(!a.empty());
    CHECK(a + b == Matrix(1, 1, {7.0}));
    CHECK(a * b == Matrix(1, 1, {12.0}));
    CHECK(a.transpose() == a);
    CHECK(a * 2.0 == Matrix(1, 1, {6.0}));
}

void test_empty() {
    testharness::section("empty matrices");

    Matrix d;  // default: 0x0
    CHECK(d.rows() == 0 && d.cols() == 0);
    CHECK(d.empty());
    CHECK(d + d == d);
    CHECK(d * d == d);
    CHECK(d.transpose() == d);

    // A 2x0 matrix holds no elements but is not the same object as 0x0 --
    // the dimensions still carry information.
    Matrix a(2, 0);
    CHECK(a.empty());
    CHECK(a.rows() == 2 && a.cols() == 0);
    CHECK(a.transpose().rows() == 0);
    CHECK(a.transpose().cols() == 2);
    CHECK(a != d);

    // 2x0 * 0x3 is legal: the inner dimension is 0, and the result is a 2x3
    // zero matrix -- an empty operand producing a non-empty result.
    Matrix b(0, 3);
    Matrix c = a * b;
    CHECK(c.rows() == 2 && c.cols() == 3);
    CHECK(!c.empty());
    CHECK(c == Matrix(2, 3));
}

void test_at_bounds() {
    testharness::section("at() bounds checking");

    Matrix m(2, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    CHECK_NEAR(m.at(1, 2), 6.0, 1e-12);

    m.at(0, 0) = 10.0;
    CHECK_NEAR(m.at(0, 0), 10.0, 1e-12);

    CHECK_THROWS(OutOfRange, m.at(2, 0));
    CHECK_THROWS(OutOfRange, m.at(0, 3));
    CHECK_THROWS(MatrixError, m.at(99, 99));

    const Matrix& cm = m;
    CHECK_NEAR(cm.at(1, 1), 5.0, 1e-12);
    CHECK_THROWS(OutOfRange, cm.at(2, 0));

    // On a 2x0 matrix every column index is out of range; at() must reject
    // rather than computing r*0 + c and indexing an empty buffer.
    Matrix e(2, 0);
    CHECK_THROWS(OutOfRange, e.at(0, 0));
    CHECK_THROWS(OutOfRange, e.at(1, 0));
}

} // namespace

int main() {
    test_construction();
    test_add_sub();
    test_multiply();
    test_dimension_mismatch();
    test_transpose();
    test_value_semantics();
    test_self_assignment();
    test_moved_from_state();
    test_1x1();
    test_empty();
    test_at_bounds();

    return testharness::summary();
}
