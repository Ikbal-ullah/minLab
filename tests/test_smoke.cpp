#include "test_harness.hpp"
#include "mini_matlab/matrix.hpp"

#include <stdexcept>

int main() {
    testharness::section("harness self-test");

    CHECK(1 + 1 == 2);
    CHECK_NEAR(1.0, 1.0000001, 1e-6);
    CHECK_THROWS(std::runtime_error, ([] { throw std::runtime_error("x"); })());

    testharness::section("matrix scalar division");

    mini_matlab::Matrix m(2, 2, {2.0, 4.0, 6.0, 8.0});

    // Test binary operator/
    mini_matlab::Matrix div_result = m / 2.0;
    CHECK(div_result == mini_matlab::Matrix(2, 2, {1.0, 2.0, 3.0, 4.0}));

    // Test compound operator/=
    m /= 2.0;
    CHECK(m == mini_matlab::Matrix(2, 2, {1.0, 2.0, 3.0, 4.0}));

    // Test division by zero exception
    CHECK_THROWS(mini_matlab::DivisionByZero, m / 0.0);

    return testharness::summary();
}