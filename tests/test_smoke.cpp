#include "test_harness.hpp"

#include <stdexcept>

int main() {
    testharness::section("harness self-test");

    CHECK(1 + 1 == 2);
    CHECK_NEAR(1.0, 1.0000001, 1e-6);
    CHECK_THROWS(std::runtime_error, ([] { throw std::runtime_error("x"); })());

    return testharness::summary();
}
