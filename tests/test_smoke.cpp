#include "test_harness.hpp"

#include <stdexcept>

int main() {
    CHECK(1 + 1 == 2);
    CHECK_NEAR(1.0, 1.0000001, 1e-6);
    CHECK_THROWS(([] { throw std::runtime_error("x"); })(), std::runtime_error);

    return testharness::failures() == 0 ? 0 : 1;
}
