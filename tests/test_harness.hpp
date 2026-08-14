#pragma once

#include <cmath>
#include <iostream>
#include <string>

namespace testharness {

inline int& failures() {
    static int count = 0;
    return count;
}

inline void report(bool ok, const std::string& expr, const char* file, int line) {
    if (!ok) {
        ++failures();
        std::cerr << file << ":" << line << ": CHECK failed: " << expr << "\n";
    }
}

} // namespace testharness

#define CHECK(cond) \
    testharness::report((cond), #cond, __FILE__, __LINE__)

#define CHECK_NEAR(a, b, eps) \
    testharness::report(std::fabs((a) - (b)) <= (eps), #a " ~= " #b, __FILE__, __LINE__)

#define CHECK_THROWS(expr, exc_type)                                                   \
    do {                                                                               \
        bool threw = false;                                                            \
        try {                                                                          \
            (void)(expr);                                                              \
        } catch (const exc_type&) {                                                    \
            threw = true;                                                              \
        }                                                                              \
        testharness::report(threw, #expr " throws " #exc_type, __FILE__, __LINE__);    \
    } while (0)
