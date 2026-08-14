#pragma once

#include <cmath>
#include <iostream>
#include <string>

namespace testharness {

struct Stats {
    int checks = 0;
    int failures = 0;
};

inline Stats& stats() {
    static Stats s;
    return s;
}

inline void report(bool ok, const std::string& expr, const char* file, int line) {
    ++stats().checks;
    if (!ok) {
        ++stats().failures;
        std::cerr << "  FAIL  " << file << ":" << line << ": " << expr << "\n";
    }
}

inline void section(const char* name) {
    std::cout << "-- " << name << "\n";
}

inline int summary() {
    const Stats& s = stats();
    std::cout << (s.failures == 0 ? "PASS  " : "FAIL  ") << (s.checks - s.failures)
              << "/" << s.checks << " checks passed\n";
    return s.failures == 0 ? 0 : 1;
}

} // namespace testharness

// Variadic so that commas inside the expression -- CHECK(f(a, b)) -- are not
// split into separate macro arguments by the preprocessor.
#define CHECK(...) \
    testharness::report(static_cast<bool>(__VA_ARGS__), #__VA_ARGS__, __FILE__, __LINE__)

// Not variadic: three fixed parameters, so an expression containing a comma
// must be parenthesised by the caller.
#define CHECK_NEAR(a, b, eps)                                    \
    testharness::report(std::fabs((a) - (b)) <= (eps),           \
                        #a " ~= " #b " (within " #eps ")",       \
                        __FILE__, __LINE__)

// Exception type comes first so the expression can be variadic (a variadic
// parameter must be last). A wrong exception type is reported as a failure
// rather than escaping and terminating the test binary.
#define CHECK_THROWS(exc_type, ...)                                   \
    do {                                                              \
        bool threw_ = false;                                          \
        try {                                                         \
            (void)(__VA_ARGS__);                                      \
        } catch (const exc_type&) {                                   \
            threw_ = true;                                            \
        } catch (...) {                                               \
            threw_ = false;                                           \
        }                                                             \
        testharness::report(threw_, #__VA_ARGS__ " throws " #exc_type, \
                            __FILE__, __LINE__);                      \
    } while (0)
