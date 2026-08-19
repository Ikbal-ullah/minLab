# NOTES.md

Interview revision notes. Written by me after each phase, in my own words: what
the design decision was, what I rejected, and why.

# Mini-MATLAB Architecture & Decision Log

## Phase 0 — Modern CMake & Build Architecture
* **Target-Centric Design:** Configured build using scoped targets (`add_library`, `add_executable`) rather than mutating global `CMAKE_CXX_FLAGS`.
* **Scope Rules:** `PRIVATE` for internal flags, `INTERFACE` for header-only targets, and `PUBLIC` for target + consumer inheritance. Attempting `PUBLIC`/`PRIVATE` on header-only targets causes a CMake configure-time error.
* **Flag Inheritance:** `add_compile_options()` appends to directory properties. Targets snapshot directory properties at creation time (targets created prior miss subsequent directory options).

## Phase 1 — Matrix Storage & Special Members
* **Flat 1D Buffer over Nested Vectors:** Used `std::vector<double>` indexed via `r * cols_ + c`. Avoids pointer-chasing double indirection, guarantees continuous L1 cache line prefetching (64-byte fetches load 8 doubles), and prevents non-rectangular ("ragged") matrix states.
* **Rule of Zero:** Outsourced memory management entirely to `std::vector`. Writing custom destructors suppresses implicit move constructor generation, silently downgrading $O(1)$ move operations to $O(n)$ deep copies across the codebase.
* **Free Function Operators:** Binary operators like `operator+` and `operator*` are free functions taking left-hand parameters by value (`Matrix lhs`). Allows compiler move optimization on temporaries in expressions like `(A + B) + C`.
* **Unchecked `operator()` vs Checked `at()`:** `operator()` is inline/noexcept for raw inner-loop speed; `at()` provides runtime bounds validation throwing `OutOfRange`.
* **Exception Hierarchy:** All engine exceptions derive from `MatrixError` so the REPL boundary catches all runtime errors via `catch (const MatrixError&)`. Diverged from IEEE 754 division by throwing `DivisionByZero` to catch invalid operations early.
## Phase 2 — Gaussian elimination + partial pivoting

## Phase 3 — lexer + recursive-descent parser

## Phase 4 — AST, evaluator, REPL

## Phase 5 — hardening, benchmarks, README
