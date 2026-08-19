cmake_minimum_required(VERSION 3.16)
project(mini_matlab CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Must precede every add_library/add_executable below: a target snapshots the
# COMPILE_OPTIONS directory property at the moment it is created, so targets
# declared above this line would silently miss these flags.
add_compile_options(-Wall -Wextra -Wpedantic -Werror)

option(MINI_MATLAB_TRACE_SPECIAL
       "Print a line on every Matrix copy/move/destroy" OFF)

add_library(mini_matlab_core STATIC
    src/matrix.cpp
    src/solve.cpp
)

# PUBLIC: matrix.cpp needs include/ to find its own header, and so does every
# consumer that links against this target. PRIVATE would build the library but
# break the tests; INTERFACE would break the library itself.
target_include_directories(mini_matlab_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# PUBLIC is load-bearing, not stylistic. This macro changes which special
# members Matrix declares. If the library saw the declarations and a test did
# not, the two would disagree on the class definition -- an ODR violation and
# undefined behaviour. PUBLIC forces every translation unit to agree.
if(MINI_MATLAB_TRACE_SPECIAL)
    target_compile_definitions(mini_matlab_core PUBLIC MINI_MATLAB_TRACE_SPECIAL)
endif()

enable_testing()

add_executable(smoke_tests tests/test_smoke.cpp)
target_link_libraries(smoke_tests PRIVATE mini_matlab_core)
add_test(NAME smoke COMMAND smoke_tests)

add_executable(matrix_tests tests/test_matrix.cpp)
target_link_libraries(matrix_tests PRIVATE mini_matlab_core)
add_test(NAME matrix COMMAND matrix_tests)

add_executable(solve_tests tests/test_solve.cpp)
target_link_libraries(solve_tests PRIVATE mini_matlab_core)
add_test(NAME solve COMMAND solve_tests)
