cmake_minimum_required(VERSION 3.16)
project(mini_matlab CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)

add_library(mini_matlab_core INTERFACE)
target_include_directories(mini_matlab_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

enable_testing()

add_executable(mini_matlab_tests tests/test_smoke.cpp)
target_link_libraries(mini_matlab_tests PRIVATE mini_matlab_core)

add_test(NAME mini_matlab_tests COMMAND mini_matlab_tests)
