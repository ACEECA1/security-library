cmake_minimum_required(VERSION 4.2)
project(SecurityLibrary)

set(CMAKE_CXX_STANDARD 20)

add_library(SecurityLibrary STATIC library.cpp
        src/Math/matrix.tpp
        include/Math/matrix.hpp
        test/test.cpp)
add_executable(test_runner
        test/test.cpp
)

target_link_libraries(test_runner PRIVATE SecurityLibrary)