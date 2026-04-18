cmake_minimum_required(VERSION 3.20)
project(SecurityLibrary)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(SecurityLibrary STATIC
        src/Math/matrix.cpp
        include/Math/matrix.hpp
)

target_include_directories(SecurityLibrary PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_runner
        test/test.cpp
)


target_link_libraries(test_runner PRIVATE SecurityLibrary)
