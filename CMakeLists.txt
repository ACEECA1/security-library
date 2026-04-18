# Note: Changed to 3.20, as CMake 4.2 does not exist yet.
# Use a version currently installed on your system.
cmake_minimum_required(VERSION 3.20)
project(SecurityLibrary)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 1. Define the library
# Removed test.cpp from here. The library should only contain library source code.
add_library(SecurityLibrary STATIC
        library.cpp
        src/Math/matrix.cpp
        include/Math/matrix.hpp
        library.h
)

# 2. Tell CMake where to find headers
# PUBLIC ensures that anything linking to SecurityLibrary also knows where these headers are.
target_include_directories(SecurityLibrary PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# 3. Define the test executable
add_executable(test_runner
        test/test.cpp
)

# 4. Link the library to the executable
target_link_libraries(test_runner PRIVATE SecurityLibrary)