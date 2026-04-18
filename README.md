# SecurityLibrary

`SecurityLibrary` is a C++ project for security-focused implementations.

The first implemented feature is the `Matrix` class. It provides the math foundation for later security modules and currently includes arithmetic, determinant/inverse logic, LU decomposition, rank, and linear-system solving through the runnable sample (`test_runner`).

## What you can do with it

- Build `Matrix` values from dimensions, `std::vector<std::vector<double>>`, or initializer lists
- Read/write elements and query shape (`get`, `set`, `getRows`, `getColumns`)
- Apply matrix transforms and checks (`transpose`, `isSquare`, `isDiagonal`, `isSymmetric`, `isZero`)
- Use arithmetic operators (`+`, `-`, `*`, scalar `*` and `/`, plus in-place variants)
- Run LU-based operations (`luDecomposition`, `determinantLU`, `rank`)
- Solve systems like `Ax = b` with `solveLinearSystem`

## What comes next

After the matrix foundation, the next features are planned as real security implementations (for example, cryptography/security-oriented modules built on top of the core math layer).

Public API: `include/Math/matrix.hpp`  
Implementation: `src/Math/matrix.cpp`

## Project layout

```text
security-library/
  CMakeLists.txt
  include/Math/matrix.hpp
  src/Math/matrix.cpp
  test/test.cpp
```

## Build and run

This project uses CMake.

- Configured default in `CMakeLists.txt`: **C++20** (`set(CMAKE_CXX_STANDARD 20)`)
- Lowest version suggested by current code features: **C++11**

If you want to build with C++11, change `set(CMAKE_CXX_STANDARD 20)` to `set(CMAKE_CXX_STANDARD 11)` in `CMakeLists.txt`.

### Option 1: Build from terminal with CMake

Requires `cmake` to be installed and available in `PATH`.

```powershell
Set-Location .
cmake -S . -B build
cmake --build build --config Debug
.\build\test_runner.exe
```

### Option 2: Run the existing CLion debug executable

If `cmake` is not available in your shell and `cmake-build-debug/test_runner.exe` already exists:

```powershell
Set-Location .\cmake-build-debug
.\test_runner.exe
```

## Example

`test/test.cpp` shows a simple `Ax = b` solve:

```cpp
#include <iostream>
#include "../include/Math/matrix.hpp"

int main() {
    Matrix A = {{3,2,1}, {1,0,2}, {4,1,3}};
    Matrix B = {{9}, {5}, {14}};

    Matrix x = Matrix::solveLinearSystem(A, B);
    std::cout << "Solution x:\n" << x;
}
```