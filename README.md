# SecurityLibrary

`SecurityLibrary` is a small C++20 matrix library.

Right now, the project is centered around one type: `Matrix`. It supports common linear algebra operations (arithmetic, determinant, inverse, LU decomposition, rank) and includes a runnable sample (`test_runner`) that solves a linear system.

## What you can do with it

- Build `Matrix` values from dimensions, `std::vector<std::vector<double>>`, or initializer lists
- Read/write elements and query shape (`get`, `set`, `getRows`, `getColumns`)
- Apply matrix transforms and checks (`transpose`, `isSquare`, `isDiagonal`, `isSymmetric`, `isZero`)
- Use arithmetic operators (`+`, `-`, `*`, scalar `*` and `/`, plus in-place variants)
- Run LU-based operations (`luDecomposition`, `determinantLU`, `rank`)
- Solve systems like `Ax = b` with `solveLinearSystem`

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

This project uses CMake and targets C++20.

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