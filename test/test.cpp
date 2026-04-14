#include <iostream>
#include "../include/Math/matrix.hpp"

int main() {
    Matrix<int> m({{2, 0, 0, 1}, {0, 1, 3, 0}, {4, 0, 1, 0}, {0, -1, 0, 2}});

    std::cout << m.toString() << " det = " << m.determinant() << std::endl;
}
