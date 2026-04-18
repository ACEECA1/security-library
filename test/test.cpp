#include <iostream>
#include "../include/Math/matrix.hpp"

int main() {
    Matrix A = {{3,2,1}, {1,0,2}, {4,1,3}};
    Matrix B = {{9}, {5}, {14}};
    try {
        Matrix x = Matrix::solveLinearSystem(A, B);
        std::cout << "Solution x:\n" << x;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}
