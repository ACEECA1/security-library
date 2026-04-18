#include "../../include/Math/matrix.hpp"
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
namespace {
constexpr double EPS = 1e-12;

bool nearZero(double x) {
    return std::abs(x) <= EPS;
}
}

Matrix::Matrix(const Matrix& original) {
    matrix = original.matrix;
    rows = original.rows;
    columns = original.columns;
}

Matrix::Matrix(size_t rows, size_t columns) : rows(rows), columns(columns) {
    matrix.assign(rows, std::vector<double>(columns, 0.0));
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) {
    rows = list.size();
    columns = (rows > 0) ? list.begin()->size() : 0;
    matrix.resize(rows, std::vector<double>(columns, 0.0));

    size_t i = 0;
    for (const auto& row_list : list) {
        if (row_list.size() != columns) {
            throw std::invalid_argument("All initializer-list rows must have the same length");
        }
        std::copy(row_list.begin(), row_list.end(), matrix[i].begin());
        i++;
    }
}

Matrix::Matrix(size_t rows, size_t columns, double value) : rows(rows), columns(columns) {
    matrix.assign(rows, std::vector<double>(columns, value));
}

Matrix::Matrix(const std::vector<std::vector<double>>& input) {
    rows = input.size();
    if (rows == 0) {
        columns = 0;
        return;
    }

    columns = input[0].size();
    for (const auto& row : input) {
        if (row.size() != columns) {
            throw std::invalid_argument("All rows must have the same length");
        }
    }

    matrix = input;
}

double Matrix::get(size_t i, size_t j) const {
    if (i >= rows || j >= columns) {
        throw std::out_of_range("Index out of range");
    }
    return matrix[i][j];
}

void Matrix::set(size_t i, size_t j, double value) {
    if (i >= rows || j >= columns) {
        throw std::out_of_range("Index out of range");
    }
    matrix[i][j] = value;
}

size_t Matrix::getRows() const {
    return rows;
}

size_t Matrix::getColumns() const {
    return columns;
}

void Matrix::setRow(const std::vector<double>& row, size_t i) {
    if (i >= rows) {
        throw std::out_of_range("Index out of range");
    }
    if (row.size() != columns) {
        throw std::length_error("Row has incorrect length");
    }
    matrix[i] = row;
}

void Matrix::setColumn(const std::vector<double>& column, size_t j) {
    if (j >= columns) {
        throw std::out_of_range("Index out of range");
    }
    if (column.size() != rows) {
        throw std::length_error("Column has incorrect length");
    }
    for (size_t i = 0; i < rows; i++) {
        matrix[i][j] = column[i];
    }
}

void Matrix::swapRow(size_t i, size_t j) {
    if (i >= rows || j >= rows) {
        throw std::out_of_range("Index out of range");
    }
    std::swap(matrix[i], matrix[j]);
}

void Matrix::swapColumn(size_t i, size_t j) {
    if (i >= columns || j >= columns) {
        throw std::out_of_range("Index out of range");
    }
    for (size_t k = 0; k < rows; k++) {
        std::swap(matrix[k][i], matrix[k][j]);
    }
}

void Matrix::transpose() {
    Matrix transposed(columns, rows);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            transposed.matrix[j][i] = matrix[i][j];
        }
    }
    matrix = std::move(transposed.matrix);
    std::swap(rows, columns);
}

bool Matrix::isSquare() const {
    return rows == columns;
}

bool Matrix::isDiagonal() const {
    if (rows != columns) return false;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < i; j++) {
            if (!nearZero(matrix[i][j])) {
                return false;
            }
        }
        for (size_t j = i + 1; j < columns; j++) {
            if (!nearZero(matrix[i][j])) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::isSymmetric() const {
    if (!isSquare()) return false;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = i + 1; j < columns; j++) {
            if (std::abs(matrix[i][j] - matrix[j][i]) > EPS) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::isZero() const {
    for (size_t i = 0; i < rows; i++) {
        const auto& row = matrix[i];
        for (size_t j = 0; j < columns; j++) {
            if (!nearZero(row[j])) {
                return false;
            }
        }
    }
    return true;
}

double Matrix::determinant() const {
    return determinant(*this);
}

Matrix Matrix::getCofactor() const {
    if (!isSquare()) {
        throw std::logic_error("Cofactor requires square matrix");
    }
    if (rows == 0) {
        return Matrix();
    }
    if (rows == 1) {
        return Matrix{{1.0}};
    }

    Matrix cofactor(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            cofactor.matrix[i][j] = determinant(chop(*this, i, j)) * ((i + j) % 2 == 0 ? 1.0 : -1.0);
        }
    }
    return cofactor;
}

Matrix Matrix::getInverse() const {
    double det = determinant(*this);
    if (nearZero(det)) {
        throw std::logic_error("The matrix is not invertible");
    }
    Matrix adjoint = getCofactor();
    adjoint.transpose();
    return adjoint * (1.0 / det);
}

void Matrix::invert() {
    double det = determinant(*this);
    if (nearZero(det)) {
        throw std::logic_error("The matrix is not invertible");
    }
    Matrix adjoint = getCofactor();
    adjoint.transpose();
    *this = adjoint * (1.0 / det);
}

double Matrix::trace() const {
    if (!isSquare()) {
        throw std::logic_error("Trace requires square matrix");
    }
    double sum = 0.0;
    for (size_t i = 0; i < rows; i++) {
        sum += matrix[i][i];
    }
    return sum;
}

Matrix Matrix::copy() const {
    return Matrix(*this);
}

Matrix Matrix::identity(size_t n) {
    Matrix id(n, n, 0.0);
    for (size_t i = 0; i < n; i++) {
        id[i][i] = 1.0;
    }
    return id;
}

double Matrix::determinantOld(const Matrix& m) {
    if (!m.isSquare()) {
        throw std::out_of_range("The matrix has to be a square matrix");
    }
    if (m.rows == 0) {
        return 1.0;
    }
    if (m.rows == 1) {
        return m[0][0];
    }
    if (m.rows == 2) {
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }

    double det = 0.0;
    int mult = 1;
    for (size_t j = 0; j < m.columns; j++) {
        det += mult * m.get(0, j) * determinantOld(chop(m, 0, j));
        mult *= -1;
    }
    return det;
}

double Matrix::determinant(const Matrix& m) {
    return m.determinantLU();
}

Matrix Matrix::chop(const Matrix& m, size_t row, size_t col) {
    if (m.rows <= 1 || m.columns <= 1) {
        throw std::logic_error("Cannot chop a matrix smaller than 2x2");
    }
    if (row >= m.rows || col >= m.columns) {
        throw std::out_of_range("Index out of range");
    }

    Matrix result(m.rows - 1, m.columns - 1);
    size_t a = 0;

    for (size_t i = 0; i < m.rows; i++) {
        if (i == row) continue;

        size_t b = 0;
        for (size_t j = 0; j < m.columns; j++) {
            if (j == col) continue;
            result.matrix[a][b] = m.matrix[i][j];
            b++;
        }
        a++;
    }

    return result;
}

Matrix Matrix::solveLinearSystem(const Matrix& A, const Matrix& b) {
    if (!A.isSquare()) {
        throw std::logic_error("Matrix A must be square");
    }
    if (A.getRows() != b.getRows()) {
        throw std::logic_error("Dimension mismatch");
    }
    if (b.getColumns() != 1) {
        throw std::logic_error("Right-hand side must be a column vector");
    }

    Matrix L, U;
    size_t swapCount = 0;
    A.luDecomposition(L, U, swapCount);

    size_t n = A.getRows();
    Matrix y(n, 1);
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        const auto& Li = L[i];
        for (size_t k = 0; k < i; k++) {
            sum += Li[k] * y[k][0];
        }
        y[i][0] = b[i][0] - sum;
    }

    Matrix x(n, 1);
    for (int i = static_cast<int>(n) - 1; i >= 0; i--) {
        double sum = 0.0;
        const auto& Ui = U[static_cast<size_t>(i)];
        for (size_t k = static_cast<size_t>(i) + 1; k < n; k++) {
            sum += Ui[k] * x[k][0];
        }
        const double pivot = Ui[static_cast<size_t>(i)];
        if (nearZero(pivot)) {
            throw std::logic_error("The matrix is singular or nearly singular");
        }
        x[static_cast<size_t>(i)][0] = (y[static_cast<size_t>(i)][0] - sum) / pivot;
    }

    return x;
}

double Matrix::determinantLU() const {
    if (!isSquare()) {
        throw std::logic_error("Matrix must be square");
    }

    Matrix L, U;
    size_t swapCount = 0;
    this->luDecomposition(L, U, swapCount);

    double det = (swapCount % 2 == 0) ? 1.0 : -1.0;
    for (size_t i = 0; i < rows; i++) {
        det *= U[i][i];
    }
    return det;
}

void Matrix::luDecomposition(Matrix& L, Matrix& U) const {
    size_t swapCount = 0;
    luDecomposition(L, U, swapCount);
}

void Matrix::luDecomposition(Matrix& L, Matrix& U, size_t& swapCount) const {
    if (!isSquare()) {
        throw std::logic_error("Matrix must be square");
    }

    const size_t n = rows;
    L = Matrix::identity(n);
    U = *this;
    swapCount = 0;

    for (size_t i = 0; i < n; i++) {
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; k++) {
            if (std::abs(U[k][i]) > std::abs(U[maxRow][i])) {
                maxRow = k;
            }
        }

        if (nearZero(U[maxRow][i])) {
            throw std::logic_error("Matrix is singular");
        }

        if (maxRow != i) {
            U.swapRow(i, maxRow);
            for (size_t col = 0; col < i; col++) {
                std::swap(L[i][col], L[maxRow][col]);
            }
            swapCount++;
        }

        const double pivot = U[i][i];
        for (size_t k = i + 1; k < n; k++) {
            const double factor = U[k][i] / pivot;
            L[k][i] = factor;
            for (size_t j = i; j < n; j++) {
                U[k][j] -= factor * U[i][j];
            }
        }
    }
}

size_t Matrix::rank() const {
    Matrix L, U;
    size_t swapCount = 0;
    luDecomposition(L, U, swapCount);

    size_t r = 0;
    for (size_t i = 0; i < U.getRows(); i++) {
        bool nonZeroRow = false;
        const auto& row = U[i];
        for (size_t j = 0; j < U.getColumns(); j++) {
            if (!nearZero(row[j])) {
                nonZeroRow = true;
                break;
            }
        }
        if (nonZeroRow) {
            r++;
        }
    }
    return r;
}

Matrix Matrix::operator+(const Matrix& rhs) const {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    Matrix result(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        const auto& lhsRow = matrix[i];
        const auto& rhsRow = rhs.matrix[i];
        auto& outRow = result.matrix[i];
        for (size_t j = 0; j < columns; j++) {
            outRow[j] = lhsRow[j] + rhsRow[j];
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& rhs) const {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    Matrix result(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        const auto& lhsRow = matrix[i];
        const auto& rhsRow = rhs.matrix[i];
        auto& outRow = result.matrix[i];
        for (size_t j = 0; j < columns; j++) {
            outRow[j] = lhsRow[j] - rhsRow[j];
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& rhs) const {
    if (columns != rhs.rows) {
        throw std::length_error("The length of the left hand matrix's cols should equal the right hand's rows");
    }
    Matrix result(rows, rhs.columns);
    for (size_t i = 0; i < rows; i++) {
        const auto& lhsRow = matrix[i];
        auto& outRow = result.matrix[i];
        for (size_t k = 0; k < columns; k++) {
            const double lhs = lhsRow[k];
            const auto& rhsRow = rhs.matrix[k];
            for (size_t j = 0; j < rhs.columns; j++) {
                outRow[j] += lhs * rhsRow[j];
            }
        }
    }
    return result;
}

Matrix Matrix::operator/(const Matrix& rhs) const {
    return (*this) * rhs.getInverse();
}

Matrix Matrix::operator/(const double& scalar) const {
    if (nearZero(scalar)) {
        throw std::domain_error("Division by zero");
    }
    return (*this) * (1.0 / scalar);
}

Matrix Matrix::operator*(const double& scalar) const {
    Matrix result(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        const auto& lhsRow = matrix[i];
        auto& outRow = result.matrix[i];
        for (size_t j = 0; j < columns; j++) {
            outRow[j] = scalar * lhsRow[j];
        }
    }
    return result;
}

void Matrix::operator+=(const Matrix& rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    for (size_t i = 0; i < rows; i++) {
        auto& lhsRow = matrix[i];
        const auto& rhsRow = rhs.matrix[i];
        for (size_t j = 0; j < columns; j++) {
            lhsRow[j] += rhsRow[j];
        }
    }
}

void Matrix::operator-=(const Matrix& rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    for (size_t i = 0; i < rows; i++) {
        auto& lhsRow = matrix[i];
        const auto& rhsRow = rhs.matrix[i];
        for (size_t j = 0; j < columns; j++) {
            lhsRow[j] -= rhsRow[j];
        }
    }
}

void Matrix::operator*=(const Matrix& rhs) {
    if (columns != rhs.rows) {
        throw std::length_error("The length of the left hand matrix's cols should equal the right hand's rows");
    }
    Matrix result(rows, rhs.columns);
    for (size_t i = 0; i < rows; i++) {
        const auto& lhsRow = matrix[i];
        auto& outRow = result.matrix[i];
        for (size_t k = 0; k < columns; k++) {
            const double lhs = lhsRow[k];
            const auto& rhsRow = rhs.matrix[k];
            for (size_t j = 0; j < rhs.columns; j++) {
                outRow[j] += lhs * rhsRow[j];
            }
        }
    }
    *this = std::move(result);
}

void Matrix::operator/=(const double& scalar) {
    if (nearZero(scalar)) {
        throw std::domain_error("Division by zero");
    }
    *this = (*this) * (1.0 / scalar);
}

void Matrix::operator*=(const double& scalar) {
    for (size_t i = 0; i < rows; i++) {
        auto& lhsRow = matrix[i];
        for (size_t j = 0; j < columns; j++) {
            lhsRow[j] *= scalar;
        }
    }
}

bool Matrix::operator==(const Matrix& rhs) const {
    if (rows != rhs.rows || columns != rhs.columns) {
        return false;
    }
    for (size_t i = 0; i < rows; i++) {
        const auto& lhsRow = matrix[i];
        const auto& rhsRow = rhs.matrix[i];
        for (size_t j = 0; j < columns; j++) {
            if (std::abs(lhsRow[j] - rhsRow[j]) > EPS) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::operator!=(const Matrix& rhs) const {
    return !(*this == rhs);
}

std::vector<double>& Matrix::operator[](size_t i) {
    if (i >= rows) {
        throw std::out_of_range("Index out of bounds");
    }
    return matrix[i];
}

const std::vector<double>& Matrix::operator[](size_t i) const {
    if (i >= rows) {
        throw std::out_of_range("Index out of bounds");
    }
    return matrix[i];
}

std::string Matrix::toString() const {
    std::string result;
    for (size_t i = 0; i < rows; i++) {
        result += " | ";
        const auto& row = matrix[i];
        for (size_t j = 0; j < columns; j++) {
            result += std::to_string(row[j]);
            result += " | ";
        }
        result += "\n";
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    os << "Matrix (" << m.rows << "x" << m.columns << "):\n";
    for (size_t i = 0; i < m.rows; i++) {
        os << " | ";
        const auto& row = m.matrix[i];
        for (size_t j = 0; j < m.columns; j++) {
            os << std::setw(5) << row[j] << " | ";
        }
        os << "\n";
    }
    return os;
}

Matrix operator*(const double& scalar, const Matrix& m) {
    return m * scalar;
}