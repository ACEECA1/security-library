#include "../../include/Math/matrix.hpp"
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>

Matrix::Matrix(const Matrix& original) {
    matrix = original.matrix;
    rows = original.rows;
    columns = original.columns;
}

Matrix::Matrix(size_t rows, size_t columns) : rows(rows), columns(columns) {
    matrix.clear();
    for (size_t i = 0; i < rows; i++) {
        matrix.push_back(std::vector<double>(columns, 0.0));
    }
}
Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) {
    rows = list.size();
    columns = (rows > 0) ? list.begin()->size() : 0;
    matrix.resize(rows, std::vector<double>(columns));

    size_t i = 0;
    for (auto& row_list : list) {
        std::copy(row_list.begin(), row_list.end(), matrix[i].begin());
        i++;
    }
}

Matrix::Matrix(size_t rows, size_t columns, double value) : rows(rows), columns(columns) {
    Matrix resized(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            resized[i][j] = value;
        }
    }
    matrix = resized.matrix;
}

Matrix::Matrix(const std::vector<std::vector<double>>& matrix) {
    rows = matrix.size();
    if (rows == 0) return;
    columns = matrix[0].size();
    std::vector<std::vector<double>> copy = matrix;
    for (size_t i = 0; i < rows; i++) {
        if (columns < copy[i].size()) {
            columns = copy[i].size();
        }
    }
    for (size_t i = 0; i < rows; i++) {
        copy[i].resize(columns, 0.0);
        this->matrix.push_back(copy[i]);
    }
}

double Matrix::get(size_t i, size_t j) const {
    if (i >= rows || j >= columns) {
        throw std::out_of_range("The Indexes are out of range");
    }
    return matrix[i][j];
}

void Matrix::set(size_t i, size_t j, double value) {
    if (i >= rows || j >= columns) {
        throw std::out_of_range("The Indexes are out of range");
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
        throw std::out_of_range("The Index is out of range");
    }
    if (row.size() != columns) {
        throw std::length_error("Not the same length");
    }
    matrix[i] = row;
}

void Matrix::setColumn(const std::vector<double>& column, size_t j) {
    if (j >= columns) {
        throw std::out_of_range("The Index is out of range");
    }
    if (column.size() != rows) {
        throw std::length_error("Not the same length");
    }
    for (size_t i = 0; i < rows; i++) {
        matrix[i][j] = column[i];
    }
}

void Matrix::swapRow(size_t i, size_t j) {
    std::swap(matrix[i], matrix[j]);
}

void Matrix::swapColumn(size_t i, size_t j) {
    if (i >= columns || j >= columns) {
        throw std::out_of_range("The Indexes are out of range");
    }
    for (size_t k = 0; k < rows; k++) {
        std::swap(matrix[k][j], matrix[k][i]);
    }
}

void Matrix::transpose() {
    Matrix transposed(columns, rows);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            transposed.matrix[j][i] = matrix[i][j];
        }
    }
    matrix = transposed.matrix;
    std::swap(rows, columns);
}

bool Matrix::isSquare() const {
    return rows == columns;
}

bool Matrix::isDiagonal() {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            if (i != j && matrix[i][j] != 0.0) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::isSymmetric() {
    if (!isSquare()) return false;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::isZero() {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            if (matrix[i][j] != 0.0) {
                return false;
            }
        }
    }
    return true;
}

double Matrix::determinant() {
    return determinant(*this);
}

Matrix Matrix::getCofactor() const{
    if (!isSquare()) {
        throw std::logic_error("Cofactor requires square matrix");
    }
    Matrix cofactor(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            cofactor.matrix[i][j] = determinant(chop(*this, i, j)) * ((i + j) % 2 == 0 ? 1.0 : -1.0);
        }
    }
    return cofactor;
}

Matrix Matrix::getInverse() const{
    if (determinant(*this) == 0.0) {
        throw std::logic_error("The matrix is not invertible");
    }
    Matrix inverse(rows, columns);
    Matrix adjoint = getCofactor();
    adjoint.transpose();
    inverse = adjoint * (1.0 / determinant(*this));
    return inverse;
}

void Matrix::invert() {
    if (determinant(*this) == 0.0) {
        throw std::logic_error("The matrix is not invertible");
    }
    Matrix inverse(rows, columns);
    Matrix adjoint = getCofactor();
    adjoint.transpose();
    inverse = adjoint * (1.0 / determinant(*this));
    matrix = inverse.matrix;
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

Matrix Matrix::copy() {
    Matrix copied(*this);
    return copied;
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
    if (m.rows == 2) {
        double a = m[0][0];
        double b = m[0][1];
        double c = m[1][0];
        double d = m[1][1];
        return (a * d - b * c);
    }
    double det = 0.0;
    int mult = 1;
    for (size_t j = 0; j < m.columns; j++) {
        det += mult * m.get(0, j) * determinantOld(chop(m, 0, j));
        mult *= -1;
    }
    return det;
}

double Matrix::determinant(const Matrix& m)
{
	return m.determinantLU();
}

size_t Matrix::rank() {
    return 0; // Keeping original implementation stub
}

Matrix Matrix::chop(const Matrix& m, size_t row, size_t col) {
    if (m.rows <= 1 || m.columns <= 1) {
        throw std::logic_error("Cannot chop a matrix smaller than 2x2");
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
    if (A.getRows() != b.getRows()) throw std::logic_error("Dimension mismatch");

    Matrix L, U;
    A.luDecomposition(L, U);
    size_t n = A.getRows();

    Matrix y(n, 1);
    for (size_t i = 0; i < n; i++) {
        double sum = 0;
        for (size_t k = 0; k < i; k++) sum += L[i][k] * y[k][0];
        y[i][0] = b[i][0] - sum;
    }

    Matrix x(n, 1);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (size_t k = i + 1; k < n; k++) sum += U[i][k] * x[k][0];
        x[i][0] = (y[i][0] - sum) / U[i][i];
    }

    return x;
}
double Matrix::determinantLU() const {
    Matrix L, U;
    try {
        this->luDecomposition(L, U);
    }
    catch (...) { return 0.0; }

    double det = 1.0;
    for (size_t i = 0; i < rows; i++) {
        det *= U[i][i];
    }
    return det;
}


void Matrix::luDecomposition(Matrix& L, Matrix& U) const {
    if (!isSquare()) throw std::logic_error("Matrix must be square");

    size_t n = rows;
    L = Matrix::identity(n);
    U = *this;

    for (size_t i = 0; i < n; i++) {
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; k++) {
            if (std::abs(U[k][i]) > std::abs(U[maxRow][i])) {
                maxRow = k;
            }
        }

        if (maxRow != i) {
            U.swapRow(i, maxRow);
        }

        for (size_t k = i + 1; k < n; k++) {
            if (std::abs(U[i][i]) < 1e-12) continue;

            double factor = U[k][i] / U[i][i];
            L[k][i] = factor;
            for (size_t j = i; j < n; j++) {
                U[k][j] -= factor * U[i][j];
            }
        }
    }
}

Matrix Matrix::operator+(const Matrix& rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    Matrix result(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            result.matrix[i][j] = matrix[i][j] + rhs.matrix[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    Matrix result(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            result.matrix[i][j] = matrix[i][j] - rhs.matrix[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& rhs) {
    if (columns != rhs.rows) {
        throw std::length_error("The length of the left hand matrix's cols should equal the right hand's rows");
    }
    Matrix result(rows, rhs.columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < rhs.columns; j++) {
            for (size_t k = 0; k < rhs.rows; k++) {
                result.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
            }
        }
    }
    return result;
}

Matrix Matrix::operator/(const Matrix& rhs) {
	Matrix inverse = rhs.getInverse();
	return (*this) * inverse;
}

Matrix Matrix::operator/(const double& scalar)
{
    return (*this) * (1.0 / scalar);
}

Matrix Matrix::operator*(const double& scalar) {
    Matrix result(rows, columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            result.matrix[i][j] = scalar * matrix[i][j];
        }
    }
    return result;
}

void Matrix::operator+=(const Matrix& rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            matrix[i][j] += rhs.matrix[i][j];
        }
    }
}

void Matrix::operator-=(const Matrix& rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            matrix[i][j] -= rhs.matrix[i][j];
        }
    }
}

void Matrix::operator*=(const Matrix& rhs) {
    if (columns != rhs.rows) {
        throw std::length_error("The length of the left hand matrix's cols should equal the right hand's rows");
    }
    Matrix result(rows, rhs.columns);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < rhs.columns; j++) {
            for (size_t k = 0; k < rhs.rows; k++) {
                result.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
            }
        }
    }
    matrix = result.matrix;
    rows = result.rows;
    columns = result.columns;
}

void Matrix::operator/=(const double& scalar)
{
	*this = (*this) * (1.0 / scalar);
}

void Matrix::operator*=(const double& scalar) {
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            matrix[i][j] *= scalar;
        }
    }
}

bool Matrix::operator==(const Matrix& rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        return false;
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            if (matrix[i][j] != rhs.matrix[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::operator!=(const Matrix& rhs) {
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
        for (size_t j = 0; j < columns; j++) {
            result += std::to_string(matrix[i][j]);
            result += " | ";
        }
        result += "\n";
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
	os << "Matrix (" << m.rows << "x" << m.columns << "):\n";
    for (size_t i = 0; i < m.rows; i++) {
        os << " | ";
        for (size_t j = 0; j < m.columns; j++) {
            os << m.matrix[i][j];
            os << " | ";
        }
        os << "\n";
    }
    return os;
}

Matrix operator*(const double& scalar, Matrix& m)
{
    return m * scalar;
}

Matrix operator/(const double& scalar, Matrix& m)
{
    return m.getInverse() * scalar;
}

