//
// Created by ASUS on 13/04/2026.
//
#include "../../include/math/matrix.hpp"
#include  <string>
#include <stdexcept>
#include <algorithm>


template<typename T>
Matrix<T>::Matrix(const Matrix &original) {
    matrix = original.matrix;
    rows = original.rows;
    columns = original.columns;
}

template<typename T>
Matrix<T>::Matrix(size_t rows, size_t columns):rows(rows), columns(columns) {
    matrix.clear();
    for (size_t i = 0 ; i < rows ; i++) {
        matrix.push_back(std::vector<T>(columns, 0));
    }
}

template<typename T>
Matrix<T>::Matrix(size_t rows, size_t columns, T value):rows(rows) , columns(columns) {
    Matrix<T> resized(rows, columns);
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            resized[i][j] = value;
        }
    }
    matrix = resized.matrix;
}

template<typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& matrix){

    rows = matrix.size();
    if (rows == 0) return;
    columns = matrix[0].size();
    std::vector<std::vector<T>> copy = matrix;
    for (size_t i = 0 ; i < rows ; i++) {
        if (columns < copy[i].size()) {
            columns = copy[i].size();
        }
    }
    for (size_t i = 0 ; i < rows ; i++) {
        copy[i].resize(columns, 0);
        this->matrix.push_back(copy[i]);
    }
}

template<typename T>
T Matrix<T>::get(size_t i, size_t j) const{
    if (i >= rows || j >= columns) {
        throw std::out_of_range("The Indexes are out of range");
    }
    return matrix[i][j];
}

template<typename T>
void Matrix<T>::set(size_t i, size_t j, T value) {
    if (i >= rows || j >= columns) {
        throw std::out_of_range("The Indexes are out of range");
    }
    matrix[i][j] = value;
}


template<typename T>
size_t Matrix<T>::getRows() const {
    return rows;
}

template<typename T>
size_t Matrix<T>::getColumns() const {
    return columns;
}

template<typename T>
void Matrix<T>::setRow(const std::vector<T> &row, size_t i) {
    if (i>=rows) {
        throw std::out_of_range("The Index is out of range");
    }
    if (row.size() != columns) {
        throw std::length_error("Not the same length");
    }
    matrix[i] = row;
}


template<typename T>
void Matrix<T>::setColumn(const std::vector<T> &column, size_t j) {
    if (j>=columns) {
        throw std::out_of_range("The Index is out of range");
    }
    if (column.size() != rows) {
        throw std::length_error("Not the same length");
    }
    for (size_t i = 0 ; i < rows ; i++) {
        matrix[i][j] = column[i];
    }
}

template<typename T>
void Matrix<T>::swapRow(size_t i, size_t j) {
    std::swap(matrix[i] , matrix[j]);
}


template<typename T>
void Matrix<T>::swapColumn(size_t i, size_t j) {
    if (i>=columns || j >= columns) {
        throw std::out_of_range("The Indexes are out of range");
    }
    for (size_t k = 0 ; k < rows ; k++) {
        std::swap(matrix[k][j] , matrix[k][i]);
    }
}


template<typename T>
void Matrix<T>::transpose() {
    Matrix<T> transposed(columns, rows);
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            transposed.matrix[j][i] = matrix[i][j];
        }
    }
    matrix = transposed;
    std::swap(rows,columns);
}

template<typename T>
bool Matrix<T>::isSquare() const {
    return rows == columns;
}

template<typename T>
bool Matrix<T>::isDiagonal() {
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            if (i != j && matrix[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::isSymmetric() {
    if (!isSquare()) return false;
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::isZero() {
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0; j < columns ; j++) {
            if (matrix[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

template<typename T>
T Matrix<T>::determinant() {
    return determinant(*this);
}

template<typename T>
Matrix<T> Matrix<T>::getInverse() {
    return Matrix<T>(rows,columns);
}

template<typename T>
void Matrix<T>::invert() {
    return;
}

template<typename T>
T Matrix<T>::trace() const {
    if (!isSquare()) {
        throw std::logic_error("Trace requires square matrix");
    }
    T sum = 0;
    for (size_t i = 0 ; i < rows ; i++) {
        sum += matrix[i][i];
    }
    return sum;
}

template<typename T>
Matrix<T> Matrix<T>::copy() {
    Matrix<T> copied(*this);
    return copied;
}

template<typename T>
Matrix<T> Matrix<T>::identity(size_t n) {
    Matrix<T> id(n,n,0);
    for (int i = 0 ; i < n ; i++) {
        id[i][i] = 1;
    }
    return id;
}

template<typename T>
T Matrix<T>::determinant(const Matrix &m) {
    if (!m.isSquare()) {
        throw std::out_of_range("The matrix has to be a square matrix");
    }
    if (m.rows == 2) {
        T a = m[0][0];
        T b = m[0][1];
        T c = m[1][0];
        T d = m[1][1];
        return (a * d - b * c);
    }
    T det = 0;
    int mult = 1;
    for (size_t j = 0 ; j < m.columns ; j++) {
        det += mult * m.get(0,j) * determinant(chop(m , 0,j));
        mult *= -1;
    }
    return det;
}

template<typename T>
size_t Matrix<T>::rank() {
    return 0;
}


template<typename T>
Matrix<T> Matrix<T>::chop(const Matrix& m, size_t row, size_t col) {
    if (m.rows <= 1 || m.columns <= 1) {
        throw std::logic_error("Cannot chop a matrix smaller than 2x2");
    }

    Matrix<T> result(m.rows - 1, m.columns - 1);
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

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    Matrix<T> result(rows , columns);
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            result.matrix[i][j] = matrix[i][j] + rhs.matrix[i][j];
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix &rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    Matrix<T> result(rows , columns);
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            result.matrix[i][j] = matrix[i][j] - rhs.matrix[i][j];
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &rhs) {
    if (columns != rhs.rows) {
        throw std::length_error("The length of the left hand matrix's cols should equal the right hand's rows");
    }
    Matrix<T> result(rows , rhs.columns);
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < rhs.columns ; j++) {
            for (size_t k = 0 ; k < rhs.rows ; k++) {
                result.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
            }
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const Matrix &rhs) {
    return Matrix<T>(rows , columns);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T &scalar) {
    Matrix<T> result(rows , columns);
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            result.matrix[i][j] = scalar * matrix[i][j] ;
        }
    }
    return result;
}

template<typename T>
void Matrix<T>::operator+=(const Matrix &rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            matrix[i][j] += rhs.matrix[i][j];
        }
    }
}

template<typename T>
void Matrix<T>::operator-=(const Matrix &rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        throw std::length_error("The length of both matrices should be the same");
    }
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            matrix[i][j] -= rhs.matrix[i][j];
        }
    }
}

template<typename T>
void Matrix<T>::operator*=(const Matrix &rhs) {
    if (columns != rhs.rows) {
        throw std::length_error("The length of the left hand matrix's cols should equal the right hand's rows");
    }
    Matrix<T> result(rows , rhs.columns);
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < rhs.columns ; j++) {
            for (size_t k = 0 ; k < rhs.rows ; k++) {
                result.matrix[i][j] += matrix[i][k] * rhs.matrix[k][j];
            }
        }
    }
    matrix = result.matrix;
    rows = result.rows;
    columns = result.columns;
}

template<typename T>
void Matrix<T>::operator*=(const T &scalar) {
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            matrix[i][j] *= scalar;
        }
    }
}

template<typename T>
bool Matrix<T>::operator==(const Matrix &rhs) {
    if (rows != rhs.rows || columns != rhs.columns) {
        return false;
    }
    for (size_t i = 0 ; i < rows ; i++) {
        for (size_t j = 0 ; j < columns ; j++) {
            if (matrix[i][j] != rhs.matrix[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template<typename T>
bool Matrix<T>::operator!=(const Matrix &rhs) {
    return !(*this == rhs);
}


// 1. Non-const version (allows modification, e.g., m[0][0] = 5)
template<typename T>
std::vector<T>& Matrix<T>::operator[](size_t i) {
    if (i >= rows) {
        throw std::out_of_range("Index out of bounds");
    }
    return matrix[i];
}

// 2. Const version (read-only, used when the Matrix itself is const)
template<typename T>
const std::vector<T>& Matrix<T>::operator[](size_t i) const {
    if (i >= rows) {
        throw std::out_of_range("Index out of bounds");
    }
    return matrix[i];
}

template<typename T>
std::string Matrix<T>::toString() const{
    std::string result;
    for (size_t i = 0 ; i < rows ; i++) {
        result += " | ";
        for (size_t j = 0 ; j < columns ; j++) {
            result += std::to_string(matrix[i][j]);
            result += " | ";
        }
        result += "\n";
    }
    return result;
}
