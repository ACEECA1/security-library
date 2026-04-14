#pragma once
#include <iosfwd>
#include <vector>
template<typename T> class Matrix {
private:
    std::vector<std::vector<T>> matrix;
    size_t rows = 0;
    size_t columns = 0;
public:
    //Constructors
    Matrix() = default;
    Matrix(const Matrix& original);
    Matrix(size_t rows, size_t columns);
    Matrix(size_t rows, size_t columns , T value);
    Matrix(const std::vector<std::vector<T>>& matrix);

    //Public Methods
    T get(size_t i, size_t j) const;
    void set(size_t i, size_t j, T value);
    size_t getRows() const;
    size_t getColumns() const;
    void setRow(const std::vector<T>& row , size_t i);
    void setColumn(const std::vector<T>& column, size_t j);
    void swapRow(size_t i, size_t j);
    void swapColumn(size_t i, size_t j);
    void transpose();
    bool isSquare() const;
    bool isDiagonal();
    bool isSymmetric();
    bool isZero();
    T determinant();
    Matrix getInverse();
    void invert();
    T trace() const;
    Matrix copy();
    static Matrix identity(size_t n);
    static T determinant(const Matrix& m);
    size_t rank();
    static Matrix chop(const Matrix& m, size_t row, size_t col);

    // Operator Overloading :
    Matrix operator+(const Matrix& rhs);
    Matrix operator-(const Matrix& rhs);
    Matrix operator*(const Matrix& rhs);
    Matrix operator/(const Matrix& rhs);
    Matrix operator*(const T& scalar);
    void operator+=(const Matrix& rhs);
    void operator-=(const Matrix& rhs);
    void operator*=(const Matrix& rhs);
    void operator*=(const T& scalar);
    bool operator==(const Matrix& rhs);
    bool operator!=(const Matrix& rhs);
    std::vector<T>& operator[](size_t i);
    const std::vector<T>& operator[](size_t i) const;
    std::string toString() const;

    //destructor
    ~Matrix() = default;
};

#include "../../src/Math/matrix.tpp"