#pragma once
#include <string>
#include <vector>
#include <ostream>

class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    size_t rows = 0;
    size_t columns = 0;

public:
    //Constructors
    Matrix() = default;
    Matrix(const Matrix& original);
    explicit Matrix(size_t rows, size_t columns, double value);
    explicit Matrix(size_t rows, size_t columns);
    Matrix(const std::vector<std::vector<double>>& matrix);
    Matrix(std::initializer_list<std::initializer_list<double>> list);

    //Public Methods
    double get(size_t i, size_t j) const;
    void set(size_t i, size_t j, double value);
    size_t getRows() const;
    size_t getColumns() const;
    void setRow(const std::vector<double>& row, size_t i);
    void setColumn(const std::vector<double>& column, size_t j);
    void swapRow(size_t i, size_t j);
    void swapColumn(size_t i, size_t j);
    void transpose();
    bool isSquare() const;
    bool isDiagonal();
    bool isSymmetric();
    bool isZero();
    double determinant();
    Matrix getInverse() const;
    Matrix getCofactor() const;
    void invert();
    double trace() const;
    Matrix copy();
    static Matrix identity(size_t n);
    static double determinantOld(const Matrix& m);
    static double determinant(const Matrix& m);
    size_t rank();
    static Matrix chop(const Matrix& m, size_t row, size_t col);
	static Matrix solveLinearSystem(const Matrix& A, const Matrix& b);
    double determinantLU() const;
    void luDecomposition(Matrix& L, Matrix& U) const;
    // Operator Overloading :
    Matrix operator+(const Matrix& rhs);
    Matrix operator-(const Matrix& rhs);
    Matrix operator*(const Matrix& rhs);
    Matrix operator/(const Matrix& rhs);
    Matrix operator/(const double& scalar);
    Matrix operator*(const double& scalar);
    void operator+=(const Matrix& rhs);
    void operator-=(const Matrix& rhs);
    void operator*=(const Matrix& rhs);
    void operator/=(const double& scalar);
    void operator*=(const double& scalar);
    bool operator==(const Matrix& rhs);
    bool operator!=(const Matrix& rhs);
    std::vector<double>& operator[](size_t i);
    const std::vector<double>& operator[](size_t i) const;
    std::string toString() const;

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
    friend Matrix operator*(const double& scalar, Matrix& m);
    friend Matrix operator/(const double& scalar, Matrix& m);

    //destructor
    ~Matrix() = default;
};