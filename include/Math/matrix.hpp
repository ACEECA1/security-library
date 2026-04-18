#pragma once
#include <string>
#include <vector>
#include <ostream>
#include <initializer_list>

/**
 * @brief A lightweight matrix class for basic linear algebra operations.
 *
 * Stores values in row-major order and supports matrix arithmetic,
 * decomposition, inversion, rank estimation, and formatting.
 */
class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    size_t rows = 0;
    size_t columns = 0;

    static constexpr double EPSILON = 1e-12;

public:
    /**
     * @brief Constructs an empty 0x0 matrix.
     */
    Matrix() = default;

    /**
     * @brief Copy-constructs a matrix.
     * @param original Matrix to copy.
     */
    Matrix(const Matrix& original);

    /**
     * @brief Constructs a matrix with the given size and fills it with a value.
     * @param rows Number of rows.
     * @param columns Number of columns.
     * @param value Fill value.
     */
    explicit Matrix(size_t rows, size_t columns, double value);

    /**
     * @brief Constructs a rows x columns matrix filled with zeroes.
     * @param rows Number of rows.
     * @param columns Number of columns.
     */
    explicit Matrix(size_t rows, size_t columns);

    /**
     * @brief Constructs a matrix from a 2D vector.
     * @param matrix Input values.
     * @throws std::invalid_argument if rows have inconsistent lengths.
     */
    explicit Matrix(const std::vector<std::vector<double>>& matrix);

    /**
     * @brief Constructs a matrix from nested initializer lists.
     * @param list Input values.
     * @throws std::invalid_argument if rows have inconsistent lengths.
     */
    Matrix(std::initializer_list<std::initializer_list<double>> list);

    /**
     * @brief Returns the value at the given position.
     * @param i Row index.
     * @param j Column index.
     * @return Value at (i, j).
     * @throws std::out_of_range if indices are invalid.
     */
    double get(size_t i, size_t j) const;

    /**
     * @brief Sets the value at the given position.
     * @param i Row index.
     * @param j Column index.
     * @param value New value.
     * @throws std::out_of_range if indices are invalid.
     */
    void set(size_t i, size_t j, double value);

    /**
     * @brief Returns the number of rows.
     */
    size_t getRows() const;

    /**
     * @brief Returns the number of columns.
     */
    size_t getColumns() const;

    /**
     * @brief Replaces a row.
     * @param row Replacement row.
     * @param i Row index.
     * @throws std::out_of_range if row index is invalid.
     * @throws std::length_error if row length is incorrect.
     */
    void setRow(const std::vector<double>& row, size_t i);

    /**
     * @brief Replaces a column.
     * @param column Replacement column.
     * @param j Column index.
     * @throws std::out_of_range if column index is invalid.
     * @throws std::length_error if column length is incorrect.
     */
    void setColumn(const std::vector<double>& column, size_t j);

    /**
     * @brief Swaps two rows.
     * @param i First row index.
     * @param j Second row index.
     * @throws std::out_of_range if either index is invalid.
     */
    void swapRow(size_t i, size_t j);

    /**
     * @brief Swaps two columns.
     * @param i First column index.
     * @param j Second column index.
     * @throws std::out_of_range if either index is invalid.
     */
    void swapColumn(size_t i, size_t j);

    /**
     * @brief Transposes the matrix in place.
     */
    void transpose();

    /**
     * @brief Returns whether the matrix is square.
     */
    bool isSquare() const;

    /**
     * @brief Returns whether the matrix is diagonal.
     */
    bool isDiagonal() const;

    /**
     * @brief Returns whether the matrix is symmetric.
     */
    bool isSymmetric() const;

    /**
     * @brief Returns whether the matrix is all zeros.
     */
    bool isZero() const;

    /**
     * @brief Returns the determinant of the matrix.
     * @throws std::logic_error if the matrix is not square.
     */
    double determinant() const;

    /**
     * @brief Returns the inverse matrix.
     * @throws std::logic_error if the matrix is not invertible.
     */
    Matrix getInverse() const;

    /**
     * @brief Returns the cofactor matrix.
     * @throws std::logic_error if the matrix is not square.
     */
    Matrix getCofactor() const;

    /**
     * @brief Replaces this matrix with its inverse.
     * @throws std::logic_error if the matrix is not invertible.
     */
    void invert();

    /**
     * @brief Returns the trace of a square matrix.
     * @throws std::logic_error if the matrix is not square.
     */
    double trace() const;

    /**
     * @brief Returns a deep copy of the matrix.
     */
    Matrix copy() const;

    /**
     * @brief Returns an identity matrix of size n.
     * @param n Matrix size.
     */
    static Matrix identity(size_t n);

    /**
     * @brief Computes the determinant using the older recursive method.
     * @param m Input matrix.
     * @throws std::out_of_range if the matrix is not square.
     */
    static double determinantOld(const Matrix& m);

    /**
     * @brief Computes the determinant using the LU-based implementation.
     * @param m Input matrix.
     */
    static double determinant(const Matrix& m);

    /**
     * @brief Returns the estimated rank of the matrix.
     */
    size_t rank() const;

    /**
     * @brief Returns the matrix with one row and one column removed.
     * @param m Source matrix.
     * @param row Row to remove.
     * @param col Column to remove.
     * @throws std::out_of_range if row or column is invalid.
     * @throws std::logic_error if matrix is too small.
     */
    static Matrix chop(const Matrix& m, size_t row, size_t col);

    /**
     * @brief Solves a linear system A x = b.
     * @param A Coefficient matrix.
     * @param b Right-hand side column vector.
     * @return Solution vector.
     * @throws std::logic_error if dimensions are invalid or matrix is singular.
     */
    static Matrix solveLinearSystem(const Matrix& A, const Matrix& b);

    /**
     * @brief Returns the determinant computed via LU decomposition.
     */
    double determinantLU() const;

    /**
     * @brief Performs LU decomposition with partial pivoting.
     * @param L Lower triangular matrix output.
     * @param U Upper triangular matrix output.
     * @param swapCount Number of row swaps performed.
     * @throws std::logic_error if the matrix is not square or singular.
     */
    void luDecomposition(Matrix& L, Matrix& U, size_t& swapCount) const;

    /**
     * @brief Performs LU decomposition with partial pivoting.
     * @param L Lower triangular matrix output.
     * @param U Upper triangular matrix output.
     * @throws std::logic_error if the matrix is not square or singular.
     */
    void luDecomposition(Matrix& L, Matrix& U) const;

    // Operator Overloading :

    /**
     * @brief Adds two matrices.
     */
    Matrix operator+(const Matrix& rhs) const;

    /**
     * @brief Subtracts two matrices.
     */
    Matrix operator-(const Matrix& rhs) const;

    /**
     * @brief Multiplies two matrices.
     */
    Matrix operator*(const Matrix& rhs) const;

    /**
     * @brief Divides by another matrix by multiplying with its inverse.
     */
    Matrix operator/(const Matrix& rhs) const;

    /**
     * @brief Divides by a scalar.
     * @throws std::domain_error if scalar is zero.
     */
    Matrix operator/(const double& scalar) const;

    /**
     * @brief Multiplies by a scalar.
     */
    Matrix operator*(const double& scalar) const;

    /**
     * @brief Adds another matrix in place.
     */
    void operator+=(const Matrix& rhs);

    /**
     * @brief Subtracts another matrix in place.
     */
    void operator-=(const Matrix& rhs);

    /**
     * @brief Multiplies by another matrix in place.
     */
    void operator*=(const Matrix& rhs);

    /**
     * @brief Divides by a scalar in place.
     * @throws std::domain_error if scalar is zero.
     */
    void operator/=(const double& scalar);

    /**
     * @brief Multiplies by a scalar in place.
     */
    void operator*=(const double& scalar);

    /**
     * @brief Compares matrices with floating-point tolerance.
     */
    bool operator==(const Matrix& rhs) const;

    /**
     * @brief Returns the negation of operator==.
     */
    bool operator!=(const Matrix& rhs) const;

    /**
     * @brief Returns a mutable row reference.
     * @throws std::out_of_range if the row index is invalid.
     */
    std::vector<double>& operator[](size_t i);

    /**
     * @brief Returns a const row reference.
     * @throws std::out_of_range if the row index is invalid.
     */
    const std::vector<double>& operator[](size_t i) const;

    /**
     * @brief Returns a formatted string representation of the matrix.
     */
    std::string toString() const;

    /**
     * @brief Streams a formatted matrix to an output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

    /**
     * @brief Multiplies a matrix by a scalar.
     */
    friend Matrix operator*(const double& scalar, const Matrix& m);

    /**
     * @brief Destroys the matrix.
     */
    ~Matrix() = default;
};