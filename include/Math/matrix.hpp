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
     * @param original The Matrix to copy.
     */
    Matrix(const Matrix& original);

    /**
     * @brief Constructs a matrix with the given size and fills it with a specific value.
     * @param rows Number of rows.
     * @param columns Number of columns.
     * @param value The value to fill all elements with.
     */
    explicit Matrix(size_t rows, size_t columns, double value);

    /**
     * @brief Constructs a matrix with the given size filled with zeroes.
     * @param rows Number of rows.
     * @param columns Number of columns.
     */
    explicit Matrix(size_t rows, size_t columns);

    /**
     * @brief Constructs a matrix from a 2D vector.
     * @param matrix Input 2D vector containing matrix values.
     * @throws std::invalid_argument if rows have inconsistent lengths.
     */
    explicit Matrix(const std::vector<std::vector<double>>& matrix);

    /**
     * @brief Constructs a matrix from nested initializer lists.
     * @param list Nested initializer lists representing rows and columns.
     * @throws std::invalid_argument if rows have inconsistent lengths.
     */
    Matrix(std::initializer_list<std::initializer_list<double>> list);

    /**
     * @brief Retrieves the value at a specific row and column.
     * @param i Row index.
     * @param j Column index.
     * @return The double value at the specified position.
     * @throws std::out_of_range if indices are out of bounds.
     */
    double get(size_t i, size_t j) const;

    /**
     * @brief Sets the value at a specific row and column.
     * @param i Row index.
     * @param j Column index.
     * @param value The value to set.
     * @throws std::out_of_range if indices are out of bounds.
     */
    void set(size_t i, size_t j, double value);

    /**
     * @brief Gets the number of rows in the matrix.
     * @return The number of rows.
     */
    size_t getRows() const;

    /**
     * @brief Gets the number of columns in the matrix.
     * @return The number of columns.
     */
    size_t getColumns() const;

    /**
     * @brief Replaces an entire row with a new vector of values.
     * @param row A vector containing the new row values.
     * @param i The index of the row to replace.
     * @throws std::out_of_range if the index is out of bounds.
     * @throws std::length_error if the row length does not match the column count.
     */
    void setRow(const std::vector<double>& row, size_t i);

    /**
     * @brief Replaces an entire column with a new vector of values.
     * @param column A vector containing the new column values.
     * @param j The index of the column to replace.
     * @throws std::out_of_range if the index is out of bounds.
     * @throws std::length_error if the column length does not match the row count.
     */
    void setColumn(const std::vector<double>& column, size_t j);

    /**
     * @brief Swaps two rows within the matrix.
     * @param i Index of the first row.
     * @param j Index of the second row.
     * @throws std::out_of_range if indices are out of bounds.
     */
    void swapRow(size_t i, size_t j);

    /**
     * @brief Swaps two columns within the matrix.
     * @param i Index of the first column.
     * @param j Index of the second column.
     * @throws std::out_of_range if indices are out of bounds.
     */
    void swapColumn(size_t i, size_t j);

    /**
     * @brief Transposes the matrix in place (swaps rows and columns).
     */
    void transpose();

    /**
     * @brief Checks if the matrix is square (rows == columns).
     * @return True if square, false otherwise.
     */
    bool isSquare() const;

    /**
     * @brief Checks if the matrix is a diagonal matrix.
     * @return True if all non-diagonal elements are close to zero, false otherwise.
     */
    bool isDiagonal() const;

    /**
     * @brief Checks if the matrix is symmetric (equals its transpose).
     * @return True if symmetric, false otherwise.
     */
    bool isSymmetric() const;

    /**
     * @brief Checks if all elements in the matrix are zero.
     * @return True if all elements are close to zero, false otherwise.
     */
    bool isZero() const;

    /**
     * @brief Calculates the determinant of the matrix.
     * @return The determinant value.
     */
    double determinant() const;

    /**
     * @brief Computes and returns the inverse of the matrix.
     * @return A new Matrix representing the inverse.
     * @throws std::logic_error if the matrix is not square or is singular (non-invertible).
     */
    Matrix getInverse() const;

    /**
     * @brief Computes and returns the cofactor matrix.
     * @return A new Matrix representing the cofactor matrix.
     * @throws std::logic_error if the matrix is not square.
     */
    Matrix getCofactor() const;

    /**
     * @brief Inverts the matrix in place.
     * @throws std::logic_error if the matrix is not square or is singular.
     */
    void invert();

    /**
     * @brief Calculates the trace of the matrix (sum of main diagonal elements).
     * @return The trace value.
     * @throws std::logic_error if the matrix is not square.
     */
    double trace() const;

    /**
     * @brief Creates and returns a deep copy of the matrix.
     * @return A copy of the current matrix.
     */
    Matrix copy() const;

    /**
     * @brief Generates an identity matrix of size n x n.
     * @param n The size (rows and columns) of the identity matrix.
     * @return An n x n identity Matrix.
     */
    static Matrix identity(size_t n);

    /**
     * @brief Calculates the determinant using the older Laplace expansion method.
     * @param m The matrix to evaluate.
     * @return The determinant value.
     * @throws std::out_of_range if the matrix is not square.
     */
    static double determinantOld(const Matrix& m);

    /**
     * @brief Static helper to calculate the determinant using LU decomposition.
     * @param m The matrix to evaluate.
     * @return The determinant value.
     */
    static double determinant(const Matrix& m);

    /**
     * @brief Estimates the rank of the matrix via LU decomposition.
     * @return The estimated rank (number of linearly independent rows).
     */
    size_t rank() const;

    /**
     * @brief Creates a submatrix by removing a specific row and column.
     * @param m The original matrix.
     * @param row The index of the row to remove.
     * @param col The index of the column to remove.
     * @return A new smaller Matrix.
     * @throws std::logic_error if the matrix is smaller than 2x2.
     * @throws std::out_of_range if row or col indices are invalid.
     */
    static Matrix chop(const Matrix& m, size_t row, size_t col);

    /**
     * @brief Solves a linear system of equations Ax = b.
     * @param A The coefficient matrix (must be square).
     * @param b The constants column vector.
     * @return A column vector matrix x representing the solution.
     * @throws std::logic_error on dimension mismatch or if A is singular.
     */
    static Matrix solveLinearSystem(const Matrix& A, const Matrix& b);

    /**
     * @brief Calculates the determinant using LU decomposition internally.
     * @return The determinant value.
     * @throws std::logic_error if the matrix is not square.
     */
    double determinantLU() const;

    /**
     * @brief Performs LU decomposition (A = LU) and tracks row swaps.
     * @param L The matrix to store the Lower triangular result.
     * @param U The matrix to store the Upper triangular result.
     * @param swapCount Reference to a counter for row swaps (used for determinant sign).
     * @throws std::logic_error if the matrix is not square or is singular.
     */
    void luDecomposition(Matrix& L, Matrix& U, size_t& swapCount) const;

    /**
     * @brief Performs LU decomposition (A = LU) without tracking row swaps.
     * @param L The matrix to store the Lower triangular result.
     * @param U The matrix to store the Upper triangular result.
     * @throws std::logic_error if the matrix is not square or is singular.
     */
    void luDecomposition(Matrix& L, Matrix& U) const;

    /**
     * @brief Adds another matrix to this one.
     * @param rhs The matrix to add.
     * @return A new Matrix containing the sum.
     * @throws std::length_error if matrix dimensions do not match.
     */
    Matrix operator+(const Matrix& rhs) const;

    /**
     * @brief Subtracts another matrix from this one.
     * @param rhs The matrix to subtract.
     * @return A new Matrix containing the difference.
     * @throws std::length_error if matrix dimensions do not match.
     */
    Matrix operator-(const Matrix& rhs) const;

    /**
     * @brief Multiplies this matrix by another matrix.
     * @param rhs The matrix to multiply by.
     * @return A new Matrix containing the product.
     * @throws std::length_error if columns of LHS do not match rows of RHS.
     */
    Matrix operator*(const Matrix& rhs) const;

    /**
     * @brief Multiplies this matrix by the inverse of another matrix (A * B^-1).
     * @param rhs The matrix to divide by.
     * @return A new Matrix containing the result.
     */
    Matrix operator/(const Matrix& rhs) const;

    /**
     * @brief Divides all elements in the matrix by a scalar.
     * @param scalar The scalar value to divide by.
     * @return A new Matrix with scaled elements.
     * @throws std::domain_error if the scalar is zero.
     */
    Matrix operator/(const double& scalar) const;

    /**
     * @brief Multiplies all elements in the matrix by a scalar.
     * @param scalar The scalar value to multiply by.
     * @return A new Matrix with scaled elements.
     */
    Matrix operator*(const double& scalar) const;

    /**
     * @brief Adds another matrix to this matrix in-place.
     * @param rhs The matrix to add.
     * @throws std::length_error if matrix dimensions do not match.
     */
    void operator+=(const Matrix& rhs);

    /**
     * @brief Subtracts another matrix from this matrix in-place.
     * @param rhs The matrix to subtract.
     * @throws std::length_error if matrix dimensions do not match.
     */
    void operator-=(const Matrix& rhs);

    /**
     * @brief Multiplies this matrix by another matrix in-place.
     * @param rhs The matrix to multiply by.
     * @throws std::length_error if columns of LHS do not match rows of RHS.
     */
    void operator*=(const Matrix& rhs);

    /**
     * @brief Divides all elements in this matrix by a scalar in-place.
     * @param scalar The scalar value to divide by.
     * @throws std::domain_error if the scalar is zero.
     */
    void operator/=(const double& scalar);

    /**
     * @brief Multiplies all elements in this matrix by a scalar in-place.
     * @param scalar The scalar value to multiply by.
     */
    void operator*=(const double& scalar);

    /**
     * @brief Checks if two matrices are equal (within epsilon bounds).
     * @param rhs The matrix to compare against.
     * @return True if dimensions and all elements match, false otherwise.
     */
    bool operator==(const Matrix& rhs) const;

    /**
     * @brief Checks if two matrices are not equal.
     * @param rhs The matrix to compare against.
     * @return True if dimensions or any elements differ, false otherwise.
     */
    bool operator!=(const Matrix& rhs) const;

    /**
     * @brief Accesses a specific row in the matrix (mutable).
     * @param i The row index.
     * @return A reference to the vector representing the row.
     * @throws std::out_of_range if the index is out of bounds.
     */
    std::vector<double>& operator[](size_t i);

    /**
     * @brief Accesses a specific row in the matrix (read-only).
     * @param i The row index.
     * @return A const reference to the vector representing the row.
     * @throws std::out_of_range if the index is out of bounds.
     */
    const std::vector<double>& operator[](size_t i) const;

    /**
     * @brief Returns a formatted string representation of the matrix.
     * @return A multi-line string showing the matrix grid.
     */
    std::string toString() const;

    /**
     * @brief Stream insertion operator to print the matrix formatting.
     * @param os The output stream.
     * @param m The matrix to print.
     * @return The modified output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

    /**
     * @brief Multiplies a scalar by a matrix (commutative multiplication).
     * @param scalar The scalar value.
     * @param m The matrix.
     * @return A new Matrix with scaled elements.
     */
    friend Matrix operator*(const double& scalar, const Matrix& m);

    /**
     * @brief Default destructor.
     */
    ~Matrix() = default;
};