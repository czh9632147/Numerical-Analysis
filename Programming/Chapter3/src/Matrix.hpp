#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <stdexcept>
#include <cstring>
#include <cmath>
#include <algorithm>

class Matrix {
protected:
    int rows, cols;
    double* data;

public:
    // Constructors and destructor
    Matrix() : rows(0), cols(0), data(nullptr) {}
    
    explicit Matrix(int n) : rows(n), cols(n) {
        data = new double[n * n]();
    }
    
    Matrix(int r, int c) : rows(r), cols(c) {
        data = new double[r * c]();
    }
    
    // Copy constructor
    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols) {
        data = new double[rows * cols];
        std::memcpy(data, other.data, sizeof(double) * rows * cols);
    }
    
    // Move constructor
    Matrix(Matrix&& other) noexcept 
        : rows(other.rows), cols(other.cols), data(other.data) {
        other.data = nullptr;
        other.rows = other.cols = 0;
    }
    
    // Destructor
    ~Matrix() {
        delete[] data;
    }
    
    // Assignment operators
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            delete[] data;
            rows = other.rows;
            cols = other.cols;
            data = new double[rows * cols];
            std::memcpy(data, other.data, sizeof(double) * rows * cols);
        }
        return *this;
    }
    
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            delete[] data;
            rows = other.rows;
            cols = other.cols;
            data = other.data;
            other.data = nullptr;
            other.rows = other.cols = 0;
        }
        return *this;
    }
    
    // Access operators
    double* operator[](int r) {
        if (r < 0 || r >= rows) {
            throw std::out_of_range("Row index out of bounds");
        }
        return data + r * cols;
    }
    
    const double* operator[](int r) const {
        if (r < 0 || r >= rows) {
            throw std::out_of_range("Row index out of bounds");
        }
        return data + r * cols;
    }
    
    // Element access
    double& at(int r, int c) {
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::out_of_range("Index out of bounds");
        }
        return data[r * cols + c];
    }
    
    const double& at(int r, int c) const {
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::out_of_range("Index out of bounds");
        }
        return data[r * cols + c];
    }
    
    // Swap rows
    void swapRows(int r1, int r2) {
        if (r1 < 0 || r1 >= rows || r2 < 0 || r2 >= rows) {
            throw std::out_of_range("Row indices out of bounds");
        }
        for (int j = 0; j < cols; j++) {
            std::swap(data[r1 * cols + j], data[r2 * cols + j]);
        }
    }
    
    // Getters
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    
    // Solve linear system Ax = b
    static Matrix solve(const Matrix& A, const Matrix& b) {
        if (A.cols != A.rows || A.rows != b.rows || b.cols != 1) {
            throw std::invalid_argument("Invalid dimensions for linear system");
        }
        
        // Create copies for manipulation
        Matrix tmpA(A);
        Matrix tmpB(b);
        int n = A.rows;
        Matrix x(n, 1);
        
        // Gaussian elimination with partial pivoting
        for (int i = 0; i < n; i++) {
            // Find pivot
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (std::fabs(tmpA[j][i]) > std::fabs(tmpA[pivot][i])) {
                    pivot = j;
                }
            }
            
            // Swap rows if necessary
            if (pivot != i) {
                tmpA.swapRows(i, pivot);
                tmpB.swapRows(i, pivot);
            }
            
            // Check if matrix is singular
            if (std::fabs(tmpA[i][i]) < 1e-10) {
                throw std::runtime_error("Matrix is singular");
            }
            
            // Eliminate column
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                double coef = tmpA[j][i] / tmpA[i][i];
                for (int k = i; k < n; k++) {
                    tmpA[j][k] -= tmpA[i][k] * coef;
                }
                tmpB[j][0] -= tmpB[i][0] * coef;
            }
        }
        
        // Back substitution
        for (int i = 0; i < n; i++) {
            x[i][0] = tmpB[i][0] / tmpA[i][i];
        }
        
        return x;
    }
};

// Column vector class
class ColVector : public Matrix {
public:
    ColVector() : Matrix() {}
    explicit ColVector(int n) : Matrix(n, 1) {}
    
    // Element access operator
    double& operator[](int x) {
        return at(x, 0);
    }
    
    const double& operator[](int x) const {
        return at(x, 0);
    }
    
    // Copy constructor
    ColVector(const ColVector& other) : Matrix(other) {}
    
    // Move constructor
    ColVector(ColVector&& other) noexcept : Matrix(std::move(other)) {}
    
    // Assignment operators
    ColVector& operator=(const ColVector& other) {
        Matrix::operator=(other);
        return *this;
    }
    
    ColVector& operator=(ColVector&& other) noexcept {
        Matrix::operator=(std::move(other));
        return *this;
    }
};

#endif