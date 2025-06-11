#include "../include/linearEquation.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <numeric> // For std::iota

// --- Public Static Methods ---

std::pair<LinearEquation::Matrix, LinearEquation::Vector>
LinearEquation::loadSystemFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error: Unable to open file " + filename);
    }

    std::string line;
    int n = 0;

    // Read size
    if (std::getline(file, line) && line.rfind("size:", 0) == 0) {
        std::istringstream(line.substr(5)) >> n;
    } else {
        throw std::runtime_error("Malformed file: 'size:' line not found or incorrect.");
    }
    if (n <= 0) {
        throw std::runtime_error("Invalid matrix size in file.");
    }

    Matrix A(n, Vector(n));
    Vector b(n);

    // Read vector b
    if (std::getline(file, line) && line.rfind("b:", 0) == 0) {
        std::istringstream ss(line.substr(2));
        for (int i = 0; i < n; ++i) ss >> b[i];
    } else {
        throw std::runtime_error("Malformed file: 'b:' line not found or incorrect.");
    }

    // Read matrix A
    for (int i = 0; i < n; ++i) {
        if (std::getline(file, line) && line.rfind("A_row" + std::to_string(i), 0) == 0) {
            std::istringstream ss(line.substr(6 + std::to_string(i).length()));
            for (int j = 0; j < n; ++j) ss >> A[i][j];
        } else {
             throw std::runtime_error("Malformed file: A_row line " + std::to_string(i) + " not found or incorrect.");
        }
    }

    return {A, b};
}


std::optional<LinearEquation::Vector>
LinearEquation::solveWithGaussianElimination(Matrix A, Vector b) {
    const int n = A.size();
    if (n == 0 || A[0].size() != n || b.size() != n) {
        return std::nullopt; // Invalid input
    }

    // Forward Elimination with Partial Pivoting
    for (int i = 0; i < n; ++i) {
        int pivot_row = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(A[j][i]) > std::abs(A[pivot_row][i])) {
                pivot_row = j;
            }
        }

        if (std::abs(A[pivot_row][i]) < 1e-9) {
            return std::nullopt; // Singular or near-singular matrix
        }

        std::swap(A[i], A[pivot_row]);
        std::swap(b[i], b[pivot_row]);

        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; ++k) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // Back Substitution
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

std::optional<LinearEquation::Vector>
LinearEquation::solveWithLU(const Matrix& A, const Vector& b) {
    auto lu_result_opt = luDecomposition(A);
    if (!lu_result_opt) {
        return std::nullopt; // Decomposition failed (singular matrix)
    }

    const auto& lu = *lu_result_opt;
    const int n = A.size();
    Vector b_permuted(n);
    for (int i = 0; i < n; ++i) {
        b_permuted[i] = b[lu.P[i]];
    }

    Vector y = forwardSubstitution(lu.L, b_permuted);
    Vector x = backwardSubstitution(lu.U, y);
    return x;
}


// --- Private Static Helper Methods ---

std::optional<LinearEquation::LUResult>
LinearEquation::luDecomposition(const Matrix& A) {
    const int n = A.size();
    if (n == 0 || A[0].size() != n) {
        return std::nullopt;
    }

    LUResult res;
    res.L = Matrix(n, Vector(n, 0.0));
    res.U = Matrix(n, Vector(n, 0.0));
    res.P.resize(n);
    std::iota(res.P.begin(), res.P.end(), 0); // Fills P with 0, 1, 2, ...

    Matrix A_copy = A;

    for (int k = 0; k < n; ++k) {
        int pivot_row = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A_copy[i][k]) > std::abs(A_copy[pivot_row][k])) {
                pivot_row = i;
            }
        }
        if (std::abs(A_copy[pivot_row][k]) < 1e-9) {
             return std::nullopt; // Singular matrix
        }
        std::swap(res.P[k], res.P[pivot_row]);
        std::swap(A_copy[k], A_copy[pivot_row]);

        for (int i = k + 1; i < n; ++i) {
            A_copy[i][k] /= A_copy[k][k];
            for (int j = k + 1; j < n; ++j) {
                A_copy[i][j] -= A_copy[i][k] * A_copy[k][j];
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j) {
                res.L[i][j] = A_copy[i][j];
            } else {
                res.U[i][j] = A_copy[i][j];
                if (i == j) res.L[i][j] = 1.0;
            }
        }
    }
    return res;
}

LinearEquation::Vector
LinearEquation::forwardSubstitution(const Matrix& L, const Vector& b) {
    const int n = L.size();
    Vector y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }
    return y;
}

LinearEquation::Vector
LinearEquation::backwardSubstitution(const Matrix& U, const Vector& y) {
    const int n = U.size();
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
    return x;
}