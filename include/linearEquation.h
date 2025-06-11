#ifndef LINEAREQUATION_H
#define LINEAREQUATION_H

#include <vector>
#include <string>
#include <optional>
#include <utility>

/**
 * A static utility class for solving systems of linear equations of the form Ax = b.
 */
class LinearEquation {
public:
    // --- Type Aliases for clarity ---
    using Matrix = std::vector<std::vector<double>>;
    using Vector = std::vector<double>;

    /**
     * Loads a linear system from a file.
     * The file format must be:
     * size: N
     * b: b0 b1 ... bN-1
     * A_row0: a00 a01 ... a0N-1
     * A_row1: a10 a11 ... a1N-1
     * ...
     * @param filename The path to the data file.
     * @return A pair containing the matrix A and the vector b.
     * @throws std::runtime_error if the file cannot be opened or is malformed.
     */
    static std::pair<Matrix, Vector> loadSystemFromFile(const std::string& filename);

    /**
     * Solves the system Ax = b using Gaussian elimination with partial pivoting.
     * @param A The coefficient matrix.
     * @param b The constant vector.
     * @return An optional containing the solution vector x if a unique solution exists,
     *         otherwise std::nullopt.
     */
    static std::optional<Vector> solveWithGaussianElimination(Matrix A, Vector b);

    /**
     * Solves the system Ax = b using LU decomposition with partial pivoting.
     * @param A The coefficient matrix.
     * @param b The constant vector.
     * @return An optional containing the solution vector x if a unique solution exists,
     *         otherwise std::nullopt.
     */
    static std::optional<Vector> solveWithLU(const Matrix& A, const Vector& b);

private:
    // Private struct to hold the results of LU decomposition
    struct LUResult {
        Matrix L;
        Matrix U;
        std::vector<int> P; // Permutation vector for pivoting
    };

    static std::optional<LUResult> luDecomposition(const Matrix& A);
    static Vector forwardSubstitution(const Matrix& L, const Vector& b);
    static Vector backwardSubstitution(const Matrix& U, const Vector& y);
};

#endif //LINEAREQUATION_H