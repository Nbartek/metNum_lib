#include <iostream>
#include <vector>
#include <cmath>
#include <optional>
#include "../include/linearEquation.h" // Adjust path if necessary

// Helper to make test output clear
void check(bool condition, const std::string& test_name) {
    if (condition) {
        std::cout << "[ PASS ] " << test_name << std::endl;
    } else {
        std::cout << "[ FAIL ] " << test_name << std::endl;
        exit(1); // Exit on first failure
    }
}

// Helper to check for vector equality
bool are_vectors_equal(const std::vector<double>& a, const std::vector<double>& b, double epsilon = 1e-7) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > epsilon) return false;
    }
    return true;
}

void test_file_loading() {
    std::cout << "\n--- Testing File Loading ---\n";
    // Test 1: Successful load
    auto system = LinearEquation::loadSystemFromFile("code.txt");
    check(system.first.size() == 3, "Load: Correct matrix size");
    check(system.second.size() == 3, "Load: Correct vector size");
    check(system.first[2][1] == 5.0, "Load: Correct matrix element loaded");
    check(system.second[1] == -4.0, "Load: Correct vector element loaded");

    // Test 2: File not found
    bool exception_thrown = false;
    try {
        LinearEquation::loadSystemFromFile("non_existent_file.txt");
    } catch (const std::runtime_error& e) {
        exception_thrown = true;
    }
    check(exception_thrown, "Load: Throws exception for missing file");
}

void test_solvers() {
    std::cout << "\n--- Testing Linear System Solvers ---\n";
    auto system = LinearEquation::loadSystemFromFile("matrix_data.txt");
    LinearEquation::Matrix A = system.first;
    LinearEquation::Vector b = system.second;
    LinearEquation::Vector expected_x = {5.0, 3.0, -2.0};

    // Test 1: Gaussian Elimination
    auto x_gauss = LinearEquation::solveWithGaussianElimination(A, b);
    check(x_gauss.has_value(), "Gauss: Solver returns a value");
    check(are_vectors_equal(x_gauss.value(), expected_x), "Gauss: Solution is correct");

    // Test 2: LU Decomposition
    auto x_lu = LinearEquation::solveWithLU(A, b);
    check(x_lu.has_value(), "LU: Solver returns a value");
    check(are_vectors_equal(x_lu.value(), expected_x), "LU: Solution is correct");

    // Test 3: Singular Matrix
    LinearEquation::Matrix singular_A = {{1, 1, 1}, {2, 2, 2}, {3, 4, 5}};
    LinearEquation::Vector singular_b = {1, 2, 3};
    auto x_singular_gauss = LinearEquation::solveWithGaussianElimination(singular_A, singular_b);
    auto x_singular_lu = LinearEquation::solveWithLU(singular_A, singular_b);
    check(!x_singular_gauss.has_value(), "Gauss: Fails correctly for singular matrix");
    check(!x_singular_lu.has_value(), "LU: Fails correctly for singular matrix");
}