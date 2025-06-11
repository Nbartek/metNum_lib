#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include "Interpolation/Interpolation.h"

// Helper to make test output clear
void check(bool condition, const std::string& test_name) {
    if (condition) {
        std::cout << "[ PASS ] " << test_name << std::endl;
    } else {
        std::cout << "[ FAIL ] " << test_name << std::endl;
        exit(1); // Exit on first failure
    }
}

// Helper to check for floating-point equality
bool are_doubles_equal(double a, double b, double epsilon = 1e-9) {
    return std::fabs(a - b) < epsilon;
}

void test_load_data() {
    std::cout << "\n--- Testing load_data ---\n";
    // Test 1: Successful load
    auto data = Interpolation::load_data("data.txt");
    check(data.first.size() == 3, "load_data: Correct number of x values");
    check(data.second.size() == 3, "load_data: Correct number of y values");
    check(are_doubles_equal(data.first[1], 2.0), "load_data: Correct x value loaded");
    check(are_doubles_equal(data.second[2], 16.0), "load_data: Correct y value loaded");

    // Test 2: File not found
    bool exception_thrown = false;
    try {
        Interpolation::load_data("non_existent_file.txt");
    } catch (const std::runtime_error& e) {
        exception_thrown = true;
    }
    check(exception_thrown, "load_data: Throws exception for missing file");
}

void test_interpolation_methods() {
    std::cout << "\n--- Testing Lagrange and Newton Interpolation ---\n";
    // Data for f(x) = x^2
    std::vector<double> x = {1, 2, 4};
    std::vector<double> y = {1, 4, 16};

    // --- Lagrange Tests ---
    // Test 1: Interpolate at an existing data point
    double lagrange_val1 = Interpolation::LagrangeInterpolation(x, y, 2.0);
    check(are_doubles_equal(lagrange_val1, 4.0), "Lagrange: Interpolating at x=2 gives f(2)=4");

    // Test 2: Interpolate at a new point (should be exact for a quadratic)
    double lagrange_val2 = Interpolation::LagrangeInterpolation(x, y, 3.0);
    check(are_doubles_equal(lagrange_val2, 9.0), "Lagrange: Interpolating at x=3 gives f(3)=9");

    // --- Newton Tests ---
    // Test 1: Interpolate at an existing data point
    double newton_val1 = Interpolation::newtonInterpolation(x, y, 2.0);
    check(are_doubles_equal(newton_val1, 4.0), "Newton: Interpolating at x=2 gives f(2)=4");

    // Test 2: Interpolate at a new point
    double newton_val2 = Interpolation::newtonInterpolation(x, y, 3.0);
    check(are_doubles_equal(newton_val2, 9.0), "Newton: Interpolating at x=3 gives f(3)=9");
}

void test_polynomial_evaluators() {
    std::cout << "\n--- Testing Polynomial Evaluators ---\n";
    // P(x) = 5 + 3x + 2x^2. Coefficients are {5, 3, 2}
    std::vector<double> coeffs = {5.0, 3.0, 2.0};

    // --- Normal Polynomial Tests ---
    // Test 1: P(0) should be the constant term, 5
    double normal_val1 = Interpolation::normalPolynomial(coeffs, 0.0);
    check(are_doubles_equal(normal_val1, 5.0), "Normal Polynomial: P(0) = 5");

    // Test 2: P(2) = 5 + 3*2 + 2*4 = 5 + 6 + 8 = 19
    double normal_val2 = Interpolation::normalPolynomial(coeffs, 2.0);
    check(are_doubles_equal(normal_val2, 19.0), "Normal Polynomial: P(2) = 19");

    // --- Horner's Method Tests ---
    // Test 1: P(0) should be 5
    double horner_val1 = Interpolation::hornerPolynomial(coeffs, 0.0);
    check(are_doubles_equal(horner_val1, 5.0), "Horner Polynomial: P(0) = 5");

    // Test 2: P(2) should be 19
    double horner_val2 = Interpolation::hornerPolynomial(coeffs, 2.0);
    check(are_doubles_equal(horner_val2, 19.0), "Horner Polynomial: P(2) = 19");
}