#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <functional>
#include "../include/numericIntegration.h" // Adjust path if necessary

// --- Test Helper Utilities ---

// Helper to make test output clear and exit on failure
void check(bool condition, const std::string& test_name) {
    if (condition) {
        std::cout << "[ PASS ] " << test_name << std::endl;
    } else {
        std::cout << "[ FAIL ] " << test_name << std::endl;
        exit(1);
    }
}

// Helper to check for floating-point equality with a given tolerance
bool are_doubles_equal(double a, double b, double epsilon = 1e-6) {
    return std::fabs(a - b) < epsilon;
}

// --- Functions and Polynomials for Testing ---

// Test function 1: A simple quadratic, f(x) = x^2
// Integral from 0 to 3 is [x^3/3]_0^3 = 9.
double quadratic_func(double x) {
    return x * x;
}

// Test function 2: A trigonometric function, f(x) = sin(x)
// Integral from 0 to PI is [-cos(x)]_0^PI = 2.
double sin_func(double x) {
    return std::sin(x);
}

// Test polynomial 1: P(x) = 2x^3 - 3x^2 + 4x - 5
// Coefficients for Horner's method are {-5, 4, -3, 2}
// Integral from 1 to 3 is [0.5x^4 - x^3 + 2x^2 - 5x]_1^3 = (40.5 - 27 + 18 - 15) - (0.5 - 1 + 2 - 5) = 16.5 - (-3.5) = 20.
std::vector<double> poly_coeffs = {-5, 4, -3, 2};


// --- Test Suites ---

void test_constructor_and_gauss() {
    std::cout << "\n--- Testing Constructor and Gauss-Legendre ---" << std::endl;

    // Test 1: Constructor with a valid `n` should not throw
    try {
        numericIntegration integrator(4, {0, 1}, quadratic_func);
        check(true, "Constructor: Accepts valid n (2-5)");
    } catch (...) {
        check(false, "Constructor: Accepts valid n (2-5)");
    }

    // Test 2: Constructor with an invalid `n` should throw std::invalid_argument
    bool exception_thrown = false;
    try {
        numericIntegration integrator(6, {0, 1}, quadratic_func);
    } catch (const std::invalid_argument& e) {
        exception_thrown = true;
    }
    check(exception_thrown, "Constructor: Throws exception for invalid n (>5)");

    // Test 3: Gauss-Legendre on a polynomial it should integrate exactly
    // A 2-point rule (n=2) integrates polynomials of degree up to 3 exactly.
    // Our quadratic (degree 2) should be exact with one step (precision=1).
    numericIntegration integrator_quad(2, {0.0, 3.0}, quadratic_func);
    double result_quad = integrator_quad.gaussLegendre(1); // 1 sub-interval
    check(are_doubles_equal(result_quad, 9.0, 1e-9), "Gauss-Legendre: Integrates x^2 exactly");

    // Test 4: Gauss-Legendre on a transcendental function
    numericIntegration integrator_sin(5, {0.0, M_PI}, sin_func);
    double result_sin = integrator_sin.gaussLegendre(10); // 10 sub-intervals for precision
    check(are_doubles_equal(result_sin, 2.0, 1e-7), "Gauss-Legendre: Integrates sin(x) accurately");
}

void test_newton_cotes_methods() {
    std::cout << "\n--- Testing Newton-Cotes Methods (Rectangle, Trapezoid, Simpson) ---" << std::endl;
    const int N_STEPS = 100; // Use an even number for Simpson's rule

    // --- Polynomial Integration Tests ---
    // The specific `f` and `n` in the constructor don't matter for `...Poli` methods, but must be valid.
    numericIntegration poly_integrator(2, {1.0, 3.0}, nullptr);
    const double poly_exact_integral = 20.0;

    // Test 5: Rectangle method for polynomials
    double rect_poly_res = poly_integrator.rectangleMethodPoli(poly_coeffs, N_STEPS);
    check(are_doubles_equal(rect_poly_res, poly_exact_integral, 1e-3), "Rectangle Method (Poly): Is reasonably accurate");

    // Test 6: Trapezoid method for polynomials
    double trap_poly_res = poly_integrator.trapezeMethodPoli(poly_coeffs, N_STEPS);
    check(are_doubles_equal(trap_poly_res, poly_exact_integral, 1e-3), "Trapezoid Method (Poly): Is reasonably accurate");

    // Test 7: Simpson's method for polynomials (should be very accurate for polynomials)
    double simp_poly_res = poly_integrator.simpsonMethodPoli(poly_coeffs, N_STEPS);
    check(are_doubles_equal(simp_poly_res, poly_exact_integral, 1e-9), "Simpson's Method (Poly): Is very accurate");

    // Test 8: Simpson's method edge case (odd number of steps)
    // The implementation prints to cerr and returns 0. We test for the return value.
    double simp_odd_res = poly_integrator.simpsonMethodPoli(poly_coeffs, 99);
    check(are_doubles_equal(simp_odd_res, 0.0, 1e-9), "Simpson's Method (Poly): Handles odd N correctly (returns 0)");

    // --- Function Integration Tests ---
    numericIntegration func_integrator(2, {0.0, M_PI}, sin_func);
    const double func_exact_integral = 2.0;

    // Test 9: Rectangle method for functions
    double rect_func_res = func_integrator.rectangleMethodFun(N_STEPS);
    check(are_doubles_equal(rect_func_res, func_exact_integral, 1e-3), "Rectangle Method (Func): Is reasonably accurate");

    // Test 10: Trapezoid method for functions
    double trap_func_res = func_integrator.trapezeMethodFun(N_STEPS);
    check(are_doubles_equal(trap_func_res, func_exact_integral, 1e-3), "Trapezoid Method (Func): Is reasonably accurate");

    // Test 11: Simpson's method for functions
    double simp_func_res = func_integrator.simpsonMethodFun(N_STEPS);
    check(are_doubles_equal(simp_func_res, func_exact_integral, 1e-6), "Simpson's Method (Func): Is very accurate");
}

void test_horner_method() {
    std::cout << "\n--- Testing Horner's Method Helper ---" << std::endl;
    // P(x) = 2x^3 - 3x^2 + 4x - 5
    // P(2) = 2(8) - 3(4) + 4(2) - 5 = 16 - 12 + 8 - 5 = 7
    // P(0) = -5
    numericIntegration dummy_integrator(2, {0, 1}, nullptr); // Dummy object to call the method

    // Test 12: Horner's method with a non-zero x
    double horner_res1 = dummy_integrator.hornerPolinolmial(poly_coeffs, 2.0);
    check(are_doubles_equal(horner_res1, 7.0), "Horner's Method: Correct for P(2)");

    // Test 13: Horner's method with x=0
    double horner_res2 = dummy_integrator.hornerPolinolmial(poly_coeffs, 0.0);
    check(are_doubles_equal(horner_res2, -5.0), "Horner's Method: Correct for P(0)");
}