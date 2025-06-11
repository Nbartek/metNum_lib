#include <iostream>
#include <cmath>
#include <optional>
#include "NonlinearSolver/NonlinearSolver.h"

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
bool are_doubles_equal(double a, double b, double epsilon = 1e-7) {
    return std::fabs(a - b) < epsilon;
}

// --- Test Functions ---
// f(x) = x^2 - 4. Roots: -2, 2
double poly_func(double x) { return x * x - 4.0; }
double poly_func_prime(double x) { return 2.0 * x; }

// f(x) = cos(x) - x. Root near 0.739
double trig_func(double x) { return std::cos(x) - x; }

void test_bisection() {
    std::cout << "\n--- Testing Bisection Method ---\n";
    // Test 1: Find a known root of x^2 - 4
    auto root1 = NonlinearSolver::bisection(poly_func, 0, 3);
    check(root1.has_value() && are_doubles_equal(root1.value(), 2.0), "Bisection: Finds root x=2 for x^2-4");

    // Test 2: Fails when no root is bracketed
    auto root2 = NonlinearSolver::bisection(poly_func, 2.1, 5);
    check(!root2.has_value(), "Bisection: Fails correctly with no sign change");
}

void test_false_position() {
    std::cout << "\n--- Testing False Position Method ---\n";
    // Test 1: Find a known root
    auto root1 = NonlinearSolver::falsePosition(poly_func, 0, 3);
    check(root1.has_value() && are_doubles_equal(root1.value(), 2.0), "False Position: Finds root x=2 for x^2-4");

    // Test 2: Fails when no root is bracketed
    auto root2 = NonlinearSolver::falsePosition(poly_func, 2.1, 5);
    check(!root2.has_value(), "False Position: Fails correctly with no sign change");
}

void test_secant() {
    std::cout << "\n--- Testing Secant Method ---\n";
    // Test 1: Find a known root
    auto root1 = NonlinearSolver::secant(poly_func, 0, 3);
    check(root1.has_value() && are_doubles_equal(root1.value(), 2.0), "Secant: Finds root x=2 for x^2-4");

    // Test 2: Find a transcendental root
    auto root2 = NonlinearSolver::secant(trig_func, 0, 1);
    check(root2.has_value() && are_doubles_equal(root2.value(), 0.739085), "Secant: Finds root for cos(x)-x");
}

void test_newton() {
    std::cout << "\n--- Testing Newton's Method ---\n";
    // Test 1: With analytical derivative
    auto root1 = NonlinearSolver::newton(poly_func, poly_func_prime, 3.0);
    check(root1.has_value() && are_doubles_equal(root1.value(), 2.0), "Newton (analytical): Finds root x=2 for x^2-4");

    // Test 2: With numerical derivative
    auto root2 = NonlinearSolver::newton(poly_func, -3.0);
    check(root2.has_value() && are_doubles_equal(root2.value(), -2.0), "Newton (numerical): Finds root x=-2 for x^2-4");

    // Test 3: Fails when derivative is zero
    auto root3 = NonlinearSolver::newton(poly_func, poly_func_prime, 0.0);
    check(!root3.has_value(), "Newton (analytical): Fails correctly when derivative is zero");
}