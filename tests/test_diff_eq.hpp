#include <iostream>
#include <vector>
#include <cmath>
#include "../include/differentialEquation.h" // Adjust path if necessary

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
bool are_doubles_equal(double a, double b, double epsilon) {
    return std::fabs(a - b) < epsilon;
}

// Test case: y' = y, with y(0) = 1. Analytical solution is y(t) = e^t.
auto derivative = [](double t, double y) {
    return y;
};

auto analytical_solution = [](double t) {
    return std::exp(t);
};

void test_solvers() {
    const double t_start = 0.0;
    const double t_end = 1.0;
    const double y0 = 1.0;
    const int N = 100; // Number of steps

    const double expected_y_final = analytical_solution(t_end);

    std::cout << "\n--- Testing ODE Solvers for y'=y, y(0)=1 over [0, 1] ---\n";
    std::cout << "Analytical solution at t=1: " << expected_y_final << std::endl;

    // Test 1: Euler Method
    auto euler_solution = DifferentialEquation::euler(derivative, t_start, t_end, y0, N);
    check(euler_solution.size() == N + 1, "Euler: Returns correct number of points");
    check(are_doubles_equal(euler_solution.back().second, expected_y_final, 0.02), "Euler: Final value is reasonably accurate");

    // Test 2: Heun's Method
    auto heun_solution = DifferentialEquation::heun(derivative, t_start, t_end, y0, N);
    check(heun_solution.size() == N + 1, "Heun: Returns correct number of points");
    check(are_doubles_equal(heun_solution.back().second, expected_y_final, 1e-4), "Heun: Final value is accurate");

    // Test 3: Midpoint Method
    auto midpoint_solution = DifferentialEquation::midpoint(derivative, t_start, t_end, y0, N);
    check(midpoint_solution.size() == N + 1, "Midpoint: Returns correct number of points");
    check(are_doubles_equal(midpoint_solution.back().second, expected_y_final, 1e-4), "Midpoint: Final value is accurate");

    // Test 4: Runge-Kutta 4th Order
    auto rk4_solution = DifferentialEquation::rungeKutta4(derivative, t_start, t_end, y0, N);
    check(rk4_solution.size() == N + 1, "RK4: Returns correct number of points");
    check(are_doubles_equal(rk4_solution.back().second, expected_y_final, 1e-8), "RK4: Final value is very accurate");
}

