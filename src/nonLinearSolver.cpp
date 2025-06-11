#include "../include/nonLinearSolver.h"
#include <cmath>
#include <limits>

// --- Public Static Methods ---

std::optional<double> NonlinearSolver::bisection(const std::function<double(double)>& f, double a, double b, double tolerance, int max_iterations) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb >= 0) {
        // Root is not bracketed or one of the endpoints is a root.
        if (std::abs(fa) < tolerance) return a;
        if (std::abs(fb) < tolerance) return b;
        return std::nullopt; // No sign change, method cannot proceed.
    }

    for (int i = 0; i < max_iterations; ++i) {
        double mid = a + (b - a) / 2.0; // Avoids potential overflow of (a+b)
        double f_mid = f(mid);

        if (std::abs(f_mid) < tolerance || (b - a) / 2.0 < tolerance) {
            return mid;
        }

        if (f(a) * f_mid < 0) {
            b = mid;
        } else {
            a = mid;
        }
    }
    return std::nullopt; // Max iterations reached
}

std::optional<double> NonlinearSolver::falsePosition(const std::function<double(double)>& f, double a, double b, double tolerance, int max_iterations) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb >= 0) {
        if (std::abs(fa) < tolerance) return a;
        if (std::abs(fb) < tolerance) return b;
        return std::nullopt;
    }

    double x = a;
    for (int i = 0; i < max_iterations; ++i) {
        double x_next = (a * f(b) - b * f(a)) / (f(b) - f(a));
        double fx_next = f(x_next);

        if (std::abs(fx_next) < tolerance || std::abs(x_next - x) < tolerance) {
            return x_next;
        }
        x = x_next;

        if (f(a) * fx_next < 0) {
            b = x_next;
        } else {
            a = x_next;
        }
    }
    return std::nullopt; // Max iterations reached
}

std::optional<double> NonlinearSolver::secant(const std::function<double(double)>& f, double x0, double x1, double tolerance, int max_iterations) {
    for (int i = 0; i < max_iterations; ++i) {
        double fx0 = f(x0);
        double fx1 = f(x1);

        if (std::abs(fx1) < tolerance) {
            return x1;
        }

        double denominator = fx1 - fx0;
        if (std::abs(denominator) < std::numeric_limits<double>::epsilon()) {
            return std::nullopt; // Division by zero or stalling
        }

        double x2 = x1 - fx1 * (x1 - x0) / denominator;
        if (std::abs(x2 - x1) < tolerance) {
            return x2;
        }

        x0 = x1;
        x1 = x2;
    }
    return std::nullopt; // Max iterations reached
}

std::optional<double> NonlinearSolver::newton(const std::function<double(double)>& f, const std::function<double(double)>& f_prime, double initial_guess, double tolerance, int max_iterations) {
    double x0 = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        double fx = f(x0);
        double f_prime_x = f_prime(x0);

        if (std::abs(f_prime_x) < std::numeric_limits<double>::epsilon()) {
            return std::nullopt; // Derivative is zero, cannot proceed.
        }

        double x1 = x0 - fx / f_prime_x;

        if (std::abs(x1 - x0) < tolerance || std::abs(f(x1)) < tolerance) {
            return x1;
        }
        x0 = x1;
    }
    return std::nullopt; // Max iterations reached
}

std::optional<double> NonlinearSolver::newton(const std::function<double(double)>& f, double initial_guess, double tolerance, int max_iterations) {
    auto f_prime = [&](double x) { return numericalDerivative(f, x); };
    return newton(f, f_prime, initial_guess, tolerance, max_iterations);
}


// --- Private Static Helper Methods ---

double NonlinearSolver::numericalDerivative(const std::function<double(double)>& f, double x) {
    const double h = 1e-7;
    return (f(x + h) - f(x - h)) / (2.0 * h);
}