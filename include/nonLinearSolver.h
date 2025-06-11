      
#ifndef NONLINEAREQUATION_H
#define NONLINEAREQUATION_H

#include <functional>
#include <optional>
#include <vector>

class NonlinearSolver {
public:
    // --- Constants ---
    static constexpr double DEFAULT_TOLERANCE = 1e-9;
    static constexpr int DEFAULT_MAX_ITERATIONS = 100;

    /**
     * Finds a root of the function f within the interval [a, b] using the bisection method.
     * @param f The function for which to find a root.
     * @param a The start of the interval.
     * @param b The end of the interval.
     * @param tolerance The desired precision for the root.
     * @param max_iterations The maximum number of iterations to perform.
     * @return An optional containing the root if found, otherwise std::nullopt.
     */
    static std::optional<double> bisection(
        const std::function<double(double)>& f, double a, double b,
        double tolerance = DEFAULT_TOLERANCE, int max_iterations = DEFAULT_MAX_ITERATIONS);

    /**
     * Finds a root of the function f within [a, b] using the false position (Regula Falsi) method.
     */
    static std::optional<double> falsePosition(
        const std::function<double(double)>& f, double a, double b,
        double tolerance = DEFAULT_TOLERANCE, int max_iterations = DEFAULT_MAX_ITERATIONS);

    /**
     * Finds a root of the function f using the secant method, starting with points x0 and x1.
     */
    static std::optional<double> secant(
        const std::function<double(double)>& f, double x0, double x1,
        double tolerance = DEFAULT_TOLERANCE, int max_iterations = DEFAULT_MAX_ITERATIONS);

    /**
     * Finds a root of the function f using Newton's method with a provided analytical derivative.
     * @param f The function.
     * @param f_prime The analytical derivative of the function.
     * @param initial_guess The starting point for the search.
     */
    static std::optional<double> newton(
        const std::function<double(double)>& f, const std::function<double(double)>& f_prime,
        double initial_guess, double tolerance = DEFAULT_TOLERANCE, int max_iterations = DEFAULT_MAX_ITERATIONS);

    /**
     * Finds a root of the function f using Newton's method with a numerically calculated derivative.
     */
    static std::optional<double> newton(
        const std::function<double(double)>& f, double initial_guess,
        double tolerance = DEFAULT_TOLERANCE, int max_iterations = DEFAULT_MAX_ITERATIONS);

private:
    /**
     * Calculates the numerical derivative of a function at a point x.
     */
    static double numericalDerivative(const std::function<double(double)>& f, double x);
};

#endif //NONLINEAREQUATION_H

    