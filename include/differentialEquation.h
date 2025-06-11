      
#ifndef DIFFERENTIALEQUATION_H
#define DIFFERENTIALEQUATION_H

#include <functional>
#include <vector>
#include <utility>

/**
 * A static utility class for solving first-order ordinary differential equations (ODEs)
 * of the form y' = f(t, y).
 */
class DifferentialEquation {
public:
    // A type alias for a point in the solution (time, value)
    using Point = std::pair<double, double>;

    /**
     * Solves the ODE using the explicit Euler method.
     * @param derivative The function f(t, y).
     * @param t_start The initial time.
     * @param t_end The final time.
     * @param y0 The initial value y(t_start).
     * @param N The number of steps to take.
     * @return A vector of (t, y) points representing the solution path.
     */
    static std::vector<Point> euler(
        const std::function<double(double, double)>& derivative,
        double t_start, double t_end, double y0, int N);

    /**
     * Solves the ODE using Heun's method (an improved Euler/predictor-corrector method).
     */
    static std::vector<Point> heun(
        const std::function<double(double, double)>& derivative,
        double t_start, double t_end, double y0, int N);

    /**
     * Solves the ODE using the classical 4th-order Runge-Kutta method (RK4).
     */
    static std::vector<Point> rungeKutta4(
        const std::function<double(double, double)>& derivative,
        double t_start, double t_end, double y0, int N);

    /**
     * Solves the ODE using the midpoint method.
     */
    static std::vector<Point> midpoint(
        const std::function<double(double, double)>& derivative,
        double t_start, double t_end, double y0, int N);
};

#endif //DIFFERENTIALEQUATION_H

    