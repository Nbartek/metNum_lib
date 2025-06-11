#include "../include/differentialEquation.h"
#include <stdexcept>

std::vector<DifferentialEquation::Point> DifferentialEquation::euler(
    const std::function<double(double, double)>& derivative,
    double t_start, double t_end, double y0, int N) {

    if (N <= 0) {
        return {{t_start, y0}};
    }

    std::vector<Point> solution;
    solution.reserve(N + 1);
    solution.push_back({t_start, y0});

    const double h = (t_end - t_start) / N;
    double t = t_start;
    double y = y0;

    for (int i = 0; i < N; ++i) {
        y = y + h * derivative(t, y);
        t = t + h;
        solution.push_back({t, y});
    }
    return solution;
}

std::vector<DifferentialEquation::Point> DifferentialEquation::heun(
    const std::function<double(double, double)>& derivative,
    double t_start, double t_end, double y0, int N) {

    if (N <= 0) {
        return {{t_start, y0}};
    }

    std::vector<Point> solution;
    solution.reserve(N + 1);
    solution.push_back({t_start, y0});

    const double h = (t_end - t_start) / N;
    double t = t_start;
    double y = y0;

    for (int i = 0; i < N; ++i) {
        double k1 = derivative(t, y);
        double y_predictor = y + h * k1;
        double k2 = derivative(t + h, y_predictor);
        y = y + (h / 2.0) * (k1 + k2);
        t = t + h;
        solution.push_back({t, y});
    }
    return solution;
}

std::vector<DifferentialEquation::Point> DifferentialEquation::rungeKutta4(
    const std::function<double(double, double)>& derivative,
    double t_start, double t_end, double y0, int N) {

    if (N <= 0) {
        return {{t_start, y0}};
    }

    std::vector<Point> solution;
    solution.reserve(N + 1);
    solution.push_back({t_start, y0});

    const double h = (t_end - t_start) / N;
    double t = t_start;
    double y = y0;

    for (int i = 0; i < N; ++i) {
        double k1 = h * derivative(t, y);
        double k2 = h * derivative(t + h / 2.0, y + k1 / 2.0);
        double k3 = h * derivative(t + h / 2.0, y + k2 / 2.0);
        double k4 = h * derivative(t + h, y + k3);
        y = y + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t = t + h;
        solution.push_back({t, y});
    }
    return solution;
}

std::vector<DifferentialEquation::Point> DifferentialEquation::midpoint(
    const std::function<double(double, double)>& derivative,
    double t_start, double t_end, double y0, int N) {

    if (N <= 0) {
        return {{t_start, y0}};
    }

    std::vector<Point> solution;
    solution.reserve(N + 1);
    solution.push_back({t_start, y0});

    const double h = (t_end - t_start) / N;
    double t = t_start;
    double y = y0;

    for (int i = 0; i < N; ++i) {
        double k1 = derivative(t, y);
        double k2 = derivative(t + h / 2.0, y + (h / 2.0) * k1);
        y = y + h * k2;
        t = t + h;
        solution.push_back({t, y});
    }
    return solution;
}