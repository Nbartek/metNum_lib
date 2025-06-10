#include <cmath>
#include <limits>

#include "../library.h" //to change when using
//
// Created by nack2 on 10.06.2025.
//
double f1(double x) {
    if (x <= -1.0) return NAN;
    if (x + 1.0 < std::numeric_limits<double>::epsilon() * 10) return -std::numeric_limits<double>::infinity(); // log10(0) -> -inf

    return log10(x + 1.0) - pow(x, 3);
}

double f2(double x) {
    return cosh(x) - std::fabs(x) - 0.55;
}

double f3(double x) {
    if (std::fabs(x + 1.0) < std::numeric_limits<double>::epsilon()) return NAN;
    return cos(3.0 * M_PI * x) / (x + 1.0);
}

int main() {
    nonLinearEquation::findAllRoots(f1, -3.0, 3.0, "function 1");
    nonLinearEquation::findAllRoots(f2, -3.0, 3.0, "function 2");
    nonLinearEquation::findAllRoots(f3, -3.0, 3.0, "function 3");
    return 0;
}