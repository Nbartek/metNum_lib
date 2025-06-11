//
// Created by nack2 on 04.06.2025.
//


#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <string>
#include <vector>
#include <utility> // Required for std::pair

class Interpolation {
public:
    /**
     * Loads x and f(x) data from a text file.
     * The file must have two lines:
     * x: 1.0 2.0 3.0
     * fx: 1.0 4.0 9.0
     * @param filename The path to the data file.
     * @return A std::pair containing the vector of x values and the vector of f(x) values.
     * @throws std::runtime_error if the file cannot be opened.
     */
    static std::pair<std::vector<double>, std::vector<double>> load_data(const std::string& filename);

    /**
     * Calculates the interpolated value at a point using the Lagrange polynomial.
     * @param dataX The vector of x coordinates.
     * @param dataY The vector of f(x) values.
     * @param pointX The x at which to interpolate.
     * @return The interpolated value L(pointX).
     */
    static double LagrangeInterpolation(const std::vector<double>& dataX, const std::vector<double>& dataY, double pointX);

    /**
     * Calculates the interpolated value at a point using Newton's divided differences polynomial.
     * @param dataX The vector of x coordinates.
     * @param dataY The vector of f(x) values.
     * @param pointX The x at which to interpolate.
     * @return The interpolated value N(pointX).
     */
    static double newtonInterpolation(const std::vector<double>& dataX, const std::vector<double>& dataY, double pointX);

    /**
     * Evaluates a polynomial P(x) = a_0 + a_1*x + a_2*x^2 + ... using a standard loop.
     * @param coefficients The vector of coefficients [a_0, a_1, a_2, ...].
     * @param x The point at which to evaluate the polynomial.
     * @return The value of the polynomial P(x).
     */
    static double normalPolynomial(const std::vector<double>& coefficients, double x);

    /**
     * Evaluates a polynomial P(x) = a_0 + a_1*x + ... using the more efficient Horner's method.
     * @param coefficients The vector of coefficients [a_0, a_1, a_2, ...].
     * @param x The point at which to evaluate the polynomial.
     * @return The value of the polynomial P(x).
     */
    static double hornerPolynomial(const std::vector<double>& coefficients, double x);

private:
    // Helper function for LagrangeInterpolation
    static double l_basis(const std::vector<double>& dataX, double pointX, int j);

    // Helper functions for newtonInterpolation
    static double p_product(const std::vector<double>& dataX, double pointX, int index);
    static double b_divided_diff(const std::vector<double>& dataX, const std::vector<double>& dataY, int index);
};

#endif //INTERPOLATION_H


