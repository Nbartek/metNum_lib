//
// Created by nack2 on 04.06.2025.
//

#include "../include/interpolation.h"

#include <fstream>
#include <iostream>
#include <__msvc_ostream.hpp>

std::pair<std::vector<double>, std::vector<double>> Interpolation::load_data(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile) {
        throw std::runtime_error("Error: Unable to open file " + filename);
    }

    std::vector<double> x_vals;
    std::vector<double> fx_vals;
    std::string line;
    double value;

    // Read X values
    if (getline(inputFile, line) && line.rfind("x:", 0) == 0) {
        std::istringstream ss_x(line.substr(2));
        while (ss_x >> value) {
            x_vals.push_back(value);
        }
    }

    // Read F(X) values
    if (getline(inputFile, line) && line.rfind("fx:", 0) == 0) {
        std::istringstream ss_fx(line.substr(3));
        while (ss_fx >> value) {
            fx_vals.push_back(value);
        }
    }

    return std::make_pair(x_vals, fx_vals);
}

double Interpolation::LagrangeInterpolation(const std::vector<double>& dataX, const std::vector<double>& dataY, double pointX) {
    double final_result = 0.0;
    for (size_t i = 0; i < dataX.size(); i++) {
        final_result += dataY[i] * l_basis(dataX, pointX, i);
    }
    return final_result;
}

double Interpolation::newtonInterpolation(const std::vector<double>& dataX, const std::vector<double>& dataY, double pointX) {
    double result = dataY[0];
    for (size_t i = 1; i < dataX.size(); i++) {
        result += b_divided_diff(dataX, dataY, i) * p_product(dataX, pointX, i);
    }
    return result;
}

double Interpolation::normalPolynomial(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (size_t i = 0; i < coefficients.size(); i++) {
        result += coefficients[i] * pow(x, i);
    }
    return result;
}

double Interpolation::hornerPolynomial(const std::vector<double>& coefficients, double x) {
    if (coefficients.empty()) {
        return 0.0;
    }
    double result = coefficients.back();
    for (int i = coefficients.size() - 2; i >= 0; i--) {
        result = result * x + coefficients[i];
    }
    return result;
}


// --- Private Static Helper Methods ---

double Interpolation::l_basis(const std::vector<double>& dataX, double pointX, int j) {
    double result = 1.0;
    for (size_t i = 0; i < dataX.size(); i++) {
        if (j == i) continue;
        result *= (pointX - dataX[i]) / (dataX[j] - dataX[i]);
    }
    return result;
}

double Interpolation::p_product(const std::vector<double>& dataX, double pointX, int index) {
    double result = 1.0;
    for (int i = 0; i < index; i++) {
        result *= (pointX - dataX[i]);
    }
    return result;
}

double Interpolation::b_divided_diff(const std::vector<double>& dataX, const std::vector<double>& dataY, int index) {
    double result = 0.0;
    for (int i = 0; i <= index; i++) {
        double temp_product = 1.0;
        for (int j = 0; j <= index; j++) {
            if (i == j) continue;
            temp_product *= (dataX[i] - dataX[j]);
        }
        result += dataY[i] / temp_product;
    }
    return result;
}

