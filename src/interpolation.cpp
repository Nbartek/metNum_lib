//
// Created by nack2 on 04.06.2025.
//

#include "../include/interpolation.h"

#include <fstream>
#include <iostream>
#include <__msvc_ostream.hpp>

Interpolation::pairLike Interpolation::load_data(const std::string& filename)
{
    std::ifstream inputFile(filename);
    pairLike func;
    std::string line;

    if (!inputFile) {
        std::cerr << "Error: Unable to open file " << filename << "\n";
        return func;
    }
    std::cout << filename << " is opened\n";


    getline(inputFile, line);
    std::istringstream ss_x(line.substr(4));
    double value;
    while (ss_x >> value) {
        func.x.push_back(value);
    }

    getline(inputFile, line);
    std::istringstream ss_fx(line.substr(6));
    while (ss_fx >> value) {
        func.fx.push_back(value);
    }

    return func;
}

double Interpolation::l_sum(int size, std::vector<double> data, double pointX, int j)
{
    double result = 1.0;
    for (size_t i = 0; i < size; i++)
    {
        if (j == i)continue;
        if (i == 0) {
            result = result * ((pointX - data[i]) / (data[j] - data[i]));
        }
        else {
            result = result * ((pointX - data[i]) / (data[j] - data[i]));
        }

    }
    return result;
}

double Interpolation::LagrangeInterpolation(std::vector<double> dataX, std::vector<double> dataY, double pointX)
{
    std::vector<double> l;
    double Final_result =0.0;
    double size = dataX.size();
    for (size_t i = 0; i < size; i++)
    {
        double result = l_sum(size, dataX, pointX,i);
        Final_result += result * dataY[i];
    }

    std::cout << "L(" << pointX << ") to: " << Final_result << "\n";
    return Final_result;
}
long double Interpolation::newtonInterpolation(double pointX, const pairLike& data) {

    long double result = data.fx[0]; // Start with the first fx value
    for (int i = 1; i < data.x.size(); i++) {
        long double temp = b_counter(pointX, data, i);
        result +=  temp* p_counter(pointX, data, i);
    }

    return result;
}
long double Interpolation::hornerPolinolmial(const pairLike& data, int x,int stop) {
    long double wynik = data.x[stop];

    for (int i = stop - 1; i >= 0; i--)
        wynik = wynik * x + data.x[i];

    return wynik;
}
long double Interpolation::normalPolynolmial(const pairLike& data,int x,int range) {
    long double wynik = 0.0;
    int stopien = data.x.size();

    for (int i = 0; i < range; i++) {
        wynik += data.x[i] * pow(x, i);
    }
    return wynik;
}
long double Interpolation::b_counter(double pointX, const pairLike& data, int index) {

    long double result = 0;
    for (int i = 0; i <= index; i++) {
        long double temp = 1;
        for (int j = 0; j <= index; j++)
        {
            if (i == j)continue;
            temp = temp * (data.x[i] - data.x[j]);
        }
        result += data.fx[i] / temp;

    }
    std::cout << "a" << index << ": " << result << std::endl;
    return result;
}
long double Interpolation::p_counter(double pointX, const pairLike& data, int index) {

    long double result = 1;
    for (int i = 0; i < index; i++) {
        result *= (pointX - data.x[i]);
    }

    return result;

}

