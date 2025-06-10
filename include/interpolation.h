//
// Created by nack2 on 04.06.2025.
//

#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <string>
#include <vector>

class Interpolation
{
    struct pairLike {
        std::vector <double> x;
        std::vector<double> fx;
    };
    public:
    pairLike load_data(const std::string& filename);
    double l_sum(int size, std::vector<double> data, double pointX,int j);
    double LagrangeInterpolation(std::vector<double> dataX,std::vector<double> dataY, double pointX);
    long double p_counter(double pointX, const pairLike& data, int index);
    long double b_counter(double pointX, const pairLike& data, int index);
    long double normalPolynolmial(const pairLike& data,int x,int range);
    long double hornerPolinolmial(const pairLike& data, int x,int stop);
    long double newtonInterpolation(double pointX, const pairLike& data);
};
#endif //INTERPOLATION_H
