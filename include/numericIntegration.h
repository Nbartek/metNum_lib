//
// Created by nack2 on 04.06.2025.
//

#ifndef NUMERICINTEGRATION_H
#define NUMERICINTEGRATION_H
#include <functional>
#include <utility>
#include <vector>

class numericIntegration
{
    std::function<double(double )> f;
    std::pair<double,double> intervals;
    double* nodes;
    double* weights;
    int n;
    std::vector<double> ai;
    public:
    void setFunction(std::function<double(double)> f)
    {
        this->f = f;
    }
    void setIntervals(std::pair<double,double> intervals)
    {
        this->intervals = intervals;
    }
    numericIntegration(int n);
    double gaussLegendre(int n, int precison);
    double hornerPolinolmial(const std::vector<double>& data, double x);
    double rectangleMethod( int n);
    double rectangleMethod( int n);
    double trapezeMethod(int precision);
    double trapezeMethod(std::function<double(double)> f,int precision);
    double simpsonMethod( int n);
    double simpsonMethod(std::function<double(double)>f, int n);

};
#endif //NUMERICINTEGRATION_H
