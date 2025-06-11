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
    numericIntegration(int n, std::pair<double,double> intervals,std::function<double(double)>f);
    double gaussLegendre(int precison);
    double hornerPolinolmial(const std::vector<double>& data, double x);
    double rectangleMethodPoli(std::vector<double>, int n);
    double rectangleMethodFun( int n);
    double trapezeMethodPoli(std::vector<double>,int n);
    double trapezeMethodFun(int n);
    double simpsonMethodPoli(std::vector<double>, int n);
    double simpsonMethodFun( int n);

};
#endif //NUMERICINTEGRATION_H
