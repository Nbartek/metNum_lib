//
// Created by nack2 on 04.06.2025.
//

#ifndef DIFFERENTIALEQUATION_H
#define DIFFERENTIALEQUATION_H
#include <functional>
#include <utility>

class DifferentialEquation
{
    struct fun {
        std::function<double(double, double)> derivative;
        std::pair<double, double> interval;
    };
    public:
    double T_s(double t,double y0);
    void heunMethod(fun f, double y0, int N);
    void rungeKutta4Method(fun f, double y0, int N);
    void midpcointMethod(fun f, double y0, int N);
    double eulerMethod(fun f, double y0, int N);


};
#endif //DIFFERENTIALEQUATION_H
