//
// Created by nack2 on 04.06.2025.
//
#include "../include/numericIntegration.h"

#include <iostream>

numericIntegration::numericIntegration(int n,std::pair<double,double> intervals,std::function<double(double)>f)
{
    this->n = n;
    this->intervals = intervals;
    this->f = f;
    static  double nodes2[] = { -1/sqrt(3), 1/sqrt(3) };
    static  double weights2[] = { 1.0, 1.0 };

    static  double nodes3[] = { -0.7745966692, 0.0, 0.7745966692 };
    static  double weights3[] = { 0.5555555556, 0.8888888889, 0.5555555556 };

    static  double nodes4[] = {-0.8611363116, -0.3399810436,0.3399810436,  0.8611363116};
    static  double weights4[] = {0.3478548451, 0.6521451549,0.6521451549, 0.3478548451};

    static  double nodes5[] = {-0.9061798459, -0.5384693101, 0.0,0.5384693101,  0.9061798459};
    static  double weights5[] = {0.2369268850, 0.4786286705, 0.5688888889,0.4786286705, 0.2369268850};

    switch (n) {
    case 2:
        nodes = nodes2;
        weights = weights2;
        break;
    case 3:
        nodes = nodes3;
        weights = weights3;
        break;
    case 4:
        nodes = nodes4;
        weights = weights4;
        break;
    case 5:
        nodes = nodes5;
        weights = weights5;
        break;
    default:
        throw std::invalid_argument("od 2 do 5");
    }
}

double numericIntegration::gaussLegendre(int precison)
{
    double result = 0.0;
    double h = (intervals.second - intervals.first) / precison;
    for (int j = 0; j < precison; ++j)
    {
        double aj = intervals.first+j*h;
        double bj = aj+h;
        for (int i = 0; i < n; ++i) {
            double xi = ((bj - aj) / 2.0) * nodes[i] + (aj + bj) / 2.0;
            result += weights[i] * f(xi);
        }
    }
    result *= h / 2.0;
    return result;
}

double numericIntegration::hornerPolinolmial(const std::vector<double>& data, double x)
{
    double wynik = data.back();
    for (int i = data.size() - 2; i >= 0; i--) {
        wynik = wynik * x + data[i];
    }
    return wynik;
}

double numericIntegration::rectangleMethodPoli(std::vector<double>ai,int n)
{
    double h = (intervals.second - intervals.first) / n;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double mid = intervals.first + i * h + h / 2.0;
        sum += hornerPolinolmial(ai, mid);
    }

    return sum * h;
}

double numericIntegration::rectangleMethodFun(int n)
{
    double h = (intervals.second - intervals.first) / n;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double mid = intervals.first + i * h + h / 2.0;
        sum += f( mid);
    }

    return sum * h;
}

double numericIntegration::trapezeMethodPoli(std::vector<double> ai,int n)
{
    // liczymy h
    double h = (intervals.second - intervals.first) / n;

    // ONa krancach przedziałów
    double suma = (hornerPolinolmial(ai, intervals.first) + hornerPolinolmial(ai, intervals.second)) / 2.0;

    // wewnątrz
    for (int i = 1; i < n; i++) {
        //liczymy xi
        double xi = intervals.first + i * h;
        suma += hornerPolinolmial(ai, xi);
    }

    // suma pól trapezów jako wynik
    return suma * h;
}

double numericIntegration::trapezeMethodFun( int n)
{
    double h = (intervals.second - intervals.first) / n;

    // ONa krancach przedziałów
    double suma = (f(intervals.first)+ f( intervals.second)) / 2.0;

    // wewnątrz
    for (int i = 1; i < n; i++) {
        //liczymy xi
        double xi = intervals.first + i * h;
        suma += f(xi);
    }

    // suma pól trapezów jako wynik
    return suma * h;
}

double numericIntegration::simpsonMethodPoli(std::vector<double> ai,int n)
{
    if (n % 2 != 0) {
        std::cerr << "Liczba kroków (n) musi być parzysta dla metody Simpsona!" << std::endl;
        return 0.0;
    }

    double h = (intervals.second - intervals.first) / n;
    double suma = hornerPolinolmial(ai, intervals.first) + hornerPolinolmial(ai, intervals.second);

    for (int i = 1; i < n; i++) {
        double xi = intervals.first + i * h;
        double fx = hornerPolinolmial(ai, xi);

        if (i % 2 == 0) {
            suma += 2 * fx;
        } else {
            suma += 4 * fx;
        }
    }

    return (h / 3.0) * suma;
}

double numericIntegration::simpsonMethodFun( int n)
{
    if (n % 2 != 0) {
        std::cerr << "Liczba kroków (n) musi być parzysta dla metody Simpsona!" << std::endl;
        return 0.0;
    }

    double h = (intervals.second - intervals.first) / n;
    double suma = f(intervals.first) + f( intervals.second);

    for (int i = 1; i < n; i++) {
        double xi = intervals.first + i * h;
        double fx = f(xi);

        if (i % 2 == 0) {
            suma += 2 * fx;
        } else {
            suma += 4 * fx;
        }
    }

    return (h / 3.0) * suma;
}

