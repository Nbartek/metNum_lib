//
// Created by nack2 on 04.06.2025.
//
#include "../include/differentialEquation.h"

double DifferentialEquation::T_s(double t, double y0)
{
    double alpha = -0.9e-12;
    double T0 = 0;
    return 1 / std::pow(3* alpha * t + T0,1.0/3.0);
}

void DifferentialEquation::heunMethod(fun f, double y0, int N)
{
    std::cout<<"Metoda Heuna\n";
    double h = (f.interval.second - f.interval.first) / N;
    double x = f.interval.first, y = y0;
    double error = 0;
    //std::cout <<0<< "\t " << x << "\t" << y <<"\t"<<0<< std::endl;
    for (int i = 0; i < N; ++i) {
        double k1 = f.derivative(x, y);
        double k2 = f.derivative(x + h, y + h * k1);
        y = y + (h / 2) * (k1 + k2);
        x = x + h;
        //double solution =T_s(x,y);
        std::cout <<i+1<< "\t " << x << "\t" << y <<"\t"<< std::endl;
        //error +=fabs(solution - y);
    }
}

void DifferentialEquation::rungeKutta4Method(fun f, double y0, int N)
{
    std::cout<<"Metoda Runge-Kutty rzedu 4\n";
    double h = (f.interval.second - f.interval.first) / N;
    double x = f.interval.first, y = y0;
    double error = 0;
    //std::cout <<0<< "\t " << x << "\t" << y <<"\t"<<0<< std::endl;
    for (int i = 0; i < N; ++i) {
        double k1 = h * f.derivative(x, y);
        double k2 = h * f.derivative(x + h / 2, y + k1 / 2);
        double k3 = h * f.derivative(x + h / 2, y + k2 / 2);
        double k4 = h * f.derivative(x + h, y + k3);
        y = y + (1.0 / 6) * (k1 + 2*k2 + 2*k3 + k4);
        x = x + h;
        //double solution =T_s(x,y);
        std::cout <<i+1<< "\t " << x << "\t" << y <<"\t"<< std::endl;
        //error +=fabs(solution - y);
    }
    //error +=fabs(solution - y);
    //std::cout<< x << "\t" << y <<"\t"<<error<< std::endl;
    //std::cout<<"Sredni blad kwadratowy: "<<error/N<<std::endl;
}

void DifferentialEquation::midpcointMethod(fun f, double y0, int N)
{
    std::cout<<"Metoda punktu srodkoweg\n";
    double h = (f.interval.second - f.interval.first) / N;
    double x = f.interval.first, y = y0;
    //std::cout <<0<< "\t " << x << "\t" << y <<"\t"<<0<< std::endl;
    double error = 0;
    for (int i = 0; i < N; ++i) {
        double k1 = f.derivative(x, y);
        double k2 = f.derivative(x + h / 2, y + h / 2 * k1);
        y = y + h * k2;
        x = x + h;
        //double solution =T_s(x,y);
        std::cout <<i+1<< "\t " << x << "\t" << y <<"\t"<< std::endl;
        //error +=fabs(solution - y);
    }
    // error +=fabs(solution - y);
    // std::cout<< x << "\t" << y <<"\t"<<error<< std::endl;
    //std::cout<<"Sredni blad kwadratowy: "<<error/N<<std::endl;
}

double DifferentialEquation::eulerMethod(fun f, double y0, int N)
{
    double h = (f.interval.second - f.interval.first) / N;
    std::vector<double> x_vals(N + 1);
    std::vector<double> y_vals(N + 1);

    x_vals[0] = f.interval.first;
    y_vals[0] = y0;


    for (int i = 0; i < N; ++i) {
        x_vals[i + 1] = x_vals[i] + h;
        y_vals[i + 1] = y_vals[i] + h * f.derivative(x_vals[i], y_vals[i]);
    }

    std::cout << std::fixed << std::setprecision(6);
    //std::cout << "N = " << N << "\n";
    double error = 0;
    error +=fabs(solution - y_vals[N]);
    std::cout<< x_vals[N] << "\t" << y_vals[N] <<"\t"<<error<< std::endl;

    //std::cout << "Temperatura po " << f.interval.second << "s: " << y_vals[N] << " K\n\n";
    for (int i = 0; i <= N; ++i) {
        double solution =T_s(x_vals[i],y_vals[i]);
        // std::cout << i << "\t" << x_vals[i] << "\t" << y_vals[i]<<"\t"<<solution;
        //std::cout << "\n";
    }
    return y_vals[N];
}
