//
// Created by nack2 on 04.06.2025.
//

#ifndef NONLINEAREQUATION_H
#define NONLINEAREQUATION_H
#include <functional>
#include <string>
using namespace std;
class nonLinearEquation
{

#define EPSILON 1E-9 /// tolerancja precyzji obliczeń
#define VAL_DIFF_EPSILON 1e-8 ///tolerancja dla uznania czy wartosci sa rowne
#define ITERATION_LIMIT 1000000 ///limit iteracji
#define MAX_STEP 0.1 ///maksymalny krok
#define DERIVATIVE_EPSILON 1E-7 ///tolerancja precyzji dla pochodnej
#define M_PI 3.14159265358979323846
    int count_newton = 0;
    int count_secant = 0;
    int count_bisection = 0;
    public:
    double metoda_bisekcji(std::function<double(double)> f, double a, double b, std::vector<double>& iteracje);
    double regulaFalsi(std::function<double(double)> f, double a, double b, std::vector<double>& iteracje);
    double metoda_siecznych(std::function<double(double)> f, double x0, double x1, std::vector<double>& iteracje);
    static void findAllRoots(std::function<double(double)> f, double start, double end, const std::string& fname);
    double metoda_newton(double (*f)(double), double (*f_prim)(double), double y0);
    double metoda_newton(double (*f)(double), double y0);

    static double  prim(double (*f)(double), double x);

    static bool testVector(double v, std::vector<double>&vect);
};

inline double nonLinearEquation::prim(double(* f)(double), double x)
{
    return (f(x + EPSILON) - f(x - EPSILON))/(2*EPSILON);
}

inline bool nonLinearEquation::testVector(double v, std::vector<double>& vect)
{
    for(double i : vect){
        if(abs(i - v) <= VAL_DIFF_EPSILON)
            return true;
    }
    return false;
}
#endif //NONLINEAREQUATION_H
