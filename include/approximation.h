//
// Created by nack2 on 04.06.2025.
//

#ifndef APROXIMATION_H
#define APROXIMATION_H
#include <functional>
#include <utility>

class Approximation
{
     double* nodes;
     double* weights;
    static const int MAXN = 100;
    int N = 0;

    struct fun
    {
        std::function<double(double )> f;
        std::pair<double,double> intervals;
    };
    public:
    Approximation(int n);
    double gaussLegendre(fun f, int n, int precison);
    static void partial_pivot(double A[MAXN][MAXN + 1], int n);
    void back_substitute(double A[MAXN][MAXN + 1], int n, double x[MAXN]);
    std::vector<double> gaussEliminationSolve(double mat[MAXN][MAXN + 1], int n);
    std::vector<double> approx(fun f, int degree, int gaussPoints = 3, int precision = 10);
    double evaluateApprox(const std::vector<double>& coeffs, double x);


};

#endif //APROXIMATION_H
