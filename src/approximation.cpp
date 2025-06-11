//
// Created by nack2 on 04.06.2025.
//
#include "../include/approximation.h"

#include <iostream>
#include <stdexcept>

Approximation::Approximation(int n)
{
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

double Approximation::gaussLegendre(fun f, int n, int precison)
{

    double result = 0.0;
    double h = (f.intervals.second - f.intervals.first) / precison;
    for (int j = 0; j < precison; ++j)
    {
        double aj = f.intervals.first+j*h;
        double bj = aj+h;
        for (int i = 0; i < n; ++i) {
            double xi = ((bj - aj) / 2.0) * nodes[i] + (aj + bj) / 2.0;
            result += weights[i] * f.f(xi);
        }
    }
    result *= h / 2.0;
    return result;
}

void Approximation::partial_pivot(double A[100][101], int n)
{
    for (int i = 0; i < n; i++) {
        int pivot_row = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[pivot_row][i])) {
                pivot_row = j;
            }
        }
        if (pivot_row != i) {
            for (int j = i; j <= n; j++) {
                std::swap(A[i][j], A[pivot_row][j]);
            }
        }
        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k <= n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }
}

void Approximation::back_substitute(double A[100][101], int n, double x[100])
{
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
    }
}

std::vector<double> Approximation::gaussEliminationSolve(double mat[100][101], int n)
{
    double A[MAXN][MAXN + 1];
    // kopiujemy oryginalną macierz, żeby jej nie nadpisać
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= n; ++j)
            A[i][j] = mat[i][j];

    double x[MAXN];
    partial_pivot(A, n);
    back_substitute(A, n, x);

    std::vector<double> result;
    for (int i = 0; i < n; ++i) {
        if (std::isnan(x[i])) {
            std::cerr << "Równanie sprzeczne lub układ osobliwy.\n";
            return {};
        }
        result.push_back(x[i]);
    }
    return result;
}

std::vector<double> Approximation::approx(fun f, int degree, int gaussPoints, int precision)
{
    N = degree + 1;
    double A[MAXN][MAXN + 1] = {0};

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fun integrand;
            integrand.f = [i, j](double x) {
                return std::pow(x, i) * std::pow(x, j);
            };
            integrand.intervals = f.intervals;
            A[i][j] = gaussLegendre(integrand, gaussPoints, precision);
        }

        fun b_integrand;
        b_integrand.f = [i, &f](double x) {
            return f.f(x) * std::pow(x, i);
        };
        b_integrand.intervals = f.intervals;
        A[i][N] = gaussLegendre(b_integrand, gaussPoints, precision);
    }

    return gaussEliminationSolve(A, N);
}

double Approximation::evaluateApprox(const std::vector<double>& coeffs, double x)
{
    double sum = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        sum += coeffs[i] * std::pow(x, i);
    }
    return sum;
}
