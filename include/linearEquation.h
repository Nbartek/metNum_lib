//
// Created by nack2 on 04.06.2025.
//

#ifndef LINEAREQUATION_H
#define LINEAREQUATION_H
#include <string>

class LinearEquation
{
    int N = NAN;
    static const int NMAX = 100;
    static int P[NMAX];
    public:
    void print(double mat[NMAX][NMAX + 1]);
    static void load_data_matrix(const std::string& filename,double result[NMAX][NMAX+1] );
    void partial_pivot(double A[NMAX][NMAX + 1], int n);
    void back_substitute(double A[NMAX][NMAX + 1], int n, double x[NMAX]);
    void backwardSubstitution(double U[NMAX][NMAX], double y[NMAX], double x[NMAX]);
    void forwardSubstitution(double L[NMAX][NMAX], double b[NMAX], double y[NMAX]);
    void test(double A[NMAX][NMAX], double x[NMAX]);
    void luDecomposition(double A[NMAX][NMAX], double L[NMAX][NMAX], double U[NMAX][NMAX]);
    void gaussElimination(double mat[NMAX][NMAX+1]);
    void luMethod(std::string fileName,bool zShow, bool doTest, bool noShowMatrix);
};
#endif //LINEAREQUATION_H
