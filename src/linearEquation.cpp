//
// Created by nack2 on 04.06.2025.
//
#include "../include/linearEquation.h"

void LinearEquation::print(double mat[100][101])
{
    for (int i = 0; i < N; i++, printf("\n"))
        for (int j = 0; j < N; j++)
            printf("%lf ", mat[i][j]);

    printf("\n");
}

void LinearEquation::load_data_matrix(const std::string& filename, double result[100][101])
{
    std::ifstream inputFile(filename);
    std::string line;

    if (!inputFile) {
        std::cerr << "Error: Unable to open file " << filename << "\n";

    }
    std::cout << filename << " is opened\n";

    getline(inputFile, line);
    std::istringstream ss1(line.substr(4));
    ss1 >> N;                           //wgranie wielkosci macierzy

    getline(inputFile, line);

    std::istringstream ss(line.substr(4));

    double value;
    for (int i = 0; i < N; i++)
    {
        ss >> value;
        result[i][N]=value;
    }

    for (int i = 0; i < N; i++)
    {
        getline(inputFile, line);
        std::istringstream ss(line.substr(6));
        for (int j = 0; j < N; j++)
        {
            ss >> value;
            result[i][j] = value;
        }
    }
}

void LinearEquation::partial_pivot(double A[100][101], int n)
{
    for (int i = 0; i < n; i++) {
        print(A);
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

void LinearEquation::back_substitute(double A[100][101], int n, double x[100])
{
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
    }
}

void LinearEquation::backwardSubstitution(double U[100][100], double y[100], double x[100])
{
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i]/=U[i][i];
    }
}

void LinearEquation::forwardSubstitution(double L[100][100], double b[100], double y[100])
{
    for (int i = 0; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
        //gdyby jedynki były gdzie indziej
        //y[i] /= L[i][i];
    }
}

void LinearEquation::test(double A[100][100], double x[100])
{
    double result[NMAX] ={0};
    for (int i = 0; i < N; i++) {
        result[i] = 0;
        for (int j = 0; j < N; j++) {
            result[i] += A[i][j] * x[j];
        }
    }
    std::cout << "Wyswietlanie b obliczonego:\n";
    for (int i = 0; i < N; i++) {
        std::cout << "b[" << i << "] = " << result[i] << std::endl;
    }
}

void LinearEquation::luDecomposition(double A[100][100], double L[100][100], double U[100][100])
{
    for (int i = 0; i < N; i++)
    {
        P[i] = i;
    }
    double A_copy[NMAX][NMAX];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A_copy[i][j] = A[i][j];

    // Pivotowanie + rozkład
    for (int i = 0; i < N; i++) {
        // Znajdź wiersz z największym elementem w kolumnie i (pivoting)
        int pivot = i;
        for (int j = i + 1; j < N; j++) {
            if (fabs(A_copy[j][i]) > fabs(A_copy[pivot][i])) {
                pivot = j;
            }
        }

        // Zamiana wierszy w A_copy
        if (pivot != i) {
            std::swap(P[i],P[pivot]);
            for (int k = 0; k < N; k++) {
                std::swap(A_copy[i][k], A_copy[pivot][k]);
                std::swap(L[i][k], L[pivot][k]);  // zamień też odpowiednie wiersze w L (częściowo zapełnione)

            }
        }

        // U[i][j]
        for (int j = i; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += L[i][k] * U[k][j];
            U[i][j] = A_copy[i][j] - sum;
        }

        // L[j][i]
        for (int j = i; j < N; j++) {
            if (i == j)
                L[i][i] = 1.0;  // diagonala L
            else {
                double sum = 0;
                for (int k = 0; k < i; k++)
                    sum += L[j][k] * U[k][i];
                L[j][i] = (A_copy[j][i] - sum) / U[i][i];
            }
        }
    }
}

void LinearEquation::gaussElimination(double mat[100][101])
{
    double x[MAXN];

    partial_pivot(mat, N);
    back_substitute(mat, N, x);
    for (size_t i = 0; i < N; i++)
    {
        if (isnan(x[i])) {
            std::cout << "Rownanie sprzeczne\n";
            return;
        }
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
}

void LinearEquation::luMethod(std::string fileName, bool zShow, bool doTest, bool noShowMatrix)
{
    double A[NMAX][NMAX];
    double b[NMAX];
    load_data_matrix(fileName,A,b);
    if (noShowMatrix)print(A);
    double L[NMAX][NMAX] = {0}, U[NMAX][NMAX] = {0};
    double z[NMAX], x[NMAX];
    luDecomposition(A, L, U);
    if (noShowMatrix)
    {
        print(L);
        print(U);
    }
    double b_permuted[N];
    for (int i = 0; i < N; i++)
        b_permuted[i] = b[P[i]];
    //dla macierzy trójkątniej dolnej
    forwardSubstitution(L, b_permuted, z);
    //dla macierzy trójkątnej górnej
    backwardSubstitution(U, z, x);

    std::cout << "Rozwiazanie dla x:\n";
    for (int i = 0; i < N; i++) {
        if (isnan(x[i]))
        {
            std::cout<<"Rownanie sprzeczne\n";
            return;
        }
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    if (zShow)
    {
        std::cout << "Rozwiazanie dla z:\n";
        for (int i = 0; i < N; i++) {
            std::cout << "z[" << i << "] = " << z[i] << std::endl;
        }
    }
    if (doTest) test(A,x);
}
