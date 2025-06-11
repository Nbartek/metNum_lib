# metNum_lib

This project is a self-made statics library for varius numerical methods and algorithms.

This readme acts as a documentation for whole project and its individual components.

Some of the comments included in the source code are polish.

---
## Organisation of code

- **examples** is a directory where you can find examples how should library functions be implemented. Each file contains code to copy to your _main.cpp_
- **include** directory has all header files for each class implemented 
- **src** directory has all source code for each and every class
- **tests** is a directory where you can find unit tests of methods included in this library

---
# What's in this library?
### This library includes 6 classes:
- ## approximation
- ## differentialEquation
- ## interpolation
- ## linearEquation
- ## nonLinearSolver
- ## numericIntegration

In this document all methods of those classes will be thoroughly described.

---
# Documentation

---
### numericIntegration

Include "metNum_lib/numericIntegration.h"

Class **numericIntegration** includes:
- Constructor with 3 parameters:
  - integer that accounts for the proper weights, and nodes set within. **IMPORTANT** this value must be between 2 and 5. This library only has those included.
  - std::pair<double,double> that holds the limits of the function.
  - std::function<double(double)> that holds the function for calculations. **IMPORTANT** even if we want to  use only methods based on polynomial this parameter must be filled.
- double gaussLegendre this metod calculates the integral using Gauss-Legendre formula. parameter decides for how many smaller parts will we divide the space inbetween limits. It returns result of integration in double.
- double hornerPolinolmial is used for calculating the value of polynomial in specified point. It needs 2 parameters:
  - std::vector<double> that holds the polynomial
  - double x that specifies the point where the result will be calculated
- 3 methods named rectangleMethodPoli, trapezeMethodPoli, simpsonMethodPoli that calculate integral from the polynomial with different formulas. Those methods need 2 parameters:
  - std::vectro<double> for the polynomial
  - int n for the 'precision', the value that the delta of intervals is divided against. The bigger, the less iterations in program = bigger possible result error.
- 3 methods named rectangleMethodFun, trapezeMethodFun, simpsonMethodFun that calculate integral from the function set up by constructor with different formulas. Those methods need 1 parameter:
    - int n for the 'precision', the value that the delta of intervals is divided against. The bigger, the less iterations in program = bigger possible result error.
**IMPORTANT** setting right intervals is crucial for proper results

### approximation

Include "metNum_lib/approximation.h"

Class **approximation** includes:
- Constructor with 1 parameter:
  - integer that accounts for the proper weights, and nodes set within. IMPORTANT this value must be between 2 and 5. This library only has those included.
- double gaussLegendre this metod calculates the integral using Gauss-Legendre formula.Parameters:
  - parameter precision decides for how many smaller parts will we divide the space inbetween limits
  - int n is **SAME AS VALUE IN CONSTRUCTOR PARAMETER**
  - std::pair<double,double> that holds the limits of the function.
  - std::function<double(double)> that holds the function for calculations.
  - It returns result of integration in double.
-  gaussEliminationSolve it's a method for solving equations with gauss elimination formula. Parameters:
  - double mat have to be a matrix of size 100x101, it is the biggest matrix that can be calculated
  - int n is used to set how big the real matrix is.
  - Return type is std::vector<double>
- std::vector<double> approx does the approximation usning gaussLegendre method. Parameters:
  - std::function for teh function approximating
  - std::pair<double,double> that holds the limits of the function.
  - int precison for setting up how gausLengerder will divide the interval
  - int gaussPoint, point where we are approximating
  - int deggre, how big the matrix is


### Interpolation

Include `metNum_lib/Interpolation.h`

The **Interpolation** class is a static utility class that provides a collection of methods for data interpolation and polynomial evaluation. Since all methods are static, you do not need to create an instance of the class; instead, you call them directly, e.g., `Interpolation::LagrangeInterpolation(...)`.

It includes:

-   `static std::pair<std::vector<double>, std::vector<double>> load_data(...)`
    This method loads a set of x and f(x) data points from a specified text file. It requires 1 parameter:
  -   `const std::string& filename`: The path to the data file.
  -   It returns a `std::pair` where the `first` element is a `std::vector<double>` of x-coordinates and the `second` is a `std::vector<double>` of the corresponding f(x) values.
  -   **IMPORTANT**: The file must be formatted with two lines, one for `x` and one for `fx`, like so:
      ```
      x: 1.0 2.0 4.0
      fx: 1.0 4.0 16.0
      ```
  -   **IMPORTANT**: This method will throw a `std::runtime_error` if the file cannot be found or opened.

-   `static double LagrangeInterpolation(...)`
    This method calculates the interpolated value at a specific point using the Lagrange polynomial formula. It requires 3 parameters:
  -   `const std::vector<double>& dataX`: A vector holding the known x-coordinates.
  -   `const std::vector<double>& dataY`: A vector holding the known f(x) values corresponding to `dataX`.
  -   `double pointX`: The point at which to calculate the interpolated value.
  -   It returns the calculated interpolated value L(pointX) as a `double`.

-   `static double newtonInterpolation(...)`
    This method calculates the interpolated value using Newton's divided differences polynomial, which often has better numerical stability. The parameters are identical to `LagrangeInterpolation`. It returns the calculated interpolated value N(pointX) as a `double`.

-   Two methods for evaluating polynomials: `normalPolynomial` and `hornerPolynomial`.
    These methods calculate the value of a polynomial for a given `x`. `hornerPolynomial` is generally more efficient and numerically stable. Both methods take 2 parameters:
  -   `const std::vector<double>& coefficients`: A vector holding the polynomial's coefficients in order of **increasing power**. For a polynomial `P(x) = 5 + 3x + 2x²`, the vector should be `{5.0, 3.0, 2.0}`.
  -   `double x`: The point at which to evaluate the polynomial.
  -   Both return the result of the polynomial P(x) as a `double`.

### NonlinearSolver

Include `metNum_lib/NonlinearSolver.h`

The **NonlinearSolver** class is a static utility class that provides a collection of robust, iterative algorithms for finding the roots of single-variable nonlinear equations (i.e., solving for `x` in `f(x) = 0`). All methods are `static` and do not require an object instance.

The methods return a `std::optional<double>`. If a root is successfully found within the given tolerance, the `optional` will contain the root's value. If the method fails (e.g., maximum iterations are reached, or a mathematical condition is not met), it will return an empty `optional` (`std::nullopt`).

It includes:

-   `static std::optional<double> bisection(...)`
    Finds a root using the Bisection method, which is guaranteed to converge if a root is bracketed. It requires at least 3 parameters:
    -   `const std::function<double(double)>& f`: The function for which to find a root.
    -   `double a`: The start of the search interval.
    -   `double b`: The end of the search interval.
    -   **IMPORTANT**: The method requires that `f(a)` and `f(b)` have opposite signs. If they do not, it will return `std::nullopt`.

-   `static std::optional<double> falsePosition(...)`
    Finds a root using the False Position (Regula Falsi) method. It generally converges faster than bisection. It requires at least 3 parameters:
    -   `const std::function<double(double)>& f`: The function.
    -   `double a`: The start of the search interval.
    -   `double b`: The end of the search interval.
    -   **IMPORTANT**: Like bisection, this method requires that `f(a)` and `f(b)` have opposite signs.

-   `static std::optional<double> secant(...)`
    Finds a root using the Secant method, which does not require the root to be bracketed but is not guaranteed to converge. It requires at least 3 parameters:
    -   `const std::function<double(double)>& f`: The function.
    -   `double x0`: The first initial guess.
    -   `double x1`: The second initial guess.

-   Two overloaded methods for `newton(...)`
    Finds a root using Newton's method, which offers very fast (quadratic) convergence but can easily diverge if the initial guess is poor or the derivative behaves unexpectedly.
    1.  `newton(const std::function<double(double)>& f, const std::function<double(double)>& f_prime, double initial_guess, ...)`
        -   This version is for when you can provide the **analytical derivative** of the function, which is faster and more accurate.
        -   `f_prime`: A function representing the exact derivative of `f`.
        -   `initial_guess`: The starting point for the search.
    2.  `newton(const std::function<double(double)>& f, double initial_guess, ...)`
        -   Use this version when you do not have the analytical derivative. It will be **calculated numerically**.
        -   `initial_guess`: The starting point for the search.
    -   **IMPORTANT**: Newton's method will fail and return `std::nullopt` if the derivative at any iteration is zero or very close to it.

### DifferentialEquation

Include `metNum_lib/differentialEquation.h`

The **DifferentialEquation** class is a static utility library for numerically solving first-order Ordinary Differential Equations (ODEs) of the form `y' = f(t, y)`. All methods are `static` and do not require creating an instance of the class.

All solver methods return a `std::vector<std::pair<double, double>>`, which represents the discrete points of the solution path. Each pair in the vector is a `(t, y)` coordinate.

-   `static std::vector<Point> euler(...)`
    Solves an ODE using the simple and fast, but less accurate, explicit Euler method. It requires 5 parameters:
    -   `const std::function<double(double, double)>& derivative`: The function `f(t, y)` that defines the ODE.
    -   `double t_start`: The initial time (or independent variable value).
    -   `double t_end`: The final time to solve until.
    -   `double y0`: The initial value of the dependent variable, `y(t_start)`.
    -   `int N`: The number of discrete steps to take across the interval.

-   `static std::vector<Point> heun(...)`
    Solves an ODE using Heun's method, a second-order predictor-corrector method that offers better accuracy than Euler's method. The parameters are identical to the `euler` method.

-   `static std::vector<Point> midpoint(...)`
    Solves an ODE using the second-order midpoint method. Its accuracy is comparable to Heun's method. The parameters are identical to the `euler` method.

-   `static std::vector<Point> rungeKutta4(...)`
    Solves an ODE using the classical fourth-order Runge-Kutta (RK4) method. This method is highly accurate and widely used for general-purpose ODE solving. The parameters are identical to the `euler` method.

**IMPORTANT**: The accuracy of all methods depends on the number of steps, `N`. A larger `N` results in a smaller step size `h = (t_end - t_start) / N`, which generally leads to a more accurate solution at the cost of more computation.

### LinearEquation

Include `metNum_lib\linearEquation.h`

The **LinearEquation** class is a static utility library for solving systems of linear equations of the form `Ax = b`. It uses modern C++ data structures and is designed to be safe and flexible. All methods are `static`.

For convenience, the following type aliases are defined within the class:
- `using Matrix = std::vector<std::vector<double>>;`
- `using Vector = std::vector<double>;`

The library includes:

-   `static std::pair<Matrix, Vector> loadSystemFromFile(...)`
    Loads a linear system from a text file. It requires 1 parameter:
    -   `const std::string& filename`: The path to the data file.
    -   It returns a `std::pair` containing the `Matrix A` and the `Vector b`.
    -   **IMPORTANT**: The file format is strict and must be as follows:
        ```
        size: 3
        b: 6 -4 27
        A_row0: 1 1 1
        A_row1: 0 2 5
        A_row2: 2 5 -1
        ```
    -   **IMPORTANT**: This method throws a `std::runtime_error` if the file cannot be opened or if its contents do not match the required format.

-   `static std::optional<Vector> solveWithGaussianElimination(...)`
    Solves the system `Ax = b` using the Gaussian elimination algorithm with partial pivoting for numerical stability. It requires 2 parameters:
    -   `Matrix A`: The N x N coefficient matrix.
    -   `Vector b`: The N-element constant vector.
    -   It returns a `std::optional<Vector>`. If a unique solution `x` is found, the optional will contain the solution vector. If the matrix `A` is singular (i.e., no unique solution exists), it will return `std::nullopt`.

-   `static std::optional<Vector> solveWithLU(...)`
    Solves the system `Ax = b` by first performing an LU decomposition with partial pivoting. This method is efficient if you need to solve for multiple `b` vectors with the same `A` matrix (though this implementation does not expose the L and U factors directly). The parameters are identical to the Gaussian elimination method.
    -   It also returns a `std::optional<Vector>`, returning `std::nullopt` if the matrix `A` is singular.