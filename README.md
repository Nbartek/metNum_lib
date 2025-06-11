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
- ## nonLinearEquation
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
