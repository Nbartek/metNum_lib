cmake_minimum_required(VERSION 3.31)
project(metNum_lib)

set(CMAKE_CXX_STANDARD 20)

add_library(metNum_lib STATIC library.cpp
        include/interpolation.h
        include/linearEquation.h
        include/approximation.h
        include/numericIntegration.h
        include/nonLinearEquation.h
        include/differentialEquation.h
        src/approximation.cpp
        examples/exampleNumericIntegration.cpp
        examples/exampleNonLinear.cpp)
