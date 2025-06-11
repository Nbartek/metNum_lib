#include <iostream>
#include <vector>
#include <string>
#include <cmath>   
#include <numeric> 
#include "../include/interpolation.hpp"

void test_interpolation_at_known_point() {
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> f = {0.0, 1.0, 4.0, 9.0};
    
    double xp = 2.0;
    double expected_result = 4.0;

    double actual_result = interpolate(x, f, xp);

    if (std::fabs(actual_result - expected_result) < 1e-9) {
        std::cout << "PASSED Test 1: Interpolation at a known point" << std::endl;
    } else {
        std::cout << "FAILED Test 1: Interpolation at a known point" 
                  << ". Expected: " << expected_result 
                  << ", Got: " << actual_result << std::endl;
    }
}

void test_interpolation_between_points() {
    std::vector<double> x = {1.0, 3.0};
    std::vector<double> f = {3.0, 7.0};
    
    double xp = 2.0;
    double expected_result = 5.0;
    double actual_result = interpolate(x, f, xp);

    if (std::fabs(actual_result - expected_result) < 1e-9) {
        std::cout << "PASSED Test 2: Interpolation between points" << std::endl;
    } else {
        std::cout << "FAILED Test 2: Interpolation between points" 
                  << ". Expected: " << expected_result 
                  << ", Got: " << actual_result << std::endl;
    }}

void test_interpolate(){
    test_interpolation_at_known_point();
    test_interpolation_between_points();
}