#include "test_approx.hpp"
#include "test_diff_eq.hpp"
#include "test_nonlin_eq.hpp"
#include "test_lin_eq.hpp"
#include "test_interpolate.hpp"
#include "test_integration.hpp"

#include "../examples/example_interpolate.hpp"
#include "../examples/example_solve_diff_eq.hpp"

int main(){
    test_solvers();

    std::cout << "\n======================" << std::endl;
    std::cout << "All tests passed!" << std::endl;
    std::cout << "======================" << std::endl;
    test_bisection();
    test_false_position();
    test_secant();
    test_newton();

    std::cout << "\n======================" << std::endl;
    std::cout << "All tests passed!" << std::endl;
    std::cout << "======================" << std::endl;

    test_integration();

    test_file_loading();
    test_solvers();

    std::cout << "\n======================" << std::endl;
    std::cout << "All tests passed!" << std::endl;
    std::cout << "======================" << std::endl;
    test_constructor_and_gauss();
    test_newton_cotes_methods();
    test_horner_method();

    std::cout << "\n======================" << std::endl;
    std::cout << "All tests passed!" << std::endl;
    std::cout << "======================" << std::endl;

    return 0;
}