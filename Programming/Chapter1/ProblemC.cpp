#include <iostream>
#include "EquationSolver.hpp"
#include "Function.hpp"
#include <cmath>

int main() {
    TanEquation tanEq;

    // Find the root near 4.5
    Newton_Method solver_4_5(tanEq, 4.5);
    double root_4_5 = solver_4_5.solve();
    std::cout << "Root near 4.5: " << root_4_5 << std::endl;

    // Find the root near 7.7
    Newton_Method solver_7_7(tanEq, 7.7);
    double root_7_7 = solver_7_7.solve();
    std::cout << "Root near 7.7: " << root_7_7 << std::endl;

    return 0;
}

