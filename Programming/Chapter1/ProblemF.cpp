#include <iostream>
#include "EquationSolver.hpp"
#include "Function.hpp"

int main() {
    // Given values
    double l = 89.0;
    double h = 49.0;
    double beta1 = 11.5;  // in degrees
    double D1 = 55.0;     // Case (a)
    double D2 = 30.0;     // Case (b)
    double D = 30.0;      // Case (c)
    
    // Newton’s Method (Part a and b)
    // Case (a) - D = 55, α ≈ 33°
    AlphaEquation alphaEq1(l, h, D1, beta1);
    Newton_Method newton1(alphaEq1, 33.0, 1e-7, 50);  // Initial guess α ≈ 33°
    double root1 = newton1.solve();
    std::cout << "Root for case (a) using Newton's Method: α ≈ " << root1 << " degrees" << std::endl;

    // Case (b) - D = 30, α with initial guess 33°
    AlphaEquation alphaEq2(l, h, D2, beta1);
    Newton_Method newton2(alphaEq2, 33.0, 1e-7, 50);  // Initial guess α = 33°
    double root2 = newton2.solve();
    std::cout << "Root for case (b) using Newton's Method: α ≈ " << root2 << " degrees" << std::endl;

    // Secant Method (Part c)
    std::cout << "case (c):" << std::endl;
    
    // Secant Method with initial guesses 10° and 20°
    AlphaEquation alphaEq(l, h, D, beta1);
    Secant_Method secant1(alphaEq, 10.0, 20.0, 1e-7, 50);
    double root_secant1 = secant1.solve();
    std::cout << "Root using Secant Method with initial guesses 10° and 20°: α ≈ " << root_secant1 << " degrees" << std::endl;
    
    // Secant Method with initial guesses 40° and 60°
    Secant_Method secant2(alphaEq, 40.0, 60.0, 1e-7, 50);
    double root_secant2 = secant2.solve();
    std::cout << "Root using Secant Method with initial guesses 40° and 60°: α ≈ " << root_secant2 << " degrees" << std::endl;

    // Try with far away initial guess (e.g. 90°)
    Secant_Method secant_far(alphaEq, 30.0, 90.0, 1e-7, 50);
    double root_secant_far = secant_far.solve();
    std::cout << "Root using Secant Method with initial guesses 30° and 90°: α ≈ " << root_secant_far << " degrees" << std::endl;

    return 0;
}

