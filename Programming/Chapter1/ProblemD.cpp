#include <iostream>
#include "EquationSolver.hpp"
#include "Function.hpp"

int main() {
    // Test Case 1: sin(x/2) - 1 with x0 = 0, x1 = pi/2
    SinEquation sinEq1;
    Secant_Method secant1_1(sinEq1, 0, M_PI / 2);
    double root1_1 = secant1_1.solve();
    std::cout << "Root for sin(x/2) - 1: " << root1_1 << std::endl;
    
    // Test Case 1: sin(x/2) - 1 with x0 = 0, x1 = 0.5
    SinEquation sinEq2;
    Secant_Method secant1_2(sinEq2, 0, 0.5);
    double root1_2 = secant1_2.solve();
    std::cout << "Root for sin(x/2) - 1: " << root1_2 << std::endl;

    // Test Case 2: e^x - tan(x) with x0 = 1, x1 = 1.4
    ExpTanEquation expTanEq1;
    Secant_Method secant2_1(expTanEq1, 1.0, 1.4);
    double root2_1 = secant2_1.solve();
    std::cout << "Root for e^x - tan(x): " << root2_1 << std::endl;
    
    // Test Case 2: e^x - tan(x) with x0 = 0.8, x1 = 2
    ExpTanEquation expTanEq2;
    Secant_Method secant2_2(expTanEq2, 0.8, 2.0);
    double root2_2 = secant2_2.solve();
    std::cout << "Root for e^x - tan(x): " << root2_2 << std::endl;

    // Test Case 3: x^3 - 12x^2 + 3x + 1 with x0 = 0, x1 = -0.5
    PolyEquation polyEq1;
    Secant_Method secant3_1(polyEq1, 0, -0.5);
    double root3_1 = secant3_1.solve();
    std::cout << "Root for x^3 - 12x^2 + 3x + 1: " << root3_1 << std::endl;

    // Test Case 3: x^3 - 12x^2 + 3x + 1 with x0 = 0.5, x1 = _1
    PolyEquation polyEq2;
    Secant_Method secant3_2(polyEq2, 0.5, -1.0);
    double root3_2 = secant3_2.solve();
    std::cout << "Root for x^3 - 12x^2 + 3x + 1: " << root3_2 << std::endl;
    
    return 0;
}
