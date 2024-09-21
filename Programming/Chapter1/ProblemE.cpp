#include <iostream>
#include "EquationSolver.hpp"
#include "Function.hpp"

int main() {
    double L = 10.0;
    double r = 1.0;
    double V = 12.4;

    // Define the function for the trough volume
    TroughVolume volumeEq(L, r, V);

    // Bisection Method
    Bisection_Method bisection(volumeEq, 0.0, 1.0, 1e-7, 1e-6, 50);
    double root_bisection = bisection.solve();
    std::cout << "Root using Bisection Method: " << root_bisection << " ft" << std::endl;

    // Newton's Method
    Newton_Method newton(volumeEq, 0.5, 1e-7, 50); // Initial guess h = 0.5
    double root_newton = newton.solve();
    std::cout << "Root using Newton's Method: " << root_newton << " ft" << std::endl;

    // Secant Method
    Secant_Method secant(volumeEq, 0.0, 0.5, 1e-7, 50); // Initial guesses h0 = 0, h1 = 0.5
    double root_secant = secant.solve();
    std::cout << "Root using Secant Method: " << root_secant << " ft" << std::endl;

    return 0;
}

