#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "Interpolation.hpp"

// Runge function f(x) = 1/(1 + 25x^2)
double rungeFunction(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}

// Function to generate x values based on Chebyshev nodes
std::vector<double> generateXValues(int n) {
    std::vector<double> x_values = ChebyshevInterpolation::chebyshevNodes(n, -1, 1);
    return x_values;
}
    

// Generate y values for given x values using the Runge function
std::vector<double> generateYValues(const std::vector<double>& x_values) {
    std::vector<double> y_values;
    for (double x : x_values) {
        y_values.push_back(rungeFunction(x));
    }
    return y_values;
}

int main() {
    std::vector<int> n_values = {5, 10, 15, 20};

    for (int n : n_values) {

        std::cout << "\nNewton interpolation at Chebyshev points for n = " << n << ":\n";

        std::vector<double> x_values = generateXValues(n);
        std::vector<double> y_values = generateYValues(x_values);

        ChebyshevInterpolation chebyshevInterpolation(x_values, y_values);
        std::cout << ": " << "P" << n << "(x) = " 
            << chebyshevInterpolation.getExpression() << std::endl;
        
        std::cout << "----------------------------------" << std::endl;
    }

    return 0;
}
