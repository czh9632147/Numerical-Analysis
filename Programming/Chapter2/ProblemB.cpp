#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "Interpolation.hpp"

// Function to generate x values based on formula: x_i = -5 + 10*i/n
std::vector<double> generateXValues(int n) {
    std::vector<double> x_values;
    for (int i = 0; i <= n; ++i) {
        x_values.push_back(-5.0 + 10.0 * i / n);
    }
    return x_values;
}

// The Runge function f(x) = 1/(1 + x^2)
double rungeFunction(double x) {
    return 1.0 / (1.0 + x * x);
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
    std::vector<int> n_values = {2, 4, 6, 8};

    for (int n : n_values) {
        std::cout << "\nNewton interpolation for n = " << n << ":\n";

        std::vector<double> x_values = generateXValues(n);
        std::vector<double> y_values = generateYValues(x_values);

        NewtonInterpolation chebyshevInterpolation(x_values, y_values);
        std::cout << "n = " << n << ": " << chebyshevInterpolation.getExpression() << std::endl;

        std::cout << "----------------------------------" << std::endl;
    }

    return 0;
}
