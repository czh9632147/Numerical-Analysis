#include <cmath>
#include <fstream>
#include "Spline.hpp"

// Define a helper function to calculate the value of f(x) = sin(pi * x)
double f(double x) {
    return sin(M_PI * x);
}

int main() {
    // Define the knots and corresponding function values
    std::vector<double> knots = {-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1}; // Knots within the interval [-1, 1]
    std::vector<double> values(knots.size());
    for (size_t i = 0; i < knots.size(); ++i) {
        values[i] = f(knots[i]);
    }

    // Create a LinearSpline object
    LinearSpline spline(knots, values);

    // Prepare the Gnuplot script
    std::ofstream gnuplot_script("plot.gp");
    gnuplot_script << "set terminal png size 800,600\n";
    gnuplot_script << "set output 'Lppspline_plot.png'\n";
    gnuplot_script << "set grid\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'f(x)'\n";
    gnuplot_script << "set title 'Linear Spline'\n";
    gnuplot_script << "plot '-' with lines title 'Linear Spline'\n";

    // Write the data points to the Gnuplot script
    for (size_t i = 0; i < knots.size(); ++i) {
        gnuplot_script << knots[i] << " " << values[i] << "\n";
    }

    // Calculate and write the linear spline interpolation results
    double x = -1.0;
    while (x <= 1.0) {
        gnuplot_script << x << " " << spline.evaluate(x) << "\n";
        x += 0.01; // The step size can be adjusted as needed
    }

    gnuplot_script.close();

    // Run the Gnuplot script
    system("gnuplot plot.gp");

    return 0;
}