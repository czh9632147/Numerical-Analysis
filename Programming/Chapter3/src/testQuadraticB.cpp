#include <vector>
#include <fstream>
#include <functional>
#include "Spline.hpp"

double f(double x) {
    return 1.0 / (1.0 + x * x);
}

int main() {
    // Define the points for the spline
    std::vector<double> prepoints = {-6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
    std::vector<double> points = {-5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
    std::vector<double> values(points.size() + 1);
    for (size_t i = 0; i < points.size() + 1; ++i) {
        values[i] = f((prepoints[i] + prepoints[i + 1]) / 2);
    }
    

    // Create a QuadraticBSpline object using the specific function
    QuadraticBSpline spline(points, values);

    // Generate a gnuplot script to plot the spline and the original function
    std::ofstream gnuplot_script("plot.gp");
    gnuplot_script << "set terminal png size 800,600\n";
    gnuplot_script << "set output 'QBspline_plot.png'\n";
    gnuplot_script << "set grid\n";
    gnuplot_script << "set xrange [-5:5]\n";
    gnuplot_script << "set yrange [0:1.1]\n"; // Adjusted yrange to accommodate the function values
    gnuplot_script << "plot '-' with lines title 'Original Function', '-' with lines title 'Quadratic B-Spline'\n";

    // Plot the original function
    for (double x = -5; x <= 5; x += 0.1) {
        gnuplot_script << x << " " << f(x) << "\n";
    }
    gnuplot_script << "e\n";  // End of the first plot

    // Plot the quadratic B-spline
    for (double x = -5; x <= 5; x += 0.1) {
        gnuplot_script << x << " " << spline.eval(x) << "\n";
    }
    gnuplot_script << "e\n";  // End of the second plot

    gnuplot_script.close();

    // Run the Gnuplot script
    system("gnuplot plot.gp");

    return 0;
}