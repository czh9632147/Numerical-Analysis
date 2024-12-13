#include <fstream>
#include <iomanip>
#include <cstdio>
#include "Spline.hpp"

// Test function: sin(M_PI * x)
double testFunction(double x) {
    return std::sin(M_PI * x);
}

// Function to generate gnuplot script and plot
void generateAndPlotGnuplot(const std::string& type, const std::vector<double>& knots, const std::vector<double>& values, const CubicSpline& spline) {
    std::ofstream gnuplot_script;
    gnuplot_script.open("plot_" + type + ".gp"); // Open the gnuplot script file

    if (!gnuplot_script.is_open()) {
        std::cerr << "Failed to open gnuplot script file." << std::endl;
        return;
    }

    // Set up gnuplot script
    gnuplot_script << "set terminal png size 800,600\n";
    gnuplot_script << "set output '" << type << "_spline_plot.png'\n";
    gnuplot_script << "set grid\n";
    gnuplot_script << "set title '" << type << " Cubic Spline'\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'y'\n";
    gnuplot_script << "set xr [-1.1:1.1]\n";
    gnuplot_script << "set yr [-1.1:1.1]\n";
    gnuplot_script << "plot '-' with lines title 'Spline', \\" << std::endl;
    gnuplot_script << "     '-' with points title 'Knots'\n";

    // Plot the spline curve
    for (double x = -1.0; x <= 1.0; x += 0.01) {
        gnuplot_script << x << " " << spline.evaluate(x) << std::endl;
    }
    gnuplot_script << "e\n";

    // Plot the knots points
    for (size_t i = 0; i < knots.size(); ++i) {
        gnuplot_script << knots[i] << " " << values[i] << std::endl;
    }
    gnuplot_script << "e\n";

    gnuplot_script.close(); // Close the gnuplot script file

    // Execute the gnuplot script
    std::string command = "gnuplot plot_" + type + ".gp";
    int ret = system(command.c_str());
    if (ret != 0) {
        std::cerr << "Failed to execute gnuplot script." << std::endl;
    }
}

int main() {
    std::vector<double> knots = {-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1};
    std::vector<double> values(knots.size());

    // Calculate function values at knots
    for (size_t i = 0; i < knots.size(); ++i) {
        values[i] = testFunction(knots[i]);
    }

    // Boundary conditions for clamped spline
    double leftDerivative = -M_PI;
    double rightDerivative = -M_PI;

    // Natural boundary condition
    std::cout << "Natural Spline:" << std::endl;
    CubicSpline naturalSpline(knots, values, 1);
    generateAndPlotGnuplot("Natural", knots, values, naturalSpline);

    // Clamped boundary condition
    std::cout << "Clamped Spline:" << std::endl;
    CubicSpline clampedSpline(knots, values, 2, leftDerivative, rightDerivative);
    generateAndPlotGnuplot("Clamped", knots, values, clampedSpline);

    // Periodic boundary condition
    std::cout << "Periodic Spline:" << std::endl;
    CubicSpline periodicSpline(knots, values, 3);
    generateAndPlotGnuplot("Periodic", knots, values, periodicSpline);

    return 0;
}