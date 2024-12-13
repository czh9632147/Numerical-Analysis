#include <fstream>
#include <iomanip>
#include <cstdio>
#include "Spline.hpp"

// Function to interpolate
double f(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}

// Generate evenly spaced nodes
std::vector<double> generateNodes(double a, double b, int N) {
    std::vector<double> nodes(N);
    double h = (b - a) / (N - 1);
    for (int i = 0; i < N; ++i) {
        nodes[i] = a + i * h;
    }
    return nodes;
}

int main() {
    int N = 41;
    double a = -1.0, b = 1.0;

    std::vector<double> nodes = generateNodes(a, b, N);
    std::vector<double> values(N);
    for (int i = 0; i < N; ++i) {
        values[i] = f(nodes[i]);
    }

    // Create degree-4 piecewise polynomial spline
    Degree_4_pp spline(nodes, values);

    // Create a data file for gnuplot
    std::ofstream dataFile("degree4_plot.dat");
    
    // Error variables
    double maxError = 0.0;
    double sumSquaredError = 0.0;
    int plotPoints = 500;
    double step = (b - a) / (plotPoints - 1);

    for (int i = 0; i < plotPoints; ++i) {
        double x = a + i * step;
        double exactValue = f(x);
        double splineValue = spline.evaluate(x);
        double error = std::abs(exactValue - splineValue);
        
        maxError = std::max(maxError, error);
        sumSquaredError += error * error;
        dataFile << x << " " << exactValue << " " << splineValue << "\n";
    }
    dataFile.close();

    // Create Gnuplot script
    std::ofstream scriptFile("plot_script.gp");
    scriptFile << "set terminal pngcairo enhanced size 800,600\n";
    scriptFile << "set output 'degree4_spline_plot.png'\n";
    scriptFile << "set title 'Degree-4 Piecewise Polynomial Spline Interpolation'\n";
    scriptFile << "set xlabel 'x'\n";
    scriptFile << "set ylabel 'y'\n";
    scriptFile << "set grid\n";
    scriptFile << "plot 'degree4_plot.dat' using 1:2 title 'Exact' with lines lt 1, \\\n";
    scriptFile << "     'degree4_plot.dat' using 1:3 title 'Degree-4 Spline' with lines lt 2\n";
    scriptFile.close();

    // Calculate and output errors
    double meanSquaredError = sumSquaredError / plotPoints;
    std::cout << "For N = " << N << ":\n";
    std::cout << "  Maximum Error: " << maxError << "\n";
    std::cout << "  Mean Squared Error: " << meanSquaredError << "\n";

    // Run Gnuplot
    system("gnuplot plot_script.gp");
    
    std::cout << "\nPlot has been generated: degree4_spline_plot.png\n";

    return 0;
}