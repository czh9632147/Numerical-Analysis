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

// Function to create Gnuplot script
void createGnuplotScript(int N) {
    std::ofstream scriptFile("cubic_plot_script_N" + std::to_string(N) + ".gp");
    scriptFile << "set terminal pngcairo enhanced size 800,600\n";
    scriptFile << "set output 'cubic_spline_plot_N" << N << ".png'\n";
    scriptFile << "set title 'Cubic Spline Interpolation (N = " << N << ")'\n";
    scriptFile << "set xlabel 'x'\n";
    scriptFile << "set ylabel 'y'\n";
    scriptFile << "set grid\n";
    scriptFile << "plot 'cubic_spline_plot_N" << N << ".dat' using 1:2 title 'Exact' with lines lt 1, \\\n";
    scriptFile << "     'cubic_spline_plot_N" << N << ".dat' using 1:3 title 'Spline' with lines lt 2\n";
    scriptFile.close();
}

// Function to run Gnuplot
void runGnuplot(int N) {
    std::string command = "gnuplot cubic_plot_script_N" + std::to_string(N) + ".gp";
    system(command.c_str());
}

int main() {
    std::vector<int> N_values = {6, 11, 21, 41, 81};
    double a = -1.0, b = 1.0;

    // Create a combined script for all plots
    std::ofstream combinedScript("cubic_plot_all.gp");
    combinedScript << "set terminal pngcairo enhanced size 1600,1200\n";
    combinedScript << "set output 'cubic_all_splines.png'\n";
    combinedScript << "set multiplot layout 3,2\n";
    combinedScript << "set grid\n";

    for (int N : N_values) {
        std::vector<double> nodes = generateNodes(a, b, N);
        std::vector<double> values(N);
        for (int i = 0; i < N; ++i) {
            values[i] = f(nodes[i]);
        }

        // Create cubic spline with natural boundary conditions
        CubicSpline spline(nodes, values, 1);

        // Create a data file for gnuplot
        std::ofstream dataFile("cubic_spline_plot_N" + std::to_string(N) + ".dat");
        if (!dataFile) {
            std::cerr << "Error opening data file!" << std::endl;
            return -1;
        }

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

        // Calculate mean squared error
        double meanSquaredError = sumSquaredError / plotPoints;

        // Output errors
        std::cout << "For N = " << N << ":\n";
        std::cout << "  Maximum Error: " << maxError << "\n";
        std::cout << "  Mean Squared Error: " << meanSquaredError << "\n";

        // Create individual plot script
        createGnuplotScript(N);
        
        // Add to combined plot script
        combinedScript << "set title 'N = " << N << "'\n";
        combinedScript << "plot 'cubic_spline_plot_N" << N << ".dat' using 1:2 title 'Exact' with lines lt 1, \\\n";
        combinedScript << "     'cubic_spline_plot_N" << N << ".dat' using 1:3 title 'Spline' with lines lt 2\n";

        // Generate individual plot
        runGnuplot(N);
    }

    combinedScript.close();

    // Generate combined plot
    system("gnuplot cubic_plot_all.gp");
    
    std::cout << "\nPlots have been generated:\n";
    std::cout << "- Individual plots: cubic_spline_plot_N<N>.png\n";
    std::cout << "- Combined plot: cubic_all_splines.png\n";

    return 0;
}