#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "Spline.hpp"

// Function to calculate the derivative at a given point
double derivativeAt(double x) {
    return -2.0 * x / std::pow(1.0 + x * x, 2);
}

// Define helper functions to calculate different f(x)
double f1(double x) {
    return 1.0 / (1.0 + x * x);
}

double f2(double x) {
    return sin(M_PI * x);
}

void plotNaturalCubicBSpline() {
    // Define the nodes and corresponding function values for f1
    std::vector<double> points = {-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<double> values(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        values[i] = f1(points[i]);
    }

    // Create a NaturalCubicBSpline object
    NaturalCubicBSpline spline(points, values);

    // Prepare the Gnuplot script
    std::ofstream gnuplot_script("NCBspline_plot.gp");
    gnuplot_script << "set terminal png size 800,600\n";
    gnuplot_script << "set output 'NCBspline_plot.png'\n";
    gnuplot_script << "set grid\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'f(x)'\n";
    gnuplot_script << "set title 'Natural Cubic BSpline'\n";
    gnuplot_script << "plot '-' with lines title 'BSpline'\n";

    // Calculate and write the B-spline interpolation results
    double x = -6.0;
    while (x <= 6.0) {
        gnuplot_script << x << " " << spline(x) << "\n";
        x += 0.1;
    }
    gnuplot_script.close();
    system("gnuplot NCBspline_plot.gp");
}

void plotClampedCubicBSpline() {
    // Define the nodes and corresponding function values for f1
    std::vector<double> points = {-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<double> values(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        values[i] = f1(points[i]);
    }

    // Define the derivative values at the start and end points
    double startDeriv = derivativeAt(points[2]); 
    double endDeriv = derivativeAt(points[points.size() - 3]); 

    // Create a ClampedCubicBSpline object
    ClampedCubicBSpline spline(points, values, startDeriv, endDeriv);

    // Prepare the Gnuplot script
    std::ofstream gnuplot_script("CCBspline_plot.gp");
    gnuplot_script << "set terminal png size 800,600\n";
    gnuplot_script << "set output 'CCBspline_plot.png'\n";
    gnuplot_script << "set grid\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'f(x)'\n";
    gnuplot_script << "set title 'Clamped Cubic BSpline'\n";
    gnuplot_script << "plot '-' with lines title 'BSpline'\n";

    // Calculate and write the B-spline interpolation results
    double x = -6.0;
    while (x <= 6.0) {
        gnuplot_script << x << " " << spline(x) << "\n";
        x += 0.1;
    }
    gnuplot_script.close();
    system("gnuplot CCBspline_plot.gp");
}

void plotPeriodicCubicBSpline() {
    // Define the nodes and corresponding function values for f2
    std::vector<double> points = {-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2};
    std::vector<double> values(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        values[i] = f2(points[i]);
    }

    // Create a PeriodicCubicBSpline object
    PeriodicCubicBSpline spline(points, values);

    // Prepare the Gnuplot script
    std::ofstream gnuplot_script("PCBspline_plot.gp");
    gnuplot_script << "set terminal png size 800,600\n";
    gnuplot_script << "set output 'PCBspline_plot.png'\n";
    gnuplot_script << "set grid\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'f(x)'\n";
    gnuplot_script << "set title 'Periodic Cubic BSpline'\n";
    gnuplot_script << "plot '-' with lines title 'BSpline'\n";

    // Calculate and write the B-spline interpolation results
    double x = -1.0;
    while (x <= 1.0) {
        gnuplot_script << x << " " << spline(x) << "\n";
        x += 0.01;
    }
    gnuplot_script.close();
    system("gnuplot PCBspline_plot.gp");
}

int main() {
    plotNaturalCubicBSpline();
    plotClampedCubicBSpline();
    plotPeriodicCubicBSpline();
    return 0;
}