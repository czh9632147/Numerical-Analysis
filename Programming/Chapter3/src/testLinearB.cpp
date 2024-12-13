#include <fstream>
#include "Spline.hpp"

int main() {
    // Define the function to interpolate
    std::vector<double> t = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
    std::vector<double> values;
    for (double x : t) {
        values.push_back(1.0 / (1.0 + x * x));
    }

    // CubicBSpline spline(t, values);
    LinearBSpline spline(t, values);

    // Create a Gnuplot script
    std::ofstream gnuplot_script("plot.gp");
    gnuplot_script << "set terminal png size 800,600\n";
    gnuplot_script << "set output 'LBspline_plot.png'\n";
    gnuplot_script << "set grid\n";
    gnuplot_script << "set xlabel 'x'\n";
    gnuplot_script << "set ylabel 'f(x)'\n";
    gnuplot_script << "set title 'Natural Cubic BSpline'\n";
    gnuplot_script << "plot '-' with lines title 'BSpline'\n";

    // Write the data to the Gnuplot script
    for (size_t i = 0; i < t.size(); ++i) {
        gnuplot_script << t[i] << " " << values[i] << "\n";
    }

    // Write the spline data to the Gnuplot script
    double x = -5.0;
    while (x <= 5.0) {
        gnuplot_script << x << " " << spline(x) << "\n";
        x += 0.1;
    }

    gnuplot_script.close();

    // Run Gnuplot
    system("gnuplot plot.gp");

    return 0;
}