#include <fstream>
#include <cmath>
#include "Spline.hpp"

// Define a helper function to calculate the value of f(x)
double f(double x) {
    return 1.0 / (1.0 + x * x);
}

// Define the specific function f(x) = 1 / (1 + x^2)
class MyFunction : public Function {
public:
    double eval(const double& x) override {
        return 1.0 / (1.0 + pow(x, 2));
    }
};

int main() {
    // Define the points for the spline
    std::vector<double> prepoints = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    std::vector<double> qpoints = {-6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
    std::vector<double> qvalues(qpoints.size() + 1);
    for (size_t i = 0; i < qpoints.size() + 1; ++i) {
        qvalues[i] = f((prepoints[i] + prepoints[i + 1]) / 2);
    }

    std::vector<double> cpoints = {-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 
                                    0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7};
    std::vector<double> cvalues(cpoints.size());
    for (size_t i = 0; i < cpoints.size(); ++i) {
        cvalues[i] = f(cpoints[i]);
    }

    // Create instances of QuadraticBSpline and NaturalCubicBSpline
    MyFunction ff;
    QuadraticBSpline quadBSpline(qpoints, qvalues);
    NaturalCubicBSpline cubicBSpline(cpoints, cvalues);

    // Generate a gnuplot script to plot the splines and the original function
    std::ofstream gnuplot_script1("plot_splines.gp");
    gnuplot_script1 << "set terminal png size 800,600\n";
    gnuplot_script1 << "set output 'QCBspline_plot.png'\n";
    gnuplot_script1 << "set grid\n";
    gnuplot_script1 << "set xrange [-5:5]\n";
    gnuplot_script1 << "set yrange [0:1.1]\n";
    gnuplot_script1 << "plot '-' with lines title 'Original Function', '-' with lines title 'Quadratic B-Spline', '-' with lines title 'Natural Cubic B-Spline'\n";

    // Plot the original function
    for (double x = -5; x <= 5; x += 0.1) {
        gnuplot_script1 << x << " " << ff.eval(x) << "\n";
    }
    gnuplot_script1 << "e\n";

    // Plot the quadratic B-spline
    for (double x = -5; x <= 5; x += 0.1) {
        gnuplot_script1 << x << " " << quadBSpline.eval(x) << "\n";
    }
    gnuplot_script1 << "e\n";

    // Plot the natural cubic B-spline
    for (double x = -5; x <= 5; x += 0.1) {
        gnuplot_script1 << x << " " << cubicBSpline.eval(x) << "\n";
    }
    gnuplot_script1 << "e\n";

    gnuplot_script1.close();

    // Calculate and print the interpolation error at the specified sites
    std::vector<double> sites = {-3.5, -3, -0.5, 0, 0.5, 3, 3.5};
    std::cout << "Error at specified sites:\n";
    for (double site : sites) {
        double quadError = fabs(quadBSpline.eval(site) - f(site));
        double cubicError = fabs(cubicBSpline.eval(site) - f(site));
        std::cout << "Site: " << site << ", Quadratic Error: " << quadError << ", Cubic Error: " << cubicError << std::endl;
    }

    // Generate a gnuplot script to plot the interpolation error
    std::ofstream gnuplot_script2("plot_error.gp");
    gnuplot_script2 << "set terminal png size 800,600\n";
    gnuplot_script2 << "set output 'error_plot.png'\n";
    gnuplot_script2 << "set grid\n";
    gnuplot_script2 << "set xrange [-5:5]\n";
    gnuplot_script2 << "set yrange [0:0.025]\n";
    gnuplot_script2 << "plot '-' with lines title 'Quadratic Error', '-' with lines title 'Cubic Error'\n";

    // Plot the quadratic error
    for (double x = -5; x <= 5; x += 0.1) {
        double error_quad = fabs(quadBSpline.eval(x) - f(x));
        gnuplot_script2 << x << " " << error_quad << "\n";
    }
    gnuplot_script2 << "e\n";

    // Plot the cubic error
    for (double x = -5; x <= 5; x += 0.1) {
        double error_cubic = fabs(cubicBSpline.eval(x) - f(x));
        gnuplot_script2 << x << " " << error_cubic << "\n";
    }
    gnuplot_script2 << "e\n";

    gnuplot_script2.close();

    system("gnuplot plot_splines.gp");
    system("gnuplot plot_error.gp");
    return 0;
}