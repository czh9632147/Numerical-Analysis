#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>

// Function to calculate truncated power function (t - x)_+
double truncated_power(double t, double x) {
    return t > x ? (t - x) * (t - x) : 0.0;
}

// Function to calculate divided differences recursively
double divided_difference(const std::vector<double>& t, int i, int n, double x) {
    if (n == 0) {
        return truncated_power(t[i], x);
    } else {
        return (divided_difference(t, i + 1, n - 1, x) - divided_difference(t, i, n - 1, x));
    }
}

// Function to generate data for truncated power functions and divided differences
void generate_truncated_and_differences(const std::vector<double>& knots, const std::string& base_filename) {
    double x_min = knots.front() - 1.0;
    double x_max = knots.back() + 1.0;
    double step = 0.01;

    // Generate truncated power functions (first column)
    for (size_t i = 0; i < knots.size() - 1; ++i) {
        std::ofstream file(base_filename + "_truncated_" + std::to_string(i) + ".dat");
        for (double x = x_min; x <= x_max; x += step) {
            double value = truncated_power(knots[i], x);
            file << x << " " << value << "\n";
        }
        file.close();
    }

    // Generate first divided differences (second column)
    for (size_t i = 0; i < knots.size() - 2; ++i) {
        std::ofstream file(base_filename + "_diff1_" + std::to_string(i) + ".dat");
        for (double x = x_min; x <= x_max; x += step) {
            double value = divided_difference(knots, i, 1, x);
            file << x << " " << value << "\n";
        }
        file.close();
    }

    // Generate second divided differences (third column)
    for (size_t i = 0; i < knots.size() - 3; ++i) {
        std::ofstream file(base_filename + "_diff2_" + std::to_string(i) + ".dat");
        for (double x = x_min; x <= x_max; x += step) {
            double value = divided_difference(knots, i, 2, x);
            file << x << " " << value << "\n";
        }
        file.close();
    }

    // Generate third divided differences (fourth column)
    for (size_t i = 0; i < knots.size() - 4; ++i) {
        std::ofstream file(base_filename + "_diff3_" + std::to_string(i) + ".dat");
        for (double x = x_min; x <= x_max; x += step) {
            double value = divided_difference(knots, i, 3, x) / 2;
            file << x << " " << value << "\n";
        }
        file.close();
    }
}

// Function to write Gnuplot script for grid visualization with uniform axis ranges
void write_gnuplot_script(const std::string& base_filename, const std::string& output_file, size_t num_truncated, size_t num_diff1, size_t num_diff2) {
    std::ofstream script("n=2_plot_script.gp");
    script << "set terminal pngcairo size 1200,800\n";
    script << "set output '" << output_file << "'\n";
    script << "set multiplot layout " << num_truncated << ",3 title 'n=2: Truncated Power Functions and Divided Differences'\n";
    script << "set xlabel 'x'\n";
    script << "set ylabel 'Value'\n";
    script << "set xrange [-1:5]  # Set the x-axis range uniformly\n";
    script << "set yrange [-0.5:2.5]  # Set the y-axis range uniformly\n";

    // Plot truncated power functions (first column)
    for (size_t i = 0; i < num_truncated; ++i) {
        script << "plot '" << base_filename << "_truncated_" << i << ".dat' with lines title 'Truncated " << i << "'\n";
    }

    // Plot first divided differences (second column)
    for (size_t i = 0; i < num_diff1; ++i) {
        script << "plot '" << base_filename << "_diff1_" << i << ".dat' with lines title 'Diff1 " << i << "'\n";
    }

    // Plot second divided differences (third column)
    for (size_t i = 0; i < num_diff2; ++i) {
        script << "plot '" << base_filename << "_diff2_" << i << ".dat' with lines title 'Diff2 " << i << "'\n";
    }

    // Plot third divided differences (fourth column)
    for (size_t i = 0; i < num_diff2; ++i) {
        script << "plot '" << base_filename << "_diff3_" << i << ".dat' with lines title 'Diff3 " << i << "'\n";
    }

    script << "unset multiplot\n";
    script.close();
}

int main() {
    // Define knots (t_i)
    std::vector<double> knots = {0.0, 1.0, 2.0, 3.0, 4.0};

    // Generate data for n=1 truncated power functions and divided differences
    generate_truncated_and_differences(knots, "basis_n1");

    // Write Gnuplot script for visualization
    write_gnuplot_script("basis_n1", "basis_functions_n2.png", knots.size() - 1, knots.size() - 2, knots.size() - 3);

    // Call Gnuplot to generate plots
    int result = std::system("gnuplot n=2_plot_script.gp");
    if (result != 0) {
        std::cerr << "Error: Gnuplot execution failed.\n";
    }

    // Notify user
    std::cout << "Plot generated: basis_functions_n1.png\n";

    return 0;
}