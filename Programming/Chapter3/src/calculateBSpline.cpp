#include <iostream>
#include <fstream>
#include "Spline.hpp"

double calculateBSpline(const std::vector<double>& knots, 
                       int n, 
                       const std::vector<double>& coefs,
                       double t) {
    // Check input parameters
    int M = knots.size();
    int N = M - 2 * n;
    if (N < 2 || coefs.empty() || n < 1) {
        throw std::invalid_argument("Invalid input parameters");
    }
    
    // Check if t is within the valid range
    if (t < knots[n] || t > knots[n+N-1]) {
        throw std::out_of_range("Parameter t is out of valid range");
    }

    double result = 0.0;
    // Calculate the value of the B-spline at t
    for (int i = 1; i <= n+N-1; ++i) {
        class TempBSpline : public BSpline {
        public:
            using BSpline::basis;  // Allow access to the basis function
            int n;  // Degree of the B-spline

            TempBSpline(const std::vector<double>& k, const std::vector<double>& c, int s) : BSpline() {
                knots = k;
                coefs = c;
                size = s;
            }
                
            double getBasis(int i, double x) const override {
                return basis(i, n, x);
            }
        } temp(knots, coefs, N);
            
        temp.n = n;  // Set the degree of the B-spline
        result += coefs[i-1] * temp.basis(i, n, t);
    }
    
    return result;
}

int main() {
    try {
        // Get user input
        int n, M, N;
        std::cout << "Please enter the degree of the B-spline n: ";
        std::cin >> n;
        if (n < 1) throw std::invalid_argument("The degree must be greater than 0");

        std::cout << "Please enter the number of knots M: ";
        std::cin >> M;
        N = M - 2 * n;
        if (N < 2) throw std::invalid_argument("The number of knots is not sufficient for the given degree");

        // Input the knot sequence
        std::vector<double> knots(M);
        std::cout << "Please enter " << M << " knot values (in non-decreasing order):\n";
        for (int i = 0; i < M; ++i) {
            std::cin >> knots[i];
            if (i > 0 && knots[i] < knots[i-1]) {
                throw std::invalid_argument("The knot sequence must be non-decreasing");
            }
        }

        // Input control point coefficients
        int numCoefs = N + n - 1;  // Number of control points
        std::vector<double> coefs(numCoefs);
        std::cout << "Please enter " << numCoefs << " control point coefficients:\n";
        for (int i = 0; i < numCoefs; ++i) {
            std::cin >> coefs[i];
        }

        // Generate curve points
        std::ofstream dataFile("bspline_data.txt");
        double start = knots[n];      // Start of the valid interval
        double end = knots[n+N-1];    // End of the valid interval
        int numPoints = 2000;          // Number of sampling points
        double step = (end - start) / (numPoints - 1);

        for (double t = start; t <= end; t += step) {
            double y = calculateBSpline(knots, n, coefs, t);
            dataFile << t << " " << y << std::endl;
        }
        dataFile.close();

        // Generate gnuplot script
        std::ofstream gnuplot_script("plot.gp");
        gnuplot_script << "set terminal png size 800,600\n"
                      << "set output 'n_degree_BSpline_plot.png'\n"
                      << "set grid\n"
                      << "set title 'B-Spline Curve (n = " << n << ")'\n"
                      << "set xlabel 't'\n"
                      << "set ylabel 'S(t)'\n"
                      << "set style line 1 lc rgb '#0060ad' lt 1 lw 2\n"
                      << "set style line 2 lc rgb '#dd181f' lt 1 pt 7 ps 1.5\n"
                      << "plot 'bspline_data.txt' with lines ls 1 title 'B-Spline curve'\n";
        gnuplot_script.close();

        // Call gnuplot
        system("gnuplot plot.gp");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}