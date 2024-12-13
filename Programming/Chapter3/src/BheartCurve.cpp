#include <fstream>
#include "Spline.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>

// Function to generate t parameter points
void generateParameterPoints(std::vector<double>& tPoints, int numPoints) {
    double step = 2 * M_PI / numPoints;
    for (int i = 0; i < numPoints * 2; ++i) {
        double t = i * step;
        tPoints.push_back(t);
    }
    tPoints.push_back(4 * M_PI);
}

// Function to generate x and y points based on the heart curve
void generateHeartCurve(const std::vector<double>& tPoints, std::vector<double>& xPoints, std::vector<double>& yPoints) {
    for (double t : tPoints) {
        double x = 16 * std::pow(std::sin(t), 3);
        double y = 13 * std::cos(t) - 5 * std::cos(2 * t) - 2 * std::cos(3 * t) - std::cos(4 * t);
        xPoints.push_back(x);
        yPoints.push_back(y);
    }
}

// Function to calculate the original curve point at any t
void getOriginalPoint(double t, double& x, double& y) {
    x = 16 * std::pow(std::sin(t), 3);
    y = 13 * std::cos(t) - 5 * std::cos(2 * t) - 2 * std::cos(3 * t) - std::cos(4 * t);
}

// Function to find the closest point on the curve to a given point (x0, y0) using the parametric equations
double findClosestT(double x0, double y0, double t_min = 0, double t_max = 2 * M_PI, double tol = 1e-6) {
    const double phi = (1 + std::sqrt(5)) / 2; // Golden ratio
    double t1 = t_max - (t_max - t_min) / phi;
    double t2 = t_min + (t_max - t_min) / phi;

    auto squaredDistance = [&](double t) {
        double x, y;
        getOriginalPoint(t, x, y);
        return std::pow(x - x0, 2) + std::pow(y - y0, 2);
    };

    while (t_max - t_min > tol) {
        if (squaredDistance(t1) > squaredDistance(t2)) {
            t_min = t1;
            t1 = t2;
            t2 = t_min + (t_max - t_min) / phi;
        } else {
            t_max = t2;
            t2 = t1;
            t1 = t_max - (t_max - t_min) / phi;
        }
    }

    return (t_min + t_max) / 2;
}

// Function to calculate error using the shortest distance to the original curve
double calculateError(double x_spline, double y_spline) {
    double t_closest = findClosestT(x_spline, y_spline);
    double x_orig, y_orig;
    getOriginalPoint(t_closest, x_orig, y_orig);
    return std::sqrt(std::pow(x_orig - x_spline, 2) + std::pow(y_orig - y_spline, 2));
}

int main() {
    // Array of numPoints values to iterate through
    std::vector<int> numPointsArray = {10, 40, 160};

    for (int numPoints : numPointsArray) {
        // Generate t parameter points
        std::vector<double> tPoints;
        generateParameterPoints(tPoints, numPoints);

        // Generate x and y points for the heart curve
        std::vector<double> xPoints;
        std::vector<double> yPoints;
        generateHeartCurve(tPoints, xPoints, yPoints);

        try {
            double xstartd = 48 * std::pow(std::sin(tPoints[0]), 2) * std::cos(tPoints[0]);
            double xendd = 48 * std::pow(std::sin(tPoints[tPoints.size()-1]), 2) * std::cos(tPoints[tPoints.size()-1]);
            double ystartd = (-13) * std::sin(tPoints[0]) + 10 * std::sin(2 * tPoints[0]) 
                           + 6 * std::sin(3 * tPoints[0]) + 4 * std::sin(4 * tPoints[0]);
            double yendd = (-13) * std::sin(tPoints[tPoints.size()-1]) + 10 * std::sin(2 * tPoints[tPoints.size()-1]) 
                           + 6 * std::sin(3 * tPoints[tPoints.size()-1]) + 4 * std::sin(4 * tPoints[tPoints.size()-1]);

            for (int i = 1; i <= 3; i++) {
                // Create Cubic B-Splines for x(t) and y(t)
                auto splineX = CubicBSplineFactory::createSpline(i, tPoints, xPoints, xstartd, xendd);
                auto splineY = CubicBSplineFactory::createSpline(i, tPoints, yPoints, ystartd, yendd);

                // Create filenames for data and scripts
                std::string dataFileName = "spline_data_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".txt";
                std::string gnuplotScriptName = "plot_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".gp";
                std::string outputImageName = "B"+ std::to_string(i) + "heartspline_plot_" + std::to_string(numPoints) + ".png";
                std::string error3DImageName = "B"+ std::to_string(i) + "heartspline_error3d_" + std::to_string(numPoints) + ".png";

                // Write spline points and error data to a file for Gnuplot
                std::ofstream dataFile(dataFileName);
                if (!dataFile) {
                    throw std::runtime_error("Unable to open file for writing");
                }

                // Write header
                dataFile << "# t x_original y_original x_spline y_spline error\n";

                // Calculate and write points with higher resolution
                double step = 2 * M_PI / 1000;
                double maxError = 0;
                double avgError = 0;
                int pointCount = 0;

                for (double t = 0; t <= 2 * M_PI; t += step) {
                    double x_orig, y_orig;
                    getOriginalPoint(t, x_orig, y_orig);
                    double x_spline = splineX->eval(t);
                    double y_spline = splineY->eval(t);
                    double error = calculateError(x_spline, y_spline);

                    dataFile << t << " " 
                            << x_orig << " " 
                            << y_orig << " "
                            << x_spline << " "
                            << y_spline << " "
                            << error << "\n";

                    maxError = std::max(maxError, error);
                    avgError += error;
                    pointCount++;
                }
                dataFile.close();

                avgError /= pointCount;

                // Calculate error at t = 0.5π and t = 1.5π
                double error_half_pi, error_three_halves_pi;
                {
                    // Calculate error at t = 0.5π
                    double x_orig_half_pi, y_orig_half_pi;
                    getOriginalPoint(0.5 * M_PI, x_orig_half_pi, y_orig_half_pi);
                    double x_spline_half_pi = splineX->eval(0.5 * M_PI);
                    double y_spline_half_pi = splineY->eval(0.5 * M_PI);
                    error_half_pi = calculateError(x_spline_half_pi, y_spline_half_pi);

                    // Calculate error at t = 1.5π
                    double x_orig_three_halves_pi, y_orig_three_halves_pi;
                    getOriginalPoint(1.5 * M_PI, x_orig_three_halves_pi, y_orig_three_halves_pi);
                    double x_spline_three_halves_pi = splineX->eval(1.5 * M_PI);
                    double y_spline_three_halves_pi = splineY->eval(1.5 * M_PI);
                    error_three_halves_pi = calculateError(x_spline_three_halves_pi, y_spline_three_halves_pi);
                }

                // Write statistics to a separate file
                std::string statsFileName = "heart_stats_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".txt";
                std::ofstream statsFile(statsFileName);
                statsFile << "Number of control points: " << numPoints << "\n";
                statsFile << "Spline type: B" << i << "\n";
                statsFile << "Maximum error: " << maxError << "\n";
                statsFile << "Average error: " << avgError << "\n";
                statsFile << "Error at t = 0.5π: " << error_half_pi << "\n";
                statsFile << "Error at t = 1.5π: " << error_three_halves_pi << "\n";
                statsFile.close();

                // Generate original 2D plot script
                std::ofstream gnuplotScript(gnuplotScriptName);
                gnuplotScript << "set terminal png size 800,600\n";
                gnuplotScript << "set output '" << outputImageName << "'\n";
                gnuplotScript << "set grid\n";
                gnuplotScript << "set title 'Heart Curve (B" << i << ", " << numPoints << " points)'\n";
                gnuplotScript << "plot '" << dataFileName << "' using 4:5 with lines linewidth 2 title 'Fitted Curve', \\\n";
                gnuplotScript << "     '" << dataFileName << "' using 2:3 with lines title 'Original Curve'\n";
                gnuplotScript.close();

                // Generate 3D error visualization script
                std::string error3DScriptName = "heart_plot_error3d_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".gp";
                std::ofstream error3DScript(error3DScriptName);
                error3DScript << "set terminal png size 1200,800\n";
                error3DScript << "set output '" << error3DImageName << "'\n";
                error3DScript << "set grid\n";
                error3DScript << "set title 'Heart Curve Error Analysis (B" << i << ", " << numPoints << " points)'\n";
                error3DScript << "set xlabel 'X'\n";
                error3DScript << "set ylabel 'Y'\n";
                error3DScript << "set zlabel 'Error'\n";
                error3DScript << "set view 60,30\n";
                error3DScript << "splot '" << dataFileName << "' using 4:5:6 with lines title 'Error Distribution'\n";
                error3DScript.close();

                // Run Gnuplot scripts
                std::string command = "gnuplot " + gnuplotScriptName;
                int result = system(command.c_str());
                if (result != 0) {
                    std::cerr << "Failed to run Gnuplot script for 2D plot." << std::endl;
                } else {
                    std::cout << "2D plot generated: " << outputImageName << std::endl;
                }

                command = "gnuplot " + error3DScriptName;
                result = system(command.c_str());
                if (result != 0) {
                    std::cerr << "Failed to run Gnuplot script for 3D error plot." << std::endl;
                } else {
                    std::cout << "3D error plot generated: " << error3DImageName << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    return 0;
}