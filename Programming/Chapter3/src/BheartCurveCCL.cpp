#include <fstream>
#include "Spline.hpp"

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

// Function to calculate cumulative chordal lengths
std::pair<double, std::vector<double>> calculateCumulativeChordalLengths(const std::vector<double>& tPoints, const std::vector<double>& xPoints, 
                                                                         const std::vector<double>& yPoints) {
    double limit;
    std::vector<double> newtPoints;
    newtPoints.push_back(0);
    for (int i = 1; i < tPoints.size(); i++) {
        double dx = xPoints[i] - xPoints[i-1];
        double dy = yPoints[i] - yPoints[i-1];
        double chordalLength = std::sqrt(dx * dx + dy * dy);
        newtPoints.push_back(newtPoints[i-1] + chordalLength);
        if (i == tPoints.size() / 2) {
            limit = newtPoints[i];
        }
    }
    return std::make_pair(limit, newtPoints);
}

int main() {
    std::vector<int> numPointsArray = {10, 40, 160};

    for (int numPoints : numPointsArray) {
        // Generate t parameter points
        std::vector<double> tPoints;
        generateParameterPoints(tPoints, numPoints);
        std::vector<double> xPoints;
        std::vector<double> yPoints;
        generateHeartCurve(tPoints, xPoints, yPoints);
        double limit;
        std::pair<double, std::vector<double>> result = calculateCumulativeChordalLengths(tPoints, xPoints, yPoints);
        limit = result.first;
        std::vector<double> newtPoints = result.second;

        try {
            double xstartd = 48 * std::pow(std::sin(tPoints[0]), 2) * std::cos(tPoints[0]);
            double xendd = 48 * std::pow(std::sin(tPoints[tPoints.size()-1]), 2) * std::cos(tPoints[tPoints.size()-1]);
            double ystartd = (-13) * std::sin(tPoints[0]) + 10 * std::sin(2 * tPoints[0]) 
                           + 6 * std::sin(3 * tPoints[0]) + 4 * std::sin(4 * tPoints[0]);
            double yendd = (-13) * std::sin(tPoints[tPoints.size()-1]) + 10 * std::sin(2 * tPoints[tPoints.size()-1]) 
                           + 6 * std::sin(3 * tPoints[tPoints.size()-1]) + 4 * std::sin(4 * tPoints[tPoints.size()-1]);

            for (int i = 1; i <= 3; i++) {
                // Create Periodic Cubic B-Splines for x(t) and y(t)
                auto splineX = CubicBSplineFactory::createSpline(i, newtPoints, xPoints, xstartd, xendd);
                auto splineY = CubicBSplineFactory::createSpline(i, newtPoints, yPoints, ystartd, yendd);

                // Create filenames for data and scripts
                std::string dataFileName = "spline_data_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".txt";
                std::string gnuplotScriptName = "plot_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".gp";
                std::string outputImageName = "B"+ std::to_string(i) + "CCLheartspline_plot_" + std::to_string(numPoints) + ".png";
                std::string error3DImageName = "B"+ std::to_string(i) + "CCLheartspline_error3d_" + std::to_string(numPoints) + ".png";

                // Write spline points and error data to a file for Gnuplot
                std::ofstream dataFile(dataFileName);
                if (!dataFile) {
                    throw std::runtime_error("Unable to open file for writing");
                }

                // Write header
                dataFile << "# t x_original y_original x_spline y_spline error\n";

                // Calculate and write points with higher resolution
                double step = limit / 1000;
                double maxError = 0;
                double avgError = 0;
                int pointCount = 0;

                for (double s = 0; s <= limit; s += step) {
                    double x_orig, y_orig;
                    getOriginalPoint(s, x_orig, y_orig);
                    double x_spline = splineX->eval(s);
                    double y_spline = splineY->eval(s);
                    double error = calculateError(x_spline, y_spline);

                    dataFile << s << " " 
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

                // Write statistics to a separate file
                std::string statsFileName = "heart_stats_CCL_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".txt";
                std::ofstream statsFile(statsFileName);
                statsFile << "Number of control points: " << numPoints << "\n";
                statsFile << "Spline type: B" << i << "\n";
                statsFile << "Maximum error: " << maxError << "\n";
                statsFile << "Average error: " << avgError << "\n";
                statsFile.close();

                // Generate original 2D plot script
                std::ofstream gnuplotScript(gnuplotScriptName);
                gnuplotScript << "set terminal png size 800,600\n";
                gnuplotScript << "set output '" << outputImageName << "'\n";
                gnuplotScript << "set grid\n";
                gnuplotScript << "set title 'Heart Curve (CCL B" << i << ", " << numPoints << " points)'\n";
                gnuplotScript << "plot '" << dataFileName << "' using 4:5 with lines title 'Fitted Curve', \\\n";
                gnuplotScript << "     '" << dataFileName << "' using 2:3 with lines title 'Original Curve'\n";
                gnuplotScript.close();

                // Generate 3D error visualization script
                std::string error3DScriptName = "heart_plot_error3d_CCL_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".gp";
                std::ofstream error3DScript(error3DScriptName);
                error3DScript << "set terminal png size 1200,800\n";
                error3DScript << "set output '" << error3DImageName << "'\n";
                error3DScript << "set grid\n";
                error3DScript << "set title 'Heart Curve Error Analysis (CCL B" << i << ", " << numPoints << " points)'\n";
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