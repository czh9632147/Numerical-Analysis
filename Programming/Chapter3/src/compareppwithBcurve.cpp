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

// Function to calculate error between the ppspline curve and the Bspline
double calculateError(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
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
                // Create Cubic ppSplines for x(t) and y(t)
                CubicSpline ppsplineX(tPoints, xPoints, i, xstartd, xendd);
                CubicSpline ppsplineY(tPoints, yPoints, i, ystartd, yendd);

                // Create Cubic B-Splines for x(t) and y(t)
                auto BsplineX = CubicBSplineFactory::createSpline(i, tPoints, xPoints, xstartd, xendd);
                auto BsplineY = CubicBSplineFactory::createSpline(i, tPoints, yPoints, ystartd, yendd);

                // Create filenames for data and scripts
                std::string dataFileName = "compare_spline_data_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".txt";
                std::string gnuplotScriptName = "plot_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".gp";
                std::string outputImageName = "compare_"+ std::to_string(i) + "heartspline_plot_" + std::to_string(numPoints) + ".png";
                std::string error3DImageName = "compare_"+ std::to_string(i) + "heartspline_error3d_" + std::to_string(numPoints) + ".png";

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
                    double pp_x = ppsplineX.evaluate(t);
                    double pp_y = ppsplineY.evaluate(t);
                    double B_x = BsplineX->eval(t);
                    double B_y = BsplineY->eval(t);
                    double error = calculateError(pp_x, pp_y, B_x, B_y);

                    dataFile << t << " " 
                            << pp_x << " " 
                            << pp_y << " "
                            << B_x << " "
                            << B_y << " "
                            << error << "\n";

                    maxError = std::max(maxError, error);
                    avgError += error;
                    pointCount++;
                }
                dataFile.close();

                avgError /= pointCount;

                // Write statistics to a separate file
                std::string statsFileName = "compare_heart_stats_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".txt";
                std::ofstream statsFile(statsFileName);
                statsFile << "Number of control points: " << numPoints << "\n";
                statsFile << "Spline type: B" << i << "\n";
                statsFile << "Maximum error: " << maxError << "\n";
                statsFile << "Average error: " << avgError << "\n";
                statsFile.close();

                // Generate 3D error visualization script
                std::string error3DScriptName = "compare_heart_plot_error3d_" + std::to_string(numPoints) + "_type" + std::to_string(i) + ".gp";
                std::ofstream error3DScript(error3DScriptName);
                error3DScript << "set terminal png size 1200,800\n";
                error3DScript << "set output '" << error3DImageName << "'\n";
                error3DScript << "set grid\n";
                error3DScript << "set title 'Compare Heart Curve Error Analysis (B" << i << ", " << numPoints << " points)'\n";
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