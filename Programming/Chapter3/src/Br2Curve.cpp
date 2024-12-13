#include <fstream>
#include "Spline.hpp"

// Function to generate t parameter points
void generateParameterPoints(std::vector<double>& tPoints, int numPoints) {
    double step = 6 * M_PI / numPoints;
    for (int i = 0; i < numPoints * 2; ++i) {
        double t = i * step;
        tPoints.push_back(t);
    }
}

// Function to generate x and y points based on the curve
void generateCurve(const std::vector<double>& tPoints, std::vector<double>& xPoints, std::vector<double>& yPoints) {
    for (double t : tPoints) {
        double x = std::sin(t) + t * std::cos(t);
        double y = std::cos(t) - t * std::sin(t);
        xPoints.push_back(x);
        yPoints.push_back(y);
    }
}

int main() {
    // Array of numPoints values to iterate through
    std::vector<int> numPointsArray = {10, 40, 160};

    for (int numPoints : numPointsArray) {
        // Generate t parameter points
        std::vector<double> tPoints;
        generateParameterPoints(tPoints, numPoints);

        // Generate x and y points for the curve
        std::vector<double> xPoints;
        std::vector<double> yPoints;
        generateCurve(tPoints, xPoints, yPoints);

        try {
            double xstartd = 2 * std::cos(tPoints[0]) - tPoints[0] * std::sin(tPoints[0]);
            double xendd = 2 * std::cos(tPoints[tPoints.size()-1]) - tPoints[tPoints.size()-1] * std::sin(tPoints[tPoints.size()-1]);
            double ystartd = (-2) * std::sin(tPoints[0]) - tPoints[0] * std::cos(tPoints[0]);
            double yendd = (-2) * std::sin(tPoints[tPoints.size()-1]) - tPoints[tPoints.size()-1] * std::cos(tPoints[tPoints.size()-1]);

            for (int i = 1; i <= 3; i++){
                // Create Periodic Cubic B-Splines for x(t) and y(t)
                auto splineX = CubicBSplineFactory::createSpline(i, tPoints, xPoints, xstartd, xendd);
                auto splineY = CubicBSplineFactory::createSpline(i, tPoints, yPoints, ystartd, yendd);

                // Create filenames for data and scripts
                std::string dataFileName = "spline_data_" + std::to_string(numPoints) + ".txt";
                std::string gnuplotScriptName = "plot_" + std::to_string(numPoints) + ".gp";
                std::string outputImageName = "B"+ std::to_string(i) + "r2spline_plot_" + std::to_string(numPoints) + ".png";

                // Write spline points to a file for Gnuplot
                std::ofstream dataFile(dataFileName);
                if (!dataFile) {
                    throw std::runtime_error("Unable to open file for writing");
                }

                double step = 6 * M_PI / 1000; // Higher resolution for smoother plot
                for (double t = 0; t <= 6 * M_PI; t += step) {
                    double x = splineX->eval(t);
                    double y = splineY->eval(t);
                    dataFile << x << " " << y << "\n";
                }
                dataFile.close();

                // Generate Gnuplot script
                std::ofstream gnuplotScript(gnuplotScriptName);
                gnuplotScript << "set terminal png size 800,600\n";
                gnuplotScript << "set output '" << outputImageName << "'\n";
                gnuplotScript << "set grid\n";
                gnuplotScript << "plot '" << dataFileName << "' with lines title 'r_2 (" << numPoints << " points), BType" << i << "'\n";
                gnuplotScript.close();

                // Run Gnuplot script to generate the plot
                std::string command = "gnuplot " + gnuplotScriptName;
                int result = system(command.c_str());
                if (result != 0) {
                    std::cerr << "Failed to run Gnuplot script." << std::endl;
                } else {
                    std::cout << "Plot generated: " << outputImageName << std::endl;
                }
            }

        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    return 0;
}