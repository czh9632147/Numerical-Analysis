#include <fstream>
#include "Spline.hpp"

// Function to generate t parameter points
void generateParameterPoints(std::vector<double>& tPoints, int numPoints) {
    double step = 2 * M_PI / numPoints;
    for (int i = 0; i < numPoints; ++i) {
        double t = i * step;
        tPoints.push_back(t);
    }
    tPoints.push_back(2 * M_PI);
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

            for (int i = 1; i <= 3; i++){
                // Create Cubic ppSplines for x(t) and y(t)
                CubicSpline splineX(tPoints, xPoints, i, xstartd, xendd);
                CubicSpline splineY(tPoints, yPoints, i, ystartd, yendd);

                // Create filenames for data and scripts
                std::string dataFileName = "spline_data_" + std::to_string(numPoints) + ".txt";
                std::string gnuplotScriptName = "plot_" + std::to_string(numPoints) + ".gp";
                std::string outputImageName = "pp"+ std::to_string(i) + "spline_plot_" + std::to_string(numPoints) + ".png";

                // Write spline points to a file for Gnuplot
                std::ofstream dataFile(dataFileName);
                if (!dataFile) {
                    throw std::runtime_error("Unable to open file for writing");
                }

                double step = 2 * M_PI / 1000; // Higher resolution for smoother plot
                for (double t = 0; t <= 2 * M_PI; t += step) {
                    double x = splineX.evaluate(t);
                    double y = splineY.evaluate(t);
                    dataFile << x << " " << y << "\n";
                }
                dataFile.close();

                // Generate Gnuplot script
                std::ofstream gnuplotScript(gnuplotScriptName);
                gnuplotScript << "set terminal png size 800,600\n";
                gnuplotScript << "set output '" << outputImageName << "'\n";
                gnuplotScript << "set grid\n";
                gnuplotScript << "plot '" << dataFileName << "' with lines linewidth 2 title 'Heart Curve (" 
                              << numPoints << " points), ppCubicType" << i << "'\n";
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
