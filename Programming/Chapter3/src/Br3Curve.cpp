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

// Function to generate x and y points based on the curve
void generateCurve(const std::vector<double>& tPoints, std::vector<double>& xPoints, 
                        std::vector<double>& yPoints, std::vector<double>& zPoints) {
    for (double t : tPoints) {
        double x = std::sin(std::cos(t)) * std::cos(std::sin(t));
        double y = std::sin(std::cos(t)) * std::sin(std::sin(t));
        double z = std::cos(std::cos(t));
        xPoints.push_back(x);
        yPoints.push_back(y);
        zPoints.push_back(z);
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
        std::vector<double> zPoints;
        generateCurve(tPoints, xPoints, yPoints, zPoints);

        try {
            double xstartd = (-1) * std::cos(std::cos(tPoints[0])) * std::sin(tPoints[0]) * std::cos(std::sin(tPoints[0])) 
                           - std::sin(std::sin(tPoints[0])) * std::cos(tPoints[0]) * std::sin(std::cos(tPoints[0]));
            double xendd = (-1) * std::cos(std::cos(tPoints[tPoints.size()-1])) * std::sin(tPoints[tPoints.size()-1]) * std::cos(std::sin(tPoints[tPoints.size()-1])) 
                         - std::sin(std::sin(tPoints[tPoints.size()-1])) * std::cos(tPoints[tPoints.size()-1]) * std::sin(std::cos(tPoints[tPoints.size()-1]));
            double ystartd = (-1) * std::cos(std::cos(tPoints[0])) * std::sin(tPoints[0]) * std::sin(std::sin(tPoints[0])) 
                           + std::sin(std::cos(tPoints[0])) * std::cos(tPoints[0]) * std::cos(std::sin(tPoints[0]));
            double yendd = (-1) * std::cos(std::cos(tPoints[tPoints.size()-1])) * std::sin(tPoints[tPoints.size()-1]) * std::sin(std::sin(tPoints[tPoints.size()-1])) 
                         + std::sin(std::cos(tPoints[tPoints.size()-1])) * std::cos(tPoints[tPoints.size()-1]) * std::cos(std::sin(tPoints[tPoints.size()-1]));
            double zstartd = std::sin(std::cos(tPoints[0])) * std::sin(tPoints[0]);
            double zendd = std::sin(std::cos(tPoints[tPoints.size()-1])) * std::sin(tPoints[tPoints.size()-1]);

            for (int i = 1; i <= 3; i++){
                // Create Periodic Cubic B-Splines for x(t) and y(t)
                auto splineX = CubicBSplineFactory::createSpline(i, tPoints, xPoints, xstartd, xendd);
                auto splineY = CubicBSplineFactory::createSpline(i, tPoints, yPoints, ystartd, yendd);
                auto splineZ = CubicBSplineFactory::createSpline(i, tPoints, zPoints, zstartd, zendd);

                std::string dataFileName = "spline_data_" + std::to_string(numPoints) + ".txt";
                std::string gnuplotScriptName = "plot_" + std::to_string(numPoints) + ".gp";
                std::string outputImageName = "B" + std::to_string(i) + "r3spline_plot_" + std::to_string(numPoints) + ".png";

                std::ofstream dataFile(dataFileName);
                if (!dataFile) {
                    throw std::runtime_error("Unable to open file for writing");
                }

                double step = 2 * M_PI / 1000; // Higher resolution for smoother plot
                for (double t = 0; t <= 2 * M_PI; t += step) {
                    double x = splineX->eval(t);
                    double y = splineY->eval(t);
                    double z = splineZ->eval(t);
                    dataFile << x << " " << y << " " << z << "\n";
                }
                dataFile.close();

                std::ofstream gnuplotScript(gnuplotScriptName);
                gnuplotScript << "set terminal png size 800,600\n";
                gnuplotScript << "set output '" << outputImageName << "'\n";
                gnuplotScript << "set view 60, 30, 1.5, 0.5\n"; // Set view angle for 3D plot
                gnuplotScript << "splot '" << dataFileName << "' with lines title 'r_3 (" << numPoints << " points), BType" << i << "'\n";
                gnuplotScript.close();

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
