#include <fstream>
#include "Spline.hpp"

void generateParameterPoints(std::vector<double>& tPoints, int numPoints) {
    double step = 2 * M_PI / numPoints;
    for (int i = 0; i < numPoints * 2; ++i) {
        double t = i * step;
        tPoints.push_back(t);
    }
    tPoints.push_back(4 * M_PI);
}

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

std::pair<double, std::vector<double>> calculateCumulativeChordalLengths(const std::vector<double>& tPoints, const std::vector<double>& xPoints, 
                                                                         const std::vector<double>& yPoints, const std::vector<double>& zPoints) {
    double limit;
    std::vector<double> newtPoints;
    newtPoints.push_back(0);
    for (int i = 1; i < tPoints.size(); i++) {
        double dx = xPoints[i] - xPoints[i-1];
        double dy = yPoints[i] - yPoints[i-1];
        double dz = zPoints[i] - zPoints[i-1];
        double chordalLength = std::sqrt(dx * dx + dy * dy + dz * dz);
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
        std::vector<double> tPoints;
        generateParameterPoints(tPoints, numPoints);
        std::vector<double> xPoints;
        std::vector<double> yPoints;
        std::vector<double> zPoints;
        generateCurve(tPoints, xPoints, yPoints, zPoints);
        double limit;
        std::pair<double, std::vector<double>> result = calculateCumulativeChordalLengths(tPoints, xPoints, yPoints, zPoints);
        limit = result.first;
        std::vector<double> newtPoints = result.second;

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

            for (int i = 1; i <= 3; i++) {
                // Create Periodic Cubic B-Splines for x(t) and y(t)
                auto splineX = CubicBSplineFactory::createSpline(i, newtPoints, xPoints, xstartd, xendd);
                auto splineY = CubicBSplineFactory::createSpline(i, newtPoints, yPoints, ystartd, yendd);
                auto splineZ = CubicBSplineFactory::createSpline(i, newtPoints, zPoints, zstartd, zendd);

                std::string dataFileName = "spline_data_" + std::to_string(numPoints) + ".txt";
                std::string gnuplotScriptName = "plot_" + std::to_string(numPoints) + ".gp";
                std::string outputImageName = "B" + std::to_string(i) + "CCLr3spline_plot_" + std::to_string(numPoints) + ".png";

                std::ofstream dataFile(dataFileName);
                if (!dataFile) {
                    throw std::runtime_error("Unable to open file for writing");
                }

                double step = limit / 1000; // Higher resolution for smoother plot
                for (double t = 0; t <= limit; t += step) {
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
                gnuplotScript << "splot '" << dataFileName << "' with lines title 'r_3 (cumulative chordal length" 
                              << numPoints << " points), BType" << i << "'\n";
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
