#include <fstream>
#include "Spline.hpp"

std::pair<std::vector<double>, std::vector<double>> RiemannBall(const std::vector<double>& xBall, 
                                                                const std::vector<double>& yBall, 
                                                                const std::vector<double>& zBall) {
    int n = xBall.size();
    std::vector<double> xPlane(n), yPlane(n);
    for (int i = 0; i < n; i++) {
        xPlane[i] = 0.5 * xBall[i] / (1 - 0.5 * (zBall[i] + 1));
        yPlane[i] = 0.5 * yBall[i] / (1 - 0.5 * (zBall[i] + 1));
    }
    return {xPlane, yPlane};
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> 
invRiemannBall(const std::vector<double>& xPlane, const std::vector<double>& yPlane) {
    int n = xPlane.size();
    std::vector<double> xBall(n), yBall(n), zBall(n);
    for (int i = 0; i < n; i++) {
        xBall[i] = 2 * xPlane[i] / (xPlane[i] * xPlane[i] + yPlane[i] * yPlane[i] + 1);
        yBall[i] = 2 * yPlane[i] / (xPlane[i] * xPlane[i] + yPlane[i] * yPlane[i] + 1);
        zBall[i] = 2 * (xPlane[i] * xPlane[i] + yPlane[i] * yPlane[i]) / (xPlane[i] * xPlane[i] + yPlane[i] * yPlane[i] + 1) - 1;
    }
    return {xBall, yBall, zBall};
}

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
        double x = std::sin(std::cos(t) + 0.5) * std::cos(std::sin(2 * t));
        double y = std::sin(std::cos(t) + 0.5) * std::sin(std::sin(2 * t));
        double z = -std::cos(std::cos(t) + 0.5);
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
        std::vector<double> xBall;
        std::vector<double> yBall;
        std::vector<double> zBall;
        generateCurve(tPoints, xBall, yBall, zBall);
        std::pair<std::vector<double>, std::vector<double>> result = RiemannBall(xBall, yBall, zBall);
        std::vector<double> xPlane;
        std::vector<double> yPlane;
        xPlane = result.first;
        yPlane = result.second;

        try {
            // Create Periodic Cubic B-Splines for x(t) and y(t)
            auto splineX = CubicBSplineFactory::createSpline(1, tPoints, xPlane);
            auto splineY = CubicBSplineFactory::createSpline(1, tPoints, yPlane);

            std::string dataFileName = "spline_ball_data_2_" + std::to_string(numPoints) + ".txt";
            std::string gnuplotScriptName = "plot_ball_2_" + std::to_string(numPoints) + ".gp";
            std::string outputImageName = "B" + std::to_string(1) + "r5spline_ball_" + std::to_string(numPoints) + ".png";

            std::ofstream dataFile(dataFileName);
            if (!dataFile) {
                throw std::runtime_error("Unable to open file for writing");
            }

            std::vector<double> xx;
            std::vector<double> yy;
            double step = 2 * M_PI / 1000; // Higher resolution for smoother plot
            for (double t = 0; t <= 2 * M_PI; t += step) {
                double x = splineX->eval(t);
                double y = splineY->eval(t);
                xx.push_back(x);
                yy.push_back(y);
            }

            std::vector<double> X;
            std::vector<double> Y;
            std::vector<double> Z;
            std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> inv = invRiemannBall(xx, yy);
            X = std::get<0>(inv);
            Y = std::get<1>(inv);
            Z = std::get<2>(inv);
            for (int i = 0; i < X.size(); i++) {
                dataFile << X[i] << " " << Y[i] << " " << Z[i] << "\n";
            }
            dataFile.close();

            std::ofstream gnuplotScript(gnuplotScriptName);
            gnuplotScript << "set terminal png size 800,600\n";
            gnuplotScript << "set output '" << outputImageName << "'\n";
            gnuplotScript << "set view 30, 30, 1.5, 0.5\n"; // Set view angle for 3D plot
            gnuplotScript << "splot '" << dataFileName << "' with lines title 'r_5 (" << numPoints << " points), BType" << 1 << "'\n";
            gnuplotScript.close();

            std::string command = "gnuplot " + gnuplotScriptName;
            int result = system(command.c_str());
            if (result != 0) {
            std::cerr << "Failed to run Gnuplot script." << std::endl;
            } else {
                std::cout << "Plot generated: " << outputImageName << std::endl;
            }

        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    return 0;
}
