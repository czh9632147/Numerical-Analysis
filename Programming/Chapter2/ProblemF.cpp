#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Interpolation.hpp"

// Generate marker points on the heart curve
std::vector<BezierInterpolation::Point> generateMarkers(int m) {
    std::vector<BezierInterpolation::Point> markers;
    for (int i = 0; i <= m; ++i) {
        double t = i * (2 * M_PI / m); // Full range for heart shape
        double x = 16 * std::pow(std::sin(t), 3);
        double y = 13 * std::cos(t) - 5 * std::cos(2 * t) - 2 * std::cos(3 * t) - std::cos(4 * t);
        markers.emplace_back(x, y);
    }
    return markers;
}

// Calculate tangents at marker points
std::vector<BezierInterpolation::Point> calculateTangents(const std::vector<BezierInterpolation::Point>& markers) {
    std::vector<BezierInterpolation::Point> tangents;
    int n = markers.size();
    for (int i = 0; i < n; ++i) {
        BezierInterpolation::Point tangent;
        if (i > 0 && i < n - 1) {
            tangent = (markers[i + 1] - markers[i - 1]) / 2.0;
        } else if (i == 0) {
            tangent = markers[1] - markers[0];
        } else {
            tangent = markers[n - 1] - markers[n - 2];
        }
        tangents.push_back(tangent);
    }
    return tangents;
}

// Generate Gnuplot script for plotting
void generateGnuplotScript(const BezierInterpolation& bezierInterp,
                            const std::vector<BezierInterpolation::CubicBezier>& bezierCurves,
                            const std::string& filename) {
    std::ofstream file(filename);
    file << "set terminal pngcairo\n";
    file << "set output 'heart_approximation.png'\n";
    file << "set title 'Heart Approximation with Cubic Bezier Curves'\n";
    file << "plot '-' with lines title 'Bezier Curves'\n";

    for (const auto& curve : bezierCurves) {
        for (double t = 0; t <= 1.0; t += 0.01) {
            auto point = bezierInterp.evaluateBezier(curve, t);
            file << point.x << " " << point.y << "\n";
        }
        file << "\n";
    }

    file << "e\n";
    file.close();
}

int main() {
    std::vector<int> m_values = {10, 40, 160}; // Different m values for testing

    for (int m : m_values) {
        auto markers = generateMarkers(m);
        auto tangents = calculateTangents(markers);

        BezierInterpolation bezierInterp;
        auto bezierCurves = bezierInterp.fitCurve(markers, tangents);

        std::string scriptFilename = "heart_script_" + std::to_string(m) + ".gp";
        generateGnuplotScript(bezierInterp, bezierCurves, scriptFilename);

        std::cout << "Generated Gnuplot script: " << scriptFilename << " for m = " << m << "\n";
    }

    return 0;
}