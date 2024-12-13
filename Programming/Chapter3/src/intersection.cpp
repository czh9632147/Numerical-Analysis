#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <tuple>

class CurveSelfIntersection {
private:
    struct Point3D {
        double x, y, z;
    };

    // Generate points on a curve in 3D space
    std::vector<Point3D> generateCurvePoints(int samples) {
        std::vector<Point3D> points;
        for (double t = 0; t <= 2 * M_PI; t += 2 * M_PI / samples) {
            double x = std::sin(std::cos(t)) * std::cos(std::sin(t));
            double y = std::sin(std::cos(t)) * std::sin(std::sin(t));
            double z = std::cos(std::cos(t));
            points.push_back({x, y, z});
        }
        return points;
    }

    // Compute the distance between two line segments
    double lineSegmentDistance(const Point3D& p1, const Point3D& p2, 
                                const Point3D& p3, const Point3D& p4) {
        const double eps = 1e-10;
        auto subtract = [](const Point3D& a, const Point3D& b) {
            return Point3D{a.x - b.x, a.y - b.y, a.z - b.z};
        };

        auto dot = [](const Point3D& a, const Point3D& b) {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        };

        auto cross = [](const Point3D& a, const Point3D& b) {
            return Point3D{
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x
            };
        };

        Point3D u = subtract(p2, p1);
        Point3D v = subtract(p4, p3);
        Point3D w = subtract(p3, p1);

        double a = dot(u, u);
        double b = dot(u, v);
        double c = dot(v, v);
        double d = dot(u, w);
        double e = dot(v, w);
        double D = a*c - b*b;

        double sc, tc;
        if (D < eps) {
            sc = 0.0;
            tc = (b > c ? d/b : e/c);
        } else {
            sc = (b*e - c*d) / D;
            tc = (a*e - b*d) / D;
        }

        Point3D dP = subtract(w, 
            Point3D{sc * u.x + tc * v.x, 
                    sc * u.y + tc * v.y, 
                    sc * u.z + tc * v.z});

        return std::sqrt(dot(dP, dP));
    }

public:
    bool hasSelfIntersection(int samples = 1000, double threshold = 1e-5, 
                              int minSegmentDistance = 10) {
        auto points = generateCurvePoints(samples);
        int n = points.size();

        // Check for self-intersection
        for (int i = 0; i < n; ++i) {
            for (int j = i + minSegmentDistance; j < n; ++j) {
                double dist = lineSegmentDistance(
                    points[i], points[(i+1)%n], 
                    points[j], points[(j+1)%n]
                );
                if (dist < threshold) {
                    return true;
                }
            }
        }
        return false;
    }
};

int main() {
    CurveSelfIntersection intersectionChecker;
    bool result = intersectionChecker.hasSelfIntersection();
    std::cout << "Curve self-intersection: " 
              << (result ? "Yes" : "No") << std::endl;
    return 0;
}