#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

// Abstract base class for interpolation methods
class Interpolation {
public:
    virtual ~Interpolation() {}

    // Pure virtual function to compute the interpolated value at x
    virtual double interpolate(double x) const = 0;

    // Pure virtual function to compute the interpolated value at x
    virtual std::string getExpression() const = 0;

protected:
    // Data members to store x and y values (f(x))
    std::vector<double> x_values;
    std::vector<double> y_values;

    Interpolation(const std::vector<double>& x, const std::vector<double>& y)
        : x_values(x), y_values(y) {}
};

// Newton Interpolation class
class NewtonInterpolation : public Interpolation {
public:
    NewtonInterpolation(const std::vector<double>& x, const std::vector<double>& y)
        : Interpolation(x, y) {}

    // Implement the interpolate function to compute the interpolated value at x.
    double interpolate(double x) const override {
        int n = x_values.size();
        std::vector<double> dd(n);

        // Calculate divided differences
        for (int i = 0; i < n; ++i) {
            dd[i] = y_values[i];
        }

        for (int j = 1; j < n; ++j) {
            for (int i = n - 1; i >= j; --i) {
                dd[i] = (dd[i] - dd[i - 1]) / (x_values[i] - x_values[i - j]);
            }
        }

        // Compute the Newton polynomial at x
        double result = dd[0];
        for (int i = 1; i < n; ++i) {
            double term = dd[i];
            for (int j = 0; j < i; ++j) {
                term *= (x - x_values[j]);
            }
            result += term;
        }

        return result;
    }

    std::string getExpression() const override {
        int n = x_values.size();
        std::vector<double> dd(n);

        // Calculate divided differences
        for (int i = 0; i < n; ++i) {
            dd[i] = y_values[i];
        }

        for (int j = 1; j < n; ++j) {
            for (int i = n - 1; i >= j; --i) {
                dd[i] = (dd[i] - dd[i - 1]) / (x_values[i] - x_values[i - j]);
            }
        }

        // Build expression string
        std::stringstream expr;
        expr << dd[0];
        for (int i = 1; i < n; ++i) {
            if (std::abs(dd[i]) < 1e-10) continue;

            // Add coefficient
            if (dd[i] >= 0) expr << " + ";
            else expr << " - ";
            expr << std::abs(dd[i]);

            // Add (x - x0)(x - x1)...(x - xi-1)
            for (int j = 0; j < i; ++j) {
                expr << "(x";
                if (x_values[j] >= 0) expr << " - ";
                else expr << " + ";
                expr << std::abs(x_values[j]) << ")";
            }
        }
        
        return expr.str();
    }
};

// Chebyshev Interpolation class
class ChebyshevInterpolation : public NewtonInterpolation {
public:
    ChebyshevInterpolation(const std::vector<double>& x, const std::vector<double>& y)
        : NewtonInterpolation(x, y) {}

    // Chebyshev nodes calculation (roots of Chebyshev polynomial of degree n)
    static std::vector<double> chebyshevNodes(int n, double a, double b) {
        std::vector<double> nodes(n);
        for (int i = 0; i < n; ++i) {
            // Modified the node calculation （froms smallest to largest）
            nodes[i] = 0.5 * (a + b) + (-0.5) * (b - a) * cos(M_PI * (2 * i + 1) / (2 * n));
        }
        return nodes;
    }
};

// Hermite Interpolation class
class HermiteInterpolation : public Interpolation {
public:
    // Constructor for Hermite Interpolation
    HermiteInterpolation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& y_prime)
        : Interpolation(x, y), derivatives(y_prime) {
        if (x.size() != y.size() || y.size() != y_prime.size()) {
            throw std::invalid_argument("Invalid input: x, y, and y_prime must have the same size.");
        }
    }

    // Implement the interpolate function to compute the interpolated value at x
    double interpolate(double x) const override {
        int n = x_values.size();
        std::vector<std::vector<double>> Q(2 * n, std::vector<double>(2 * n, 0.0));

        // Initialize Q table
        for (int i = 0; i < n; ++i) {
            Q[2 * i][0] = y_values[i];        // f(x_i)
            Q[2 * i + 1][0] = y_values[i];    // f(x_i) repeated
            Q[2 * i + 1][1] = derivatives[i]; // f'(x_i)
            if (i != 0) {
                Q[2 * i][1] = (Q[2 * i][0] - Q[2 * i - 1][0]) / (x_values[i] - x_values[i - 1]); // Divided difference
            }
        }

        // Fill the rest of the Q table
        for (int j = 2; j < 2 * n; ++j) {
            for (int i = j; i < 2 * n; ++i) {
                Q[i][j] = (Q[i][j - 1] - Q[i - 1][j - 1]) / (x_values[i / 2] - x_values[(i - j) / 2]);
            }
        }

        // Compute the Hermite polynomial at x
        double result = Q[0][0];
        for (int i = 1; i < 2 * n; ++i) {
            double term = Q[i][i];
            for (int j = 0; j < i; ++j) {
                term *= (x - x_values[j / 2]);
            }
            result += term;
        }

        return result;
    }

    // Compute the derivative of the Hermite interpolating polynomial at x
    double derivativeofx(double x) const {
        double delta = 1e-3;
        if (x == x_values.front()) {               // Forward difference
            return (interpolate(x + delta) - interpolate(x)) / delta;
        }
        else if (x == x_values.back()) {         // Backward difference
            return (interpolate(x) - interpolate(x - delta)) / delta;
        }
        else {                      // Central difference
            return (interpolate(x + delta) - interpolate(x - delta)) / (2 * delta);
        }
    }

    // Implement the getExpression function to return the Hermite polynomial expression
    std::string getExpression() const override {
        int n = x_values.size();
        std::vector<std::vector<double>> Q(2 * n, std::vector<double>(2 * n, 0.0));

        // Initialize Q table
        for (int i = 0; i < n; ++i) {
            Q[2 * i][0] = y_values[i];
            Q[2 * i + 1][0] = y_values[i];
            Q[2 * i + 1][1] = derivatives[i];
            if (i != 0) {
                Q[2 * i][1] = (Q[2 * i][0] - Q[2 * i - 1][0]) / (x_values[i] - x_values[i - 1]);
            }
        }

        // Fill the Q table
        for (int j = 2; j < 2 * n; ++j) {
            for (int i = j; i < 2 * n; ++i) {
                Q[i][j] = (Q[i][j - 1] - Q[i - 1][j - 1]) / (x_values[i / 2] - x_values[(i - j) / 2]);
            }
        }

        // Construct the polynomial expression as a string
        std::ostringstream expr;
        expr << std::fixed << std::setprecision(6);
        expr << Q[0][0];

        for (int i = 1; i < 2 * n; ++i) {
            if (Q[i][i] >= 0) {
                expr << " + " << Q[i][i];
            } else {
                expr << " - " << -Q[i][i];
            }

            for (int j = 0; j < i; ++j) {
                expr << "(x - " << x_values[j / 2] << ")";
            }
        }

        return expr.str();
    }

private:
    std::vector<double> derivatives;  // Derivatives (f'(x)) for Hermite interpolation
};

// Bezier Interpolation Class
class BezierInterpolation {
public:
    struct Point {
        double x, y;
        Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {}
        Point operator+(const Point& other) const {
            return Point(x + other.x, y + other.y);
        }
        Point operator-(const Point& other) const {
            return Point(x - other.x, y - other.y);
        }
        Point operator*(double scalar) const {
            return Point(x * scalar, y * scalar);
        }
        Point operator/(double scalar) const {
            return Point(x / scalar, y / scalar);
        }
    };

    struct CubicBezier {
        Point q0, q1, q2, q3; // Control points for the Bézier curve
    };

    std::vector<CubicBezier> fitCurve(const std::vector<Point>& markers, const std::vector<Point>& tangents) {
        int m = markers.size() - 1;
        std::vector<CubicBezier> bezierCurves;

        for (int j = 0; j < m; ++j) {
            Point q0 = markers[j];
            Point q1 = markers[j] + (tangents[j] / 3.0);
            Point q2 = markers[j + 1] - (tangents[j + 1] / 3.0);
            Point q3 = markers[j + 1];

            bezierCurves.push_back({q0, q1, q2, q3});
        }

        return bezierCurves;
    }

    Point evaluateBezier(const CubicBezier& bezier, double t) const {
        double u = 1 - t;
        return bezier.q0 * (u * u * u) +
               bezier.q1 * (3 * u * u * t) +
               bezier.q2 * (3 * u * t * t) +
               bezier.q3 * (t * t * t);
    }
};

#endif

