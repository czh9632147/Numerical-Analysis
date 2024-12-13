#ifndef SPLINE_HPP
#define SPLINE_HPP

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cassert>
#include "Function.hpp"
#include "Matrix.hpp"
#include "../libs/eigen/Eigen/Dense"


// Virtual base class representing a piecewise polynomial spline
class PiecewisePolynomialSpline {
protected:
    // Vector of knots, representing the segment points of the spline
    std::vector<double> knots;

public:
    // Constructor, accepts a vector of knots
    PiecewisePolynomialSpline(const std::vector<double>& knots) : knots(knots) {}

    // Virtual destructor
    virtual ~PiecewisePolynomialSpline() {}

    // Pure virtual function: compute the value of the spline at a given x
    virtual double evaluate(double x) const = 0;

    // Pure virtual function: print the expression of the spline
    virtual void printExpression() const = 0;

    // Return the knots of the spline
    const std::vector<double>& getKnots() const {
        return knots;
    }
};

// Linear Spline class that inherits from PiecewisePolynomialSpline
class LinearSpline : public PiecewisePolynomialSpline {
private:
    // Private member: stores the knots and the corresponding function values
    std::vector<double> values;

public:
    // Constructor that accepts knots and corresponding function values
    LinearSpline(const std::vector<double>& knots, const std::vector<double>& values)
        : PiecewisePolynomialSpline(knots), values(values) {
        if (knots.size() != values.size()) {
            throw std::invalid_argument("Knots and values vectors must have the same size");
        }
        if (knots.size() < 2) {
            throw std::invalid_argument("At least two knots are required for a linear spline");
        }
    }

    // Implementation of the evaluate method
    double evaluate(double x) const override {
        if (x < knots.front() || x > knots.back()) {
            throw std::out_of_range("x is out of the bounds of the spline");
        }

        // Find the interval where x lies
        for (size_t i = 0; i < knots.size() - 1; ++i) {
            if (x >= knots[i] && x <= knots[i + 1]) {
                double dx = knots[i + 1] - knots[i];
                double dy = values[i + 1] - values[i];
                double slope = dy / dx;
                double intercept = values[i] - slope * knots[i];
                return slope * x + intercept;
            }
        }
        throw std::runtime_error("Failed to evaluate spline at x");
    }

    double operator()(double x) const {
        return evaluate(x);
    }

    // Implementation of the printExpression method
    void printExpression() const override {
        for (size_t i = 0; i < knots.size() - 1; ++i) {
            double dx = knots[i + 1] - knots[i];
            double dy = values[i + 1] - values[i];
            double slope = dy / dx;
            double intercept = values[i] - slope * knots[i];

            std::cout << "Segment " << i + 1 << ": "
                      << slope << " * x + " << intercept
                      << " for x in [" << knots[i] << ", " << knots[i + 1] << "]\n";
        }
    }
};

// Cubic Spline class that inherits from PiecewisePolynomialSpline
class CubicSpline : public PiecewisePolynomialSpline {
private:
    int type;
    double leftDerivative, rightDerivative;
    std::vector<double> values;       // Function values at the knots
    std::vector<double> a, b, c, d;   // Coefficients a, b, c, d for each segment

public:
    // Constructor
    CubicSpline(const std::vector<double>& knots, const std::vector<double>& values, int type, 
                double leftDerivative = 0.0, double rightDerivative = 0.0)
        : PiecewisePolynomialSpline(knots), values(values) {
        if (knots.size() != values.size() || knots.size() < 2) {
            throw std::invalid_argument("Invalid number of knots or function values.");
        }
    
        if (type == 1) {
            computeNaturalSpline();
        } else if (type == 2) {
            computeClampedSpline(leftDerivative, rightDerivative);
        } else if (type == 3) {
            computePeriodicSpline();
        } else {
            throw std::invalid_argument("Invalid type.");
        }
    }

    // Implementation of the evaluate method
    double evaluate(double x) const override {
        if (x < knots.front() || x > knots.back()) {
            throw std::out_of_range("x is out of the bounds of the spline.");
        }
        
        // Find the interval that x is in and calculate the polynomial value
        for (size_t i = 0; i < knots.size() - 1; ++i) {
            if (x >= knots[i] && x <= knots[i + 1]) {
                double dx = x - knots[i];
                return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
            }
        }
        throw std::runtime_error("Failed to evaluate spline at x.");
    }

    double operator()(double x) const {
        return evaluate(x);
    }

    // Implementation of the printExpression method
    void printExpression() const override {
        for (size_t i = 0; i < knots.size() - 1; ++i) {
            std::cout << "Segment " << i + 1 << ": "
                      << std::setprecision(4) << a[i] << " + "
                      << b[i] << " * (x - " << knots[i] << ") + "
                      << c[i] << " * (x - " << knots[i] << ")^2 + "
                      << d[i] << " * (x - " << knots[i] << ")^3 "
                      << "for x in [" << knots[i] << ", " << knots[i + 1] << "]\n";
        }
    }

private:
    // Coefficient calculation for a natural spline
    void computeNaturalSpline() {
        size_t n = knots.size() - 1;
        a = values;
        b.resize(n);
        c.resize(n + 1);
        d.resize(n);

        // Create a tridiagonal matrix
        std::vector<double> h(n);
        for (size_t i = 0; i < n; ++i) {
            h[i] = knots[i + 1] - knots[i];
        }

        // Build the system of equations
        std::vector<double> alpha(n);
        for (size_t i = 1; i < n; ++i) {
            alpha[i] = (3.0 / h[i]) * (a[i + 1] - a[i]) - (3.0 / h[i - 1]) * (a[i] - a[i - 1]);
        }

        // Solve the tridiagonal matrix
        std::vector<double> l(n + 1, 1.0), mu(n), z(n + 1, 0.0);
        for (size_t i = 1; i < n; ++i) {
            l[i] = 2.0 * (knots[i + 1] - knots[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        // Set boundary conditions
        l[n] = 1.0;
        z[n] = 0.0;
        c[n] = 0.0;

        // Back substitution to calculate b, c, d
        for (int j = n - 1; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
    }

    // Coefficient calculation for a clamped cubic spline
    void computeClampedSpline(double leftDerivative, double rightDerivative) {
        size_t n = knots.size() - 1;
        a = values;
        b.resize(n);
        c.resize(n + 1);
        d.resize(n);

        // Create a tridiagonal matrix
        std::vector<double> h(n);
        for (size_t i = 0; i < n; ++i) {
            h[i] = knots[i + 1] - knots[i];
        }

        // Build the system of equations, set boundary conditions
        std::vector<double> alpha(n + 1);
        alpha[0] = 3.0 * ((a[1] - a[0]) / h[0] - leftDerivative);
        alpha[n] = 3.0 * (rightDerivative - (a[n] - a[n - 1]) / h[n - 1]);
        for (size_t i = 1; i < n; ++i) {
            alpha[i] = (3.0 / h[i]) * (a[i + 1] - a[i]) - (3.0 / h[i - 1]) * (a[i] - a[i - 1]);
        }

        // Solve the tridiagonal matrix
        std::vector<double> l(n + 1), mu(n), z(n + 1);
        l[0] = 2.0 * h[0];
        mu[0] = 0.5;
        z[0] = alpha[0] / l[0];
        for (size_t i = 1; i < n; ++i) {
            l[i] = 2.0 * (knots[i + 1] - knots[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[n] = h[n - 1] * (2.0 - mu[n - 1]);
        z[n] = (alpha[n] - h[n - 1] * z[n - 1]) / l[n];
        c[n] = z[n];

        // Back substitution to calculate b, c, d
        for (int j = n - 1; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
    }

    // Coefficient calculation for a periodic spline
    void computePeriodicSpline() {
        size_t n = knots.size();
        a = values;
        b.resize(n - 1);
        c.resize(n - 1);
        d.resize(n - 1);
        
        std::vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = knots[i + 1] - knots[i];
        }

        // Build the right-hand side of the system of equations
        std::vector<double> alpha(n - 1);
        for (size_t i = 0; i < n - 2; ++i) {
            alpha[i] = 6.0 * ((a[i + 2] - a[i + 1]) / h[i + 1] - (a[i + 1] - a[i]) / h[i]) / (h[i + 1] + h[i]);
        }
        alpha[n - 2] = 6.0 * ((a[1] - a[0]) / h[0] - (a[n - 1] - a[n - 2]) / h[n - 2]) / (h[n - 2] + h[0]);

        std::vector<double> l(n - 1), mu(n - 1);
        for (size_t i = 0; i < n - 2; ++i) {
            l[i] = h[i + 1] / (h[i] + h[i + 1]);
            mu[i] = h[i] / (h[1] + h[i + 1]);
        }
        l[n - 2] = h[0] / (h[n - 2] + h[0]);
        mu[n - 2] = h[n - 2] / (h[n - 2] + h[0]);

        // Create a tridiagonal matrix
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n - 1, n - 1);
        A(0, 0) = 2.0;
        A(0, 1) = l[0];
        A(0, n - 2) = mu[0];
        for (size_t i = 1; i <= n - 3; ++i) {
            for (size_t j = i - 1; j <= i + 1; ++j) {
                if (j == i - 1) {
                    A(i, j) = mu[i];
                } else if (j == i) {
                    A(i, j) = 2.0;
                } else if (j == i + 1) {
                    A(i, j) = l[i];
                }
            }
        }
        A(n - 2, 0) = l[n - 2];
        A(n - 2, n - 3) = mu[n - 2];
        A(n - 2, n - 2) = 2.0;

        // Solve the system of equations
        Eigen::VectorXd x = A.fullPivLu().solve(Eigen::Map<Eigen::VectorXd>(alpha.data(), alpha.size()));
        std::vector<double> M(n);
        M[0] = x[n - 2];
        for (size_t i = 1; i <= n - 1; ++i) {
            M[i] = x[i - 1];
        }

        // Calculate b, c, d
        for (size_t i = 0; i < n - 1; ++i) {
            b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (M[i + 1] + 2.0 * M[i]) / 6.0;
            c[i] = M[i] / 2.0;
            d[i] = (M[i + 1] - M[i]) / (6.0 * h[i]);
        } 
    }
};

// Degree-4 piecewise polynomial spline class (leftFirstDerivative, leftSecondDerivative, rightFirstDerivative)
class Degree_4_pp : public PiecewisePolynomialSpline {
private:
    std::vector<double> values;
    std::vector<double> a, b, c, d, e;

public:
    Degree_4_pp(const std::vector<double>& knots, const std::vector<double>& values)
        : PiecewisePolynomialSpline(knots), values(values) {
        if (knots.size() != values.size() || knots.size() < 2) {
            throw std::invalid_argument("Invalid number of knots or function values.");
        }
        computeCoefficients();
    }

    double evaluate(double x) const override {
        if (x < knots.front() || x > knots.back()) {
            throw std::out_of_range("x is out of the bounds of the spline.");
        }

        for (size_t i = 0; i < knots.size() - 1; ++i) {
            if (x >= knots[i] && x <= knots[i + 1]) {
                double dx = x - knots[i];
                return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx + e[i] * dx * dx * dx * dx;
            }
        }
        throw std::runtime_error("Failed to evaluate spline at x.");
    }

    double operator()(double x) const {
        return evaluate(x);
    }

    void printExpression() const override {
        for (size_t i = 0; i < knots.size() - 1; ++i) {
            std::cout << "Segment " << i + 1 << ": "
                      << std::setprecision(4) << a[i] << " + "
                      << b[i] << " * (x - " << knots[i] << ") + "
                      << c[i] << " * (x - " << knots[i] << ")^2 + "
                      << d[i] << " * (x - " << knots[i] << ")^3 + "
                      << e[i] << " * (x - " << knots[i] << ")^4 "
                      << "for x in [" << knots[i] << ", " << knots[i + 1] << "]\n";
        }
    }

private:
    void computeCoefficients() {
        size_t n = knots.size();
        a = values;
        b.resize(n - 1);
        c.resize(n - 1);
        d.resize(n);
        e.resize(n - 1);

        std::vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = knots[i + 1] - knots[i];
        }

        // Build the right-hand side vector
        std::vector<double> alpha(2 * n - 1);
        alpha[0] = 0.0;
        for (size_t i = 1; i <= n - 2; ++i) {
            alpha[i] = (a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1];
        }
        for (size_t i = n - 1; i <= 2 * n - 2; ++i) {
            alpha[i] = 0.0;
        }

        Eigen::MatrixXd G1 = Eigen::MatrixXd::Zero(n - 1, n - 1);
        G1(0, 0) = 1.0;
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i; ++j) {
                G1(i, j) = h[j];
            }
        }

        Eigen::MatrixXd G2 = Eigen::MatrixXd::Zero(n - 1, n);
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i + 1; ++j) {
                if (j == i - 1) {
                    G2(i, j) = 1.25 * std::pow(h[i - 1], 2);
                } else if (j == i) {
                    G2(i, j) = 0.75 * std::pow(h[i - 1], 2) + 0.75 * std::pow(h[i], 2);
                } else if (j == i + 1) {
                    G2(i, j) = 0.25 * std::pow(h[i], 2);
                }
            }
        }

        Eigen::MatrixXd G3 = Eigen::MatrixXd::Zero(n, n - 1);
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i; ++j) {
                if (j == i - 1) {
                    G3(i, j) = 2.0;
                } else if (j == i) {
                    G3(i, j) = -2.0;
                }
            }
        }

        Eigen::MatrixXd G4 = Eigen::MatrixXd::Zero(n, n);
        G4(0, 0) = 1.0;
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i; ++j) {
                G4(i, j) = h[i - 1];
            }
        }
        G4(n - 1, n - 1) = 1.0;

        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(2 * n - 1, 2 * n - 1);
        G << G1, G2,
             G3, G4;

        // Solve the system of equations
        Eigen::VectorXd x = G.fullPivLu().solve(Eigen::Map<Eigen::VectorXd>(alpha.data(), alpha.size()));

        // Extract the coefficients
        for (size_t i = 0; i < n - 1; ++i) {
            c[i] = x[i];
            d[i] = x[i + n - 1];
        }
        d[n - 1] = x[2 * n - 2];

        for (size_t i = 0; i < n - 2; ++i) {
            b[i] = (a[i + 1] - a[i]) / h[i] - c[i] * h[i] - 0.25 * d[i + 1] * std::pow(h[i], 2) 
                 - 0.75 * d[i] * std::pow(h[i], 2);
            e[i] = 0.25 * (d[i + 1] - d[i]) / h[i];
        }
        b[n - 2] = b[n - 3] + 2 * c[n - 3] * h[n - 3] + 3 * d[n - 3] * std::pow(h[n - 3], 2) 
                 + 4 * e[n - 3] * std::pow(h[n - 3], 3);
        e[n - 2] = (a[n - 1] - a[n - 2] - b[n - 2] * h[n - 2] - c[n - 2] * std::pow(h[n - 2], 2) 
                 - d[n - 2] * std::pow(h[n - 2], 3)) / std::pow(h[n - 2], 4);
    }
};


// Degree-5 piecewise polynomial spline class
class Degree_5_pp : public PiecewisePolynomialSpline {
private:
    std::vector<double> values;
    std::vector<double> a, b, c, d, e, f;

public:
    Degree_5_pp(const std::vector<double>& knots, const std::vector<double>& values)
        : PiecewisePolynomialSpline(knots), values(values) {
        if (knots.size() != values.size() || knots.size() < 2) {
            throw std::invalid_argument("Invalid number of knots or function values.");
        }
        computeCoefficients();
    }

    double evaluate(double x) const override {
        if (x < knots.front() || x > knots.back()) {
            throw std::out_of_range("x is out of the bounds of the spline.");
        }

        for (size_t i = 0; i < knots.size() - 1; ++i) {
            if (x >= knots[i] && x <= knots[i + 1]) {
                double dx = x - knots[i];
                return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx 
                     + e[i] * dx * dx * dx * dx + f[i] * dx * dx * dx * dx * dx;
            }
        }
        throw std::runtime_error("Failed to evaluate spline at x.");
    }

    double operator()(double x) const {
        return evaluate(x);
    }

    void printExpression() const override {
        for (size_t i = 0; i < knots.size() - 1; ++i) {
            std::cout << "Segment " << i + 1 << ": "
                      << std::setprecision(4) 
                      << a[i] << " + "
                      << b[i] << " * (x - " << knots[i] << ") + "
                      << c[i] << " * (x - " << knots[i] << ")^2 + "
                      << d[i] << " * (x - " << knots[i] << ")^3 + "
                      << e[i] << " * (x - " << knots[i] << ")^4 + "
                      << f[i] << " * (x - " << knots[i] << ")^5 "
                      << "for x in [" << knots[i] << ", " << knots[i + 1] << "]\n";
        }
    }

private:
    void computeCoefficients() {
        size_t n = knots.size();
        a = values;
        b.resize(n - 1);
        c.resize(n - 1);
        d.resize(n);
        e.resize(n);
        f.resize(n - 1);

        std::vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = knots[i + 1] - knots[i];
        }

        // Build the right-hand side vector
        std::vector<double> alpha(3 * n - 1);
        alpha[0] = 0.0;
        for (size_t i = 1; i <= n - 2; ++i) {
            alpha[i] = (a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1];
        }
        for (size_t i = n - 1; i <= 3 * n - 2; ++i) {
            alpha[i] = 0.0;
        }

        Eigen::MatrixXd G1 = Eigen::MatrixXd::Zero(n - 1, n - 1);
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i; ++j) {
                G1(i, j) = h[j];
            }
        }

        Eigen::MatrixXd G2 = Eigen::MatrixXd::Zero(n - 1, n);
        G2(0, 0) = 1.0;
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i; ++j) {
                if (j == i - 1) {
                    G2(i, j) = 2 * std::pow(h[i - 1], 2);
                } else if (j == i) {
                    G2(i, j) = std::pow(h[i], 2);
                }
            }
        }

        Eigen::MatrixXd G3 = Eigen::MatrixXd::Zero(n - 1, n);
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i + 1; ++j) {
                if (j == i - 1) {
                    G3(i, j) = 2.2 * std::pow(h[i - 1], 3);
                } else if (j == i) {
                    G3(i, j) = 0.8 * std::pow(h[i - 1], 3) + 0.8 * std::pow(h[i], 3);
                } else if (j == i + 1) {
                    G3(i, j) = 0.2 * std::pow(h[i], 3);
                }
            }
        }

        Eigen::MatrixXd G4 = Eigen::MatrixXd::Zero(n + 1, n - 1);
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i; ++j) {
                if (j == i - 1) {
                    G4(i, j) = 2.0;
                } else if (j == i) {
                    G4(i, j) = -2.0;
                }
            }
        }

        Eigen::MatrixXd G5 = Eigen::MatrixXd::Zero(n + 1, n);
        for (int i = 1; i <= n - 2; ++i) {
            G5(i, i - 1) = 6 * h[i - 1];
        }
        G5(n - 1, n - 1) = 1.0;
        
        Eigen::MatrixXd G6 = Eigen::MatrixXd::Zero(n + 1, n);
        G6(0, 0) = 1.0;
        for (int i = 1; i <= n - 2; ++i) {
            for (int j = i - 1; j <= i; ++j) {
                if (j == i - 1) {
                    G6(i, j) = 8 * std::pow(h[i - 1], 2);
                } else if (j == i) {
                    G6(i, j) = 4 * std::pow(h[i - 1], 2);
                }
            }
        }
        G6(n, n - 1) = 1.0;

        Eigen::MatrixXd G7 = Eigen::MatrixXd::Zero(n - 1, n - 1);

        Eigen::MatrixXd G8 = Eigen::MatrixXd::Zero(n - 1, n);
        for (int i = 0; i <= n - 2; ++i) {
            for (int j = i; j <= i + 1; ++j) {
                if (j == i) {
                    G8(i, j) = 6.0;
                } else if (j == i + 1) {
                    G8(i, j) = -6.0;
                }
            }
        }

        Eigen::MatrixXd G9 = Eigen::MatrixXd::Zero(n - 1, n);
        for (int i = 0; i <= n - 2; ++i) {
            for (int j = i; j <= i + 1; ++j) {
                G9(i, j) = 12 * h[i];
            }
        }

        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(3 * n - 1, 3 * n - 1);
        G << G1, G2, G3,
             G4, G5, G6,
             G7, G8, G9;

        // Solve the system of equations
        Eigen::VectorXd x = G.fullPivLu().solve(Eigen::Map<Eigen::VectorXd>(alpha.data(), alpha.size()));

        // Extract the coefficients
        for (size_t i = 0; i < n - 1; ++i) {
            c[i] = x[i];
            d[i] = x[i + n - 1];
            e[i] = x[i + 2 * n - 1];
        }
        d[n - 1] = x[2 * n - 2];
        e[n - 1] = x[3 * n - 2];

        for (size_t i = 0; i < n - 2; ++i) {
            b[i] = (a[i + 1] - a[i]) / h[i] - c[i] * h[i] - d[i] * std::pow(h[i], 2) 
                 - 0.8 * e[i] * std::pow(h[i], 3) - 0.2 * e[i + 1] * std::pow(h[i], 3);
            f[i] = 0.2 * (e[i + 1] - e[i]) / h[i];
        }
        b[n - 2] = b[n - 3] + 2 * c[n - 3] * h[n - 3] + 3 * d[n - 3] * std::pow(h[n - 3], 2) 
                 + 4 * e[n - 3] * std::pow(h[n - 3], 3) + 5 * f[n - 3] * std::pow(h[n - 3], 4);
        f[n - 2] = (a[n - 1] - a[n - 2] - b[n - 2] * h[n - 2] - c[n - 2] * std::pow(h[n - 2], 2)
                 - d[n - 2] * std::pow(h[n - 2], 3) - e[n - 2] * std::pow(h[n - 2], 4)) / std::pow(h[n - 2], 5);

    }
};


// BSpline class
class BSpline {
protected:
    int size;
    std::vector<double> knots;
    std::vector<double> coefs;

    // Basis function
    double basis(int i, int k, double x) const {
        if (k == 0) {
            return (knots[i-1] < x && x <= knots[i]) ? 1.0 : 0.0;
        }
        
        double left = 0.0, right = 0.0;
        
        if (knots[i+k-1] - knots[i-1] != 0) {
            left = (x - knots[i-1]) / (knots[i+k-1] - knots[i-1]) 
                  * basis(i, k-1, x);
        }
        
        if (knots[i+k] - knots[i] != 0) {
            right = (knots[i+k] - x) / (knots[i+k] - knots[i]) 
                   * basis(i+1, k-1, x);
        }
        
        return left + right;
    }

    // The first derivative of the basis function
    double basisDeriv(int i, int k, double x) const {
        if (k == 0) return 0.0;
        
        return k * (basis(i, k-1, x) / (knots[i+k-1] - knots[i-1]) 
               - basis(i+1, k-1, x) / (knots[i+k] - knots[i]));
    }

    // The second derivative of the basis function
    double basisDeriv2(int i, int k, double x) const {
        if (k <= 1) return 0.0;
        
        return k * (basisDeriv(i, k-1, x) / (knots[i+k-1] - knots[i-1])
               - basisDeriv(i+1, k-1, x) / (knots[i+k] - knots[i]));
    }

public:
    BSpline() : size(0) {}
    virtual ~BSpline() = default;

    // Get the basis function value
    virtual double getBasis(int i, double x) const = 0;
    
    // Compute the value of the spline at x
    double eval(double x) const {
        double result = 0.0;
        for(int i = 0; i < size; i++) {
            result += coefs[i] * getBasis(i+1, x);
        }
        return result;
    }
    
    double operator()(double x) const {
        return eval(x);
    }

    // Get the knots vector
    const std::vector<double>& getKnots() const { return knots; }
    
    // Get the coefficients vector
    const std::vector<double>& getCoefs() const { return coefs; }
};

// Linear B-Spline
class LinearBSpline : public BSpline {
public:
    LinearBSpline(const std::vector<double>& points, 
                  const std::vector<double>& values) {
        if (points.size() != values.size()) {
            throw std::invalid_argument("Points and values must have same size");
        }
        
        size = points.size();
        knots.clear();
        knots.push_back(points.front() - 1);
        knots.insert(knots.end(), points.begin(), points.end());
        knots.push_back(points.back() + 1);
        coefs = values;
    }

    LinearBSpline(const std::vector<double>& points, const Function& f) {
        std::vector<double> values;
        for (int i = 0; i < points.size(); i++) {
            values.push_back(eval(points[i]));
        }
        *this = LinearBSpline(points, values);
    }

    double getBasis(int i, double x) const override {
        return basis(i, 1, x);
    }
};

// Quadratic B-Spline
class QuadraticBSpline : public BSpline {
public:
    QuadraticBSpline(const std::vector<double>& points, 
                     const std::vector<double>& values) {
        if (points.size() < 2) {
            throw std::invalid_argument("Need at least 2 points");
        }
        if (values.size() != points.size() + 1) {
            throw std::invalid_argument("Invalid number of values, need one more than points");
        }
        
        size = points.size();
        setupKnots(points);
        computeCoefficients(points, values);
    }

    QuadraticBSpline(const std::vector<double>& points, const Function& f) {
        std::vector<double> values;
        values.push_back(eval(points.front()));
        for (int i = 0; i < points.size() - 1; i++) {
            values.push_back(eval((points[i] + points[i+1]) / 2));
        }
        values.push_back(eval(points.back()));
        *this = QuadraticBSpline(points, values);
    }

    double getBasis(int i, double x) const override {
        return basis(i, 2, x);
    }

private:
    void setupKnots(const std::vector<double>& points) {
        knots.clear();
        knots.push_back(points.front() - 2);
        knots.push_back(points.front() - 1);
        knots.insert(knots.end(), points.begin(), points.end());
        knots.push_back(points.back() + 1);
        knots.push_back(points.back() + 2);
    }

    void computeCoefficients(const std::vector<double>& points, 
                           const std::vector<double>& values) {
        Matrix A(size + 1);
        ColVector b(size + 1);
        
        // Set the first point condition
        A[0][0] = getBasis(1, points.front());
        A[0][1] = getBasis(2, points.front());
        b[0] = values.front();
        
        // Set the middle points conditions
        for (int i = 0; i < size - 1; i++) {
            double mid = (points[i] + points[i+1]) / 2;
            A[i+1][i] = getBasis(i+1, mid);
            A[i+1][i+1] = getBasis(i+2, mid);
            A[i+1][i+2] = getBasis(i+3, mid);
            b[i+1] = values[i+1];
        }
        
        A[size][size-1] = getBasis(size, points.back());
        A[size][size] = getBasis(size+1, points.back());
        b[size] = values.back();
        
        Matrix result = Matrix::solve(A, b);
        coefs.clear();
        for (int i = 0; i <= size; i++) {
            coefs.push_back(result[i][0]);
        }
    }
};

// Cubic B-Spline
class CubicBSpline : public BSpline {
protected:
    void setupKnots(const std::vector<double>& points) {
        knots.clear();
        knots.push_back(points.front() - 3);
        knots.push_back(points.front() - 2);
        knots.push_back(points.front() - 1);
        knots.insert(knots.end(), points.begin(), points.end());
        knots.push_back(points.back() + 1);
        knots.push_back(points.back() + 2);
        knots.push_back(points.back() + 3);
    }

public:
    double getBasis(int i, double x) const override {
        return basis(i, 3, x);
    }
    
    virtual double getDerivative(int i, double x) const {
        return basisDeriv(i, 3, x);
    }
    
    virtual double getSecondDerivative(int i, double x) const {
        return basisDeriv2(i, 3, x);
    }
};

// Natural Cubic B-Spline
class NaturalCubicBSpline : public CubicBSpline {
public:
    NaturalCubicBSpline(const std::vector<double>& points, 
                        const std::vector<double>& values) {
        if (points.size() < 2) {
            throw std::invalid_argument("Need at least 2 points");
        }
        if (values.size() != points.size()) {
            throw std::invalid_argument("Points and values must have same size");
        }
        
        size = points.size();
        setupKnots(points);
        computeCoefficients(points, values);
    }

private:
    void computeCoefficients(const std::vector<double>& points, 
                           const std::vector<double>& values) {
        Matrix A(size + 2);
        ColVector b(size + 2);
        
        // Set the interpolation conditions
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 3; j++) {
                A[i][i+j] = getBasis(i+1+j, points[i]);
            }
            b[i] = values[i];
        }
        
        // Set the second derivative boundary conditions
        for (int j = 0; j < 3; j++) {
            A[size][j] = getSecondDerivative(1+j, points.front());
            A[size+1][size-1+j] = getSecondDerivative(size+j, points.back());
        }
        b[size] = b[size+1] = 0;
        
        Matrix result = Matrix::solve(A, b);
        coefs.clear();
        for (int i = 0; i <= size + 1; i++) {
            coefs.push_back(result[i][0]);
        }
    }
};

// Clamped Cubic B-Spline
class ClampedCubicBSpline : public CubicBSpline {
public:
    ClampedCubicBSpline(const std::vector<double>& points, 
                        const std::vector<double>& values,
                        double startDeriv, double endDeriv) {
        if (points.size() < 2) {
            throw std::invalid_argument("Need at least 2 points");
        }
        if (values.size() != points.size()) {
            throw std::invalid_argument("Points and values must have same size");
        }
        
        size = points.size();
        setupKnots(points);
        computeCoefficients(points, values, startDeriv, endDeriv);
    }

private:
    void computeCoefficients(const std::vector<double>& points, 
                           const std::vector<double>& values,
                           double startDeriv, double endDeriv) {
        Matrix A(size + 2);
        ColVector b(size + 2);
        
        // Set the interpolation conditions
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 3; j++) {
                A[i][i+j] = getBasis(i+1+j, points[i]);
            }
            b[i] = values[i];
        }
        
        // Set the derivative boundary conditions
        for (int j = 0; j < 3; j++) {
            A[size][j] = getDerivative(1+j, points.front());
            A[size+1][size-1+j] = getDerivative(size+j, points.back());
        }
        b[size] = startDeriv;
        b[size+1] = endDeriv;
        
        Matrix result = Matrix::solve(A, b);
        coefs.clear();
        for (int i = 0; i <= size + 1; i++) {
            coefs.push_back(result[i][0]);
        }
    }
};

// Periodic Cubic B-Spline
class PeriodicCubicBSpline : public CubicBSpline {
public:
    PeriodicCubicBSpline(const std::vector<double>& points, const std::vector<double>& values) {
        if (points.size() < 2) {
            throw std::invalid_argument("Need at least 2 points");
        }
        if (values.size() != points.size()) {
            throw std::invalid_argument("Points and values must have same size");
        }
        if (std::abs(values.front() - values.back()) > 1e-10) {
            throw std::invalid_argument("First and last values must be equal for periodic spline");
        }
        
        size = points.size();
        setupKnots(points);
        computeCoefs(points, values);
    }

private:
    void computeCoefs(const std::vector<double>& points, const std::vector<double>& values) {
        Matrix A(size + 2);
        ColVector b(size + 2);
        
        // Set interpolation conditions
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 3; j++) {
                A[i][i+j] = getBasis(i+1+j, points[i]);
            }
            b[i] = values[i];
        }
        
        // Set periodic conditions for first derivatives
        for (int j = 0; j < 3; j++) {
            A[size][j] = getDerivative(1+j, points.front());
            A[size][size-1+j] = -getDerivative(size+j, points.back());
        }
        b[size] = 0;
        
        // Set periodic conditions for second derivatives
        for (int j = 0; j < 3; j++) {
            A[size+1][j] = getSecondDerivative(1+j, points.front());
            A[size+1][size-1+j] = -getSecondDerivative(size+j, points.back());
        }
        b[size+1] = 0;
        
        Matrix result = Matrix::solve(A, b);
        for (int i = 0; i <= size + 1; i++) {
            coefs.push_back(result[i][0]);
        }
    }
};


class CubicBSplineFactory {
public:
    // type: 1 - Natural, 2 - Clamped, 3 - Periodic
    static std::unique_ptr<CubicBSpline> createSpline(
        int type,
        const std::vector<double>& points,
        const std::vector<double>& values,
        double startDeriv = 0.0,
        double endDeriv = 0.0
    ) {
        switch (type) {
            case 1:  // Natural
                return std::make_unique<NaturalCubicBSpline>(points, values);
                
            case 2:  // Clamped
                return std::make_unique<ClampedCubicBSpline>(
                    points, 
                    values,
                    startDeriv,
                    endDeriv
                );
                
            case 3:  // Periodic
                return std::make_unique<PeriodicCubicBSpline>(points, values);
                
            default:
                throw std::invalid_argument("Invalid B-spline type. Use 1 for Natural, 2 for Clamped, or 3 for Periodic");
        }
    }
};


#endif // SPLINE_HPP