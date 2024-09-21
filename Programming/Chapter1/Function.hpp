#ifndef FUNCTION
#define FUNCTION

#include <cmath>

class Function {
public:
    // Pure virtual function for evaluating the function at x
    virtual double operator() (double x) const = 0;
    
    // Default implementation for numerical derivative using central difference method
    virtual double derivative(double x) const {
        double h = 1e-7; // Small step size for finite differences
        return ( (*this)(x + h) - (*this)(x - h) ) / (2 * h);
    }
};

class TanEquation : public Function {
public:
    virtual double operator() (double x) const override {
        return x - std::tan(x);
    }

    // Define f'(x) = 1 - sec^2(x) = 1 - 1/cos^2(x)
    virtual double derivative(double x) const override {
        return 1 - 1.0 / (std::cos(x) * std::cos(x));
    }
};

// sin(x/2) - 1
class SinEquation : public Function {
public:
    virtual double operator() (double x) const override {
        return std::sin(x / 2) - 1;
    }
};

// e^x - tan(x)
class ExpTanEquation : public Function {
public:
    virtual double operator() (double x) const override {
        return std::exp(x) - std::tan(x);
    }
};

// x^3 - 12x^2 + 3x + 1
class PolyEquation : public Function {
public:
    virtual double operator() (double x) const override {
        return std::pow(x, 3) - 12 * std::pow(x, 2) + 3 * x + 1;
    }
};

// Function for the volume of water in the trough
class TroughVolume : public Function {
private:
    double L, r, V;

public:
    TroughVolume(double L, double r, double V) : L(L), r(r), V(V) {}

    virtual double operator() (double h) const override {
        return L * (0.5 * M_PI * r * r - r * r * std::asin(h / r) - h * std::sqrt(r * r - h * h)) - V;
    }

    // Derivative of the function for Newton's method
    virtual double derivative(double h) const override {
        return L * (-r * std::sqrt(r * r - h * h) - std::sqrt(r * r - h * h) + h * h / std::sqrt(r * r - h * h));
    }
};

// ProblemF
class AlphaEquation : public Function {
private:
    double A, B, C, E;

public:
    AlphaEquation(double l, double h, double D, double beta1_deg) {
        double beta1 = beta1_deg * M_PI / 180.0;  // Convert degrees to radians
        A = l * std::sin(beta1);
        B = l * std::cos(beta1);
        C = (h + 0.5 * D) * std::sin(beta1) - 0.5 * D * std::tan(beta1);
        E = (h + 0.5 * D) * std::cos(beta1) - 0.5 * D;
    }

    virtual double operator() (double alpha_deg) const override {
        double alpha = alpha_deg * M_PI / 180.0;  // Convert degrees to radians
        return A * std::sin(alpha) * std::cos(alpha) + B * std::sin(alpha) * std::sin(alpha)
            - C * std::cos(alpha) - E * std::sin(alpha);
    }

    virtual double derivative(double alpha_deg) const override {
        double alpha = alpha_deg * M_PI / 180.0;
        return A * std::cos(2 * alpha) + 2 * B * std::sin(alpha) * std::cos(alpha)
            + C * std::sin(alpha) - E * std::cos(alpha);
    }
};

#endif

