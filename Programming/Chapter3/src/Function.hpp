#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <vector>
#include <stdexcept>
#include <cmath>
#include <sstream>

class Function {
protected:
    std::vector<double> points;
    std::vector<double> values;
    size_t size;
public:
    virtual double eval(const double& x) = 0;  // operator() renamed to eval
    
    // First derivative
    double deriv(const double& x) {
        const double eps = 1e-6;
        return (eval(x + eps) - eval(x - eps)) / (2 * eps);
    }
    
    // Second derivative
    double deriv2(const double& x) {
        const double eps = 1e-6;
        return (eval(x + 2*eps) + eval(x - 2*eps) - 2*eval(x)) / (4 * eps * eps);
    }

    virtual ~Function() = default;
};


class DiscreteFunction : public Function {
private:
    std::vector<double> derivs;
    std::vector<double> derivs2;

public:
    DiscreteFunction(const std::vector<double>& pts, 
                     const std::vector<double>& vals,
                     const std::vector<double>& ders, 
                     const std::vector<double>& ders2) {
        if (pts.size() != vals.size() || 
            vals.size() != ders.size() || 
            ders.size() != ders2.size()) {
            throw std::invalid_argument("Vectors must have equal size");
        }
        
        points = pts;
        values = vals;
        derivs = ders;
        derivs2 = ders2;
        size = points.size();
    }

    // Evaluate function at x
    double eval(const double& x) override {
        for (size_t i = 0; i < size; i++) {
            if (std::abs(points[i] - x) < 1e-10) {
                return values[i];
            }
        }
        throw std::invalid_argument("Point not found in discrete set");
    }
    
    // Get first derivative at x
    double deriv(const double& x) {
        for (size_t i = 0; i < size; i++) {
            if (std::abs(points[i] - x) < 1e-10) {
                return derivs[i];
            }
        }
        throw std::invalid_argument("Point not found in discrete set");
    }
    
    // Get second derivative at x
    double deriv2(const double& x) {
        for (size_t i = 0; i < size; i++) {
            if (std::abs(points[i] - x) < 1e-10) {
                return derivs2[i];
            }
        }
        throw std::invalid_argument("Point not found in discrete set");
    }
};


class Polynomial : public Function {
private:
    std::vector<double> coeffs;
    size_t degree;

public:
    Polynomial() : degree(0) {}
    
    Polynomial(const std::vector<double>& coefficients) 
        : coeffs(coefficients), degree(coefficients.size()) {}
    
    // Evaluate polynomial at x
    double eval(const double& x) override {
        double result = 0.0;
        double power = 1.0;
        
        for (size_t i = 0; i < degree; i++) {
            result += coeffs[i] * power;
            power *= x;
        }
        return result;
    }
    
    // Get derivative polynomial
    Polynomial getDerivative() const {
        if (degree <= 1) {
            return Polynomial({0.0});
        }
        
        std::vector<double> deriv_coeffs;
        for (size_t i = 1; i < degree; i++) {
            deriv_coeffs.push_back(coeffs[i] * i);
        }
        return Polynomial(deriv_coeffs);
    }
    
    // Get coefficients
    std::vector<double> getCoeffs() const {
        return coeffs;
    }
    
    // Get degree
    size_t getDegree() const {
        return degree;
    }
    
    // Convert polynomial to string
    std::string toString() const {
        if (degree == 0) return "0";
        
        std::stringstream ss;
        bool first = true;
        
        for (int i = degree - 1; i >= 0; i--) {
            if (coeffs[i] == 0) continue;
            
            if (coeffs[i] > 0 && !first) {
                ss << "+";
            }
            
            if (i == 0) {
                ss << coeffs[i];
            } else if (i == 1) {
                if (coeffs[i] == 1) ss << "x";
                else if (coeffs[i] == -1) ss << "-x";
                else ss << coeffs[i] << "x";
            } else {
                if (coeffs[i] == 1) ss << "x^" << i;
                else if (coeffs[i] == -1) ss << "-x^" << i;
                else ss << coeffs[i] << "x^" << i;
            }
            
            first = false;
        }
        
        return ss.str();
    }
};

#endif