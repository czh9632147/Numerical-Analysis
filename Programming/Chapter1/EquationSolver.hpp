#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"
#include <cmath>
#include <iostream>

class EquationSolver {
protected:
    const Function &F;
public:
    EquationSolver(const Function &F) : F(F) {}
    virtual double solve() = 0;
};

class Bisection_Method : public EquationSolver {
private:
    double a, b;
    double eps, delta;
    int Maxiter;

public:
    Bisection_Method(const Function &F, double a, double b, 
        double eps = 1e-7, double delta = 1e-6, int Maxiter = 50) :
        EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}
    
    virtual double solve() override {
        double fa = F(a), fb = F(b);

        // Ensure that there is a sign change between a and b
        if (fa * fb > 0) {
            std::cerr << "Bisection method error: no sign change on the interval." << std::endl;
            return NAN;
        }

        for (int iter = 0; iter < Maxiter; ++iter) {
            double c = (a + b) / 2;
            double fc = F(c);

            // Check if the result is within the tolerance limits
            if (std::fabs(fc) < delta || (b - a) / 2 < eps) {
                return c;
            }

            // Narrow the interval based on the sign of the function at midpoint
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }
        }

        // Return the midpoint after reaching the maximum number of iterations
        return (a + b) / 2;
    }
};

class Newton_Method : public EquationSolver {
private:
    double x0;
    double eps;
    int Maxiter;

public:
    Newton_Method(const Function &F, double x0, 
        double eps = 1e-7, int Maxiter = 8) :
        EquationSolver(F), x0(x0), eps(eps), Maxiter(Maxiter) {}
    
    virtual double solve() override {
        double x = x0;

        for (int iter = 0; iter < Maxiter; ++iter) {
            double fx = F(x);
            double fprime_x = F.derivative(x); // Assuming `F` has a `derivative` method
            
            if (std::fabs(fprime_x) < eps) {
                std::cerr << "Newton's method error: derivative too small." << std::endl;
                return NAN;
            }

            double x_new = x - fx / fprime_x;

            if (std::fabs(x_new - x) < eps) {
                return x_new; // Root found
            }

            x = x_new;
        }

        return x; // Return the last computed x value
    }
};


class Secant_Method : public EquationSolver {
private:
    double x0, x1;
    double eps;
    int Maxiter;

public:
    Secant_Method(const Function &F, double x0, double x1, 
                  double eps = 1e-7, int Maxiter = 50) :
        EquationSolver(F), x0(x0), x1(x1), eps(eps), Maxiter(Maxiter) {}
    
    virtual double solve() override {
        double x_prev = x0;
        double x_curr = x1;
        for (int iter = 0; iter < Maxiter; ++iter) {
            double f_curr = F(x_curr);
            double f_prev = F(x_prev);
            
            if (std::fabs(f_curr) < eps) {
                return x_curr;  // Found the root
            }
            
            // Secant method update
            double denominator = f_curr - f_prev;
            if (std::fabs(denominator) < eps) {
                std::cerr << "Secant method error: denominator too small." << std::endl;
                return NAN;
            }
            
            double x_next = x_curr - f_curr * (x_curr - x_prev) / denominator;
            
            if (std::fabs(x_next - x_curr) < eps) {
                return x_next;  // Convergence criterion met
            }
            
            x_prev = x_curr;
            x_curr = x_next;
        }
        
        std::cerr << "Secant method error: did not converge within maximum iterations." << std::endl;
        return NAN;
    }
};

#endif

