#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);

// Define F1: 1/x - tan(x)
class F1 : public Function {
public:
    double operator() (double x) const override {
        return 1.0/x - tan(x);
    }
};

// Define F2: x^{-1} - 2^x
class F2 : public Function {
public:
    double operator() (double x) const override {
        return 1 / x - pow(2, x);
    }
};

// Define F3: 2^{-x} + e^x + 2cos(x) - 6
class F3 : public Function {
public:
    double operator() (double x) const override {
        return pow(2, -x) + exp(x) + 2 * cos(x) - 6;
    }
};

// Define F4: (x^3 + 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2)
class F4 : public Function {
public:
    double operator() (double x) const override {
        return (pow(x, 3) + 4 * pow(x, 2) + 3 * x + 5) / (2 * pow(x, 3) - 9 * pow(x, 2) + 18 * x - 2);
    }
};

// Function to solve F1 using bisection
void solve_f1() {
    std::cout << "Solving 1/x - tan(x) on [0, Pi/2]" << std::endl;
    Bisection_Method solver_f1(F1(), 0, Pi/2);
    double x = solver_f1.solve();
    std::cout << "A root is: " << x << std::endl;
}

// Function to solve F2 using bisection
void solve_f2() {
    std::cout << "Solving 1/x - 2^x on [0, 1]" << std::endl;
    Bisection_Method solver_f2(F2(), 0, 1);
    double x = solver_f2.solve();
    std::cout << "A root is: " << x << std::endl;
}

// Function to solve F3 using bisection
void solve_f3() {
    std::cout << "Solving 2^{-x} + e^x + 2cos(x) - 6 on [1, 3]" << std::endl;
    Bisection_Method solver_f3(F3(), 1, 3);
    double x = solver_f3.solve();
    std::cout << "A root is: " << x << std::endl;
}

// Function to solve F4 using bisection
void solve_f4() {
    std::cout << "Solving (x^3 + 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2) on [0, 4]" << std::endl;
    Bisection_Method solver_f4(F4(), 0, 4);
    double x = solver_f4.solve();
    std::cout << "A root is: " << x << std::endl;
}

int main() {
    solve_f1();
    solve_f2();
    solve_f3();
    solve_f4();
    return 0;
}

