#include <iostream>
#include <vector>
#include <iomanip>
#include "Interpolation.hpp"

int main() {
    // Data for Day and Sp1
    std::vector<double> days = {0, 6, 10, 13, 17, 20, 28};
    std::vector<double> sp1_weights = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    std::vector<double> sp2_weights = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};

    // Create Newton Interpolation objects for Sp1 and Sp2
    NewtonInterpolation sp1_interpolation(days, sp1_weights);
    NewtonInterpolation sp2_interpolation(days, sp2_weights);

    // (a) Approximate the average weight curve
    std::cout << "Newton's Interpolating Polynomial Expression for Sp1: " << std::endl;
    std::cout << "f1" << "(x) = " << sp1_interpolation.getExpression() << std::endl;

    std::cout << "Newton's Interpolating Polynomial Expression for Sp2: " << std::endl;
    std::cout << "f2" << "(x) = " << sp2_interpolation.getExpression() << std::endl;

    // (b) Predict weight for the next 15 days
    double day_43 = 28 + 15;

    double sp1_prediction = sp1_interpolation.interpolate(day_43);
    double sp2_prediction = sp2_interpolation.interpolate(day_43);

    std::cout << "\nPrediction for Sp1 weight at day 43: " << sp1_prediction << std::endl;
    if (sp1_prediction < 0) {
        std::cout << "Sample1 of larvae will die after another 15 days." << std::endl;
    }
    std::cout << "Prediction for Sp2 weight at day 43: " << sp2_prediction << std::endl;
    if (sp2_prediction < 0) {
        std::cout << "Sample2 of larvae will die after another 15 days." << std::endl;
    }

    return 0;
}

