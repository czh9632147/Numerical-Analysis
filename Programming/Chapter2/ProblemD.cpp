#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "Interpolation.hpp"

// Assume HermiteInterpolation class is included here.

int main() {
    // Time values (in seconds)
    std::vector<double> time = {0, 3, 5, 8, 13};

    // Displacement values (in feet)
    std::vector<double> displacement = {0, 225, 383, 623, 993};

    // Velocity values (in feet per second), which are the derivatives f'(x)
    std::vector<double> velocity = {75, 77, 80, 74, 72};

    // Create an instance of HermiteInterpolation
    HermiteInterpolation hermite(time, displacement, velocity);

    // Part (a) - Use Hermite polynomial to predict position and speed at t = 10
    double t_predict = 10.0;
    double predicted_position = hermite.interpolate(t_predict);
    double predicted_velocity = abs(hermite.derivativeofx(t_predict));

    std::cout << std::fixed << std::setprecision(6);  // Set precision for output
    std::cout << "Part (a):" << std::endl;
    std::cout << "Predicted position of the car at t = " << t_predict << " seconds: " 
              << predicted_position << " feet" << std::endl;
    std::cout << "Predicted velocity of the car at t = " << t_predict << " seconds: " 
              << predicted_velocity << " feet/second" << std::endl;

    // Part (b) - Check if the car exceeds 55 miles per hour (i.e., 81 feet per second)
    const double speed_limit = 81.0;  // Speed limit in feet per second

    bool exceeds_speed_limit = false;
    for (double t = 0; t <= 13; t += 0.01) {  // Finer step for accuracy
        double speed = abs(hermite.derivativeofx(t));
        if (speed > speed_limit) {
            exceeds_speed_limit = true;
            std::cout << "Part (b): The car exceeds the speed limit of 81 ft/sec at t = " 
                      << std::fixed << std::setprecision(2) << t << " seconds." << std::endl;
            break;
        }
    }

    if (!exceeds_speed_limit) {
        std::cout << "Part (b): The car never exceeds the speed limit of 81 ft/sec." << std::endl;
    }

    return 0;
}

