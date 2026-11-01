#pragma once

#include <functional>

// Custom force function
// Give the force function to be used in the MD simulation
// For example here are the Lennard-Jones potential and force functions

std::function<double(double)> ljForceFunction = [](double r) {
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    return 24.0 * (2.0 / r12 - 1.0 / r6) / r; 
};

std::function<double(double)> ljPotentialFunction = [](double r) {
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    return 4.0 * (1.0 / r12 - 1.0 / r6);
};