#pragma once
#include <vector>

struct Solution {
    std::vector<std::vector<double>> u;  // u-velocity
    std::vector<std::vector<double>> v;  // v-velocity
    std::vector<std::vector<double>> T;  // temperature
    std::vector<std::vector<double>> phi;// Phi
    std::vector<std::vector<double>> Qu, Qv, QT;

    // You can add more in the future, e.g.:
    // std::vector<std::vector<double>> p;  // pressure
    // std::vector<std::vector<double>> phi; // scalar field
};