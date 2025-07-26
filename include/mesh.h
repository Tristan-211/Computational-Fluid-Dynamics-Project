#pragma once
#include <vector>

struct Mesh {
    std::vector<double> xf;  // x-face positions (u-velocity)
    std::vector<double> yc;  // y-center positions (u-velocity)
    std::vector<double> yf;  // y-face positions (v-velocity)
    std::vector<double> xc;  // x-center positions (v-velocity)
    std::vector<double> xc3; // x-center positions 3 ghost cells (u-velocity)
    std::vector<double> yc3; // y-center positions 3 ghost cells (u-velocity)
    double h;                // grid spacing
};

Mesh generateMesh(double Lx, double Ly, int M, int N);