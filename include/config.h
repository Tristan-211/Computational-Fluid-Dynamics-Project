#pragma once
#include <vector>
#include <string>
#include <functional>
#include <array>

struct BoundaryCondition {
    std::string side;  // "top", "bottom", "left", "right"
    std::pair<double, double> range;  // e.g., {0.5, 1.25}
    std::string type;  // "inlet" or "outlet"
    std::function<double(double, double, double)> uFunction;  // u(x, y, t)
    std::function<double(double, double, double)> vFunction;  // v(x, y, t)
    std::function<double(double, double, double)> TFunction;  // T(x, y, t)
};

struct Rotor {
    std::array<double, 2> center;   // x, y
    double omega;
    double theta0;
    double temp;
    double width;
    double length;

    // Half-dimensions (cached)
    double halfWidth;
    double halfLength;

    // Precomputed target fields and transformed grids
    std::vector<std::vector<double>> IBvelu;
    std::vector<std::vector<double>> IBvelv;
    std::vector<std::vector<double>> IBtempT;

    std::vector<std::array<double, 2>> meshIBut_flat;
    std::vector<std::array<double, 2>> meshIBvt_flat;
    std::vector<std::array<double, 2>> meshIBTt_flat;

};


struct SolverConfig {
    int M, N;
    double Lx, Ly;
    double Re, Pr;
    double dt;
    std::vector<double> totalTime;
    double CFL;
    double h;
    std::vector<BoundaryCondition> boundaries;
    std::vector<Rotor> rotors;

    // solver type
    bool enableEllipticSolver = false;
    bool enableHyperbolicSolver = false;
};
