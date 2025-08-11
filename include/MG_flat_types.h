#pragma once
#include <vector>
#include <algorithm>
#include <cmath>

struct MGLevel {
    int nx{}, ny{};   // interior sizes (no ghosts)
    double h{};
    // All vectors sized to (nx+2)*(ny+2) including 1 ghost cell each side.
    std::vector<double> phi, f, res, tmp, err;
};

struct MG {
    std::vector<MGLevel> L; // 0 = finest
};

// Row-major idx; i=0..nx+1, j=0..ny+1 (ghosts included)
inline int idx(int i, int j, int nx) { return i + (nx + 2) * j; }
