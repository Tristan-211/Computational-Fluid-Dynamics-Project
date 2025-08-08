#include "IB.h"
#include "utils.h"
#include "config.h"
#include <sstream>
#include "solution.h"
#include <vector>
#include <string>
#include <cassert>


std::vector<std::vector<bool>> makeMask(const std::vector<std::vector<double>>& XR,
                                        const std::vector<std::vector<double>>& YR,
                                        double IBlh,
                                        double IBwh) {
    size_t M = XR.size();
    size_t N = XR[0].size();
    std::vector<std::vector<bool>> mask(M, std::vector<bool>(N, false));

    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            if (std::fabs(XR[i][j]) <= IBlh && std::fabs(YR[i][j]) < IBwh)
                mask[i][j] = true;

    return mask;
}


void applyMaskSet(std::vector<std::vector<double>>& target,
                  const std::vector<std::vector<bool>>& mask,
                  double value) {
    size_t M = target.size();
    size_t N = target[0].size();
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            if (mask[i][j])
                target[i][j] = value;
}


void checkRotorBounds(const std::vector<Rotor>& rotors, const Mesh& mesh) {
    double xmin = mesh.xf.front();
    double xmax = mesh.xf.back();
    double ymin = mesh.yf.front();
    double ymax = mesh.yf.back();

    for (size_t i = 0; i < rotors.size(); ++i) {
        const Rotor& r = rotors[i];

        double x0 = r.center[0];
        double y0 = r.center[1];

        double left   = x0 - r.halfLength;
        double right  = x0 + r.halfLength;
        double bottom = y0 - r.halfWidth;
        double top    = y0 + r.halfWidth;

        if (left < xmin || right > xmax || bottom < ymin || top > ymax) {
            std::ostringstream msg;
            msg << "Rotor " << i << " exceeds mesh bounds:\n"
                << "  Rotor bounds: ["
                << left << ", " << right << "] x ["
                << bottom << ", " << top << "]\n"
                << "  Mesh bounds: ["
                << xmin << ", " << xmax << "] x ["
                << ymin << ", " << ymax << "]";
            throw std::runtime_error(msg.str());
        }
    }
}


std::vector<std::vector<double>> reshape(const std::vector<double>& flat, size_t rows, size_t cols) {
    if (flat.size() != rows * cols) throw std::runtime_error("reshapeColumnMajor: size mismatch");
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));
    for (size_t j = 0; j < cols; ++j)
        for (size_t i = 0; i < rows; ++i)
            result[i][j] = flat[j * rows + i];  // column-major layout like MATLAB
    return result;
}



void initIBMeshU(const Mesh& mesh,SolverConfig& config,const  Solution& sol) {
    int M = config.M;
    int N = config.N;
    int numU = (M + 1) * (N + 2);
    int numV = (M + 2) * (N + 1);
    int rowsT = static_cast<int>(sol.T.size());
    if (rowsT == 0) throw std::runtime_error("initIBMeshT: sol.T is empty");
    int colsT = static_cast<int>(sol.T[0].size());
    // Choose coordinates based on T layout
    const std::vector<double>* X; // x centers
    const std::vector<double>* Y; // y centers
    if (rowsT == M + 2 && colsT == N + 2) {
        X = &mesh.xc;  Y = &mesh.yc;
    } else if (rowsT == M + 6 && colsT == N + 6) {
        X = &mesh.xc3; Y = &mesh.yc3;
    } else {
        throw std::runtime_error("initIBMeshT: Unexpected T dimensions");
    }
    int numT = rowsT * colsT;




    for (auto& r : config.rotors) {
        // Allocate storage
        r.IBvelu.assign(M + 1, std::vector<double>(N + 2, 0.0));
        r.meshIBut_flat.assign(2, std::vector<double>(numU, 0.0));

        //u flat mesh
        int k = 0;
        for (int j = 0; j < N + 2; ++j) {
            double y = mesh.yc[j];
            for (int i = 0; i < M + 1; ++i) {
                double x = mesh.xf[i];
                // IB velocity for this rotor (like meshyc - center)*omega
                r.IBvelu[i][j] = (y - r.center[1]) * r.omega;
                // Flattened mesh offsets
                r.meshIBut_flat[0][k] = x - r.center[0];
                r.meshIBut_flat[1][k] = y - r.center[1];
                ++k;
            }
        }
    }

    //v flast mesh
    for (auto& r : config.rotors) {
        r.IBvelv.assign(M + 2, std::vector<double>(N + 1, 0.0));
        r.meshIBvt_flat.assign(2, std::vector<double>(numV, 0.0));

        int k = 0; // column-major: j fast, i slow
        for (int j = 0; j < N + 1; ++j) {
            double y = mesh.yf[j];
            for (int i = 0; i < M + 2; ++i) {
                double x = mesh.xc[i];
                r.IBvelv[i][j] = -(x - r.center[0]) * r.omega; // depends only on i
                r.meshIBvt_flat[0][k] = x - r.center[0];
                r.meshIBvt_flat[1][k] = y - r.center[1];
                ++k;
            }
        }
    }

    for (auto& r : config.rotors) {
        r.IBtempT.assign(rowsT, std::vector<double>(colsT, r.temp)); // constant field per rotor
        r.meshIBTt_flat.assign(2, std::vector<double>(numT, 0.0));

        int k = 0; // column-major
        for (int j = 0; j < colsT; ++j) {
            double y = (*Y)[j];
            for (int i = 0; i < rowsT; ++i) {
                double x = (*X)[i];
                r.meshIBTt_flat[0][k] = x - r.center[0];
                r.meshIBTt_flat[1][k] = y - r.center[1];
                ++k;
            }
        }
    }


}


void calcSourceIB(Solution& sol,
                  const Mesh& mesh,
                  const SolverConfig& config,
                  double t)
{
    int M = config.M;
    int N = config.N;
    double dt = config.dt;

    // Prepare outputs with correct sizes
    sol.Qu.assign(M + 1, std::vector<double>(N + 2, 0.0));
    sol.Qv.assign(M + 2, std::vector<double>(N + 1, 0.0));
    int rowsT = static_cast<int>(sol.T.size());
    int colsT = rowsT ? static_cast<int>(sol.T[0].size()) : 0;
    sol.QT.assign(rowsT, std::vector<double>(colsT, 0.0));

    // Per-rotor accumulation
    for (const auto& r : config.rotors) {
        // Rotation matrix for this rotor
        double angle = r.theta0 + r.omega * t;
        std::vector<std::vector<double>> R = {
            { std::cos(angle), -std::sin(angle) },
            { std::sin(angle),  std::cos(angle) }
        };

        // ---- U component ----
        {
            int numU = (M + 1) * (N + 2);
            auto gridFlat = matMul(R, r.meshIBut_flat); // 2 x numU
            auto XR = reshape(gridFlat[0], M + 1, N + 2);
            auto YR = reshape(gridFlat[1], M + 1, N + 2);
            auto mask = makeMask(XR, YR, r.halfLength, r.halfWidth);

            for (int i = 0; i <= M; ++i) {
                for (int j = 0; j < N + 2; ++j) {
                    if (mask[i][j]) {
                        sol.Qu[i][j] += (r.IBvelu[i][j] - sol.u[i][j]) / dt;
                    }
                }
            }
        }

        // ---- V component ----
        {
            int numV = (M + 2) * (N + 1);
            auto gridFlat = matMul(R, r.meshIBvt_flat); // 2 x numV
            auto XR = reshape(gridFlat[0], M + 2, N + 1);
            auto YR = reshape(gridFlat[1], M + 2, N + 1);
            auto mask = makeMask(XR, YR, r.halfLength, r.halfWidth);

            for (int i = 0; i < M + 2; ++i) {
                for (int j = 0; j < N + 1; ++j) {
                    if (mask[i][j]) {
                        sol.Qv[i][j] += (r.IBvelv[i][j] - sol.v[i][j]) / dt;
                    }
                }
            }
        }

        // ---- T component ----
        if (rowsT && colsT) {
            auto gridFlat = matMul(R, r.meshIBTt_flat); // 2 x numT
            auto XR = reshape(gridFlat[0], rowsT, colsT);
            auto YR = reshape(gridFlat[1], rowsT, colsT);
            auto mask = makeMask(XR, YR, r.halfLength, r.halfWidth);

            for (int i = 0; i < rowsT; ++i) {
                for (int j = 0; j < colsT; ++j) {
                    if (mask[i][j]) {
                        sol.QT[i][j] += (r.IBtempT[i][j] - sol.T[i][j]) / dt;
                    }
                }
            }
        }
    }
}









