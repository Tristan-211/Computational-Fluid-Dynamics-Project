#include <vector>
#include <cmath>
#include <utility> // for std::pair
#include "config.h"
#include "solution.h"
#include "hyperbolic.h"
#include "bc.h"
#include "utils.h"

using std::vector;
using Array2D = std::vector<std::vector<double>>;

std::pair<Array2D, Array2D> hyperbolic_uv_2D(
    Array2D& u,
    Array2D& v,
    const SolverConfig& config)
{
    double h = config.h;

    int Mu = static_cast<int>(u.size());
    int Nu = static_cast<int>(u[0].size());
    int Mv = static_cast<int>(v.size());
    int Nv = static_cast<int>(v[0].size());

    // === Interpolated velocities ===
    Array2D ucen(Mu - 1, std::vector<double>(Nu));
    Array2D ucor(Mu, std::vector<double>(Nu - 1));
    Array2D vcen(Mv, std::vector<double>(Nv - 1));
    Array2D vcor(Mv - 1, std::vector<double>(Nv));

    for (int i = 0; i < Mu - 1; ++i)
        for (int j = 0; j < Nu; ++j)
            ucen[i][j] = 0.5 * (u[i][j] + u[i + 1][j]);

    for (int i = 0; i < Mu; ++i)
        for (int j = 0; j < Nu - 1; ++j)
            ucor[i][j] = 0.5 * (u[i][j] + u[i][j + 1]);

    for (int i = 0; i < Mv; ++i)
        for (int j = 0; j < Nv - 1; ++j)
            vcen[i][j] = 0.5 * (v[i][j] + v[i][j + 1]);

    for (int i = 0; i < Mv - 1; ++i)
        for (int j = 0; j < Nv; ++j)
            vcor[i][j] = 0.5 * (v[i][j] + v[i + 1][j]);

    // === Compute duu and duvx ===
    Array2D duu(Mu - 2, std::vector<double>(Nu - 2));
    for (int i = 0; i < Mu - 2; ++i)
        for (int j = 1; j < Nu - 1; ++j)
            duu[i][j - 1] = (std::pow(ucen[i + 1][j], 2) - std::pow(ucen[i][j], 2)) / h;

    Array2D duvx(Mu - 2, std::vector<double>(Nu - 2));
    for (int i = 0; i < Mu - 2; ++i)
        for (int j = 0; j < Nu - 2; ++j)
            duvx[i][j] = (ucor[i + 1][j + 1] * vcor[i + 1][j + 1] - ucor[i + 1][j] * vcor[i + 1][j]) / h;

    // === Hu: size (Mu, Nu), interior only updated ===
    Array2D Hu(Mu, std::vector<double>(Nu, 0.0));
    for (int i = 1; i < Mu - 1; ++i)
        for (int j = 1; j < Nu - 1; ++j)
            Hu[i][j] = -duu[i - 1][j - 1] - duvx[i - 1][j - 1];

    // === Compute dvv and duvy ===
    Array2D dvv(Mv - 2, std::vector<double>(Nv - 2));
    for (int i = 1; i < Mv - 1; ++i)
        for (int j = 0; j < Nv - 2; ++j)
            dvv[i - 1][j] = (std::pow(vcen[i][j + 1], 2) - std::pow(vcen[i][j], 2)) / h;

    Array2D duvy(Mv - 2, std::vector<double>(Nv - 2));
    for (int i = 0; i < Mv - 2; ++i)
        for (int j = 0; j < Nv - 2; ++j)
            duvy[i][j] = (ucor[i + 1][j + 1] * vcor[i + 1][j + 1] - ucor[i][j + 1] * vcor[i][j + 1]) / h;

    // === Hv: size (Mv, Nv), interior only updated ===
    Array2D Hv(Mv, std::vector<double>(Nv, 0.0));
    for (int i = 1; i < Mv - 1; ++i) 
        for (int j = 1; j < Nv - 1; ++j) 
            Hv[i][j] = -duvy[i - 1][j - 1] - dvv[i - 1][j - 1];

    return {Hu, Hv};
}


double psiWENO(double a, double b, double c, double d) {
    constexpr double eps = 1e-6;

    // Smoothness indicators
    double IS0 = 13.0 * std::pow(a - b, 2) + 3.0 * std::pow(a - 3.0 * b, 2);
    double IS1 = 13.0 * std::pow(b - c, 2) + 3.0 * std::pow(b + c, 2);
    double IS2 = 13.0 * std::pow(c - d, 2) + 3.0 * std::pow(3.0 * c - d, 2);

    // Alpha weights
    double a0 = 1.0 / std::pow(eps + IS0, 2);
    double a1 = 6.0 / std::pow(eps + IS1, 2);
    double a2 = 3.0 / std::pow(eps + IS2, 2);

    // Nonlinear weights
    double wSum = a0 + a1 + a2;
    double w0 = a0 / wSum;
    double w2 = a2 / wSum;

    // Final WENO psi value
    return (1.0 / 3.0) * w0 * (a - 2.0 * b + c) + (1.0 / 6.0) * (w2 - 0.5) * (b - 2.0 * c + d);
}



Array2D calc_adTdx_WENO_2D(const Array2D& T, const Array2D& a, const SolverConfig& config) {
    const int M = static_cast<int>(T.size());
    const int N = static_cast<int>(T[0].size());
    double h = config.h;

    Array2D f(M, std::vector<double>(N, 0.0));

    for (int i = 3; i < M - 3; ++i) {
        for (int j = 3; j < N - 3; ++j) {
            if (a[i][j] > 0) {
                double ap = (T[i - 1][j] - 2 * T[i - 2][j] + T[i - 3][j]) / h;
                double bp = (T[i][j]     - 2 * T[i - 1][j] + T[i - 2][j]) / h;
                double cp = (T[i + 1][j] - 2 * T[i][j]     + T[i - 1][j]) / h;
                double dp = (T[i + 2][j] - 2 * T[i + 1][j] + T[i][j]) / h;
                double psi = psiWENO(ap, bp, cp, dp);

                f[i][j] = ((T[i - 2][j] - 8 * T[i - 1][j] + 8 * T[i + 1][j] - T[i + 2][j]) / (12 * h) - psi) * a[i][j];
            } else if (a[i][j] < 0) {
                double ap = (T[i + 3][j] - 2 * T[i + 2][j] + T[i + 1][j]) / h;
                double bp = (T[i + 2][j] - 2 * T[i + 1][j] + T[i][j]) / h;
                double cp = (T[i + 1][j] - 2 * T[i][j]     + T[i - 1][j]) / h;
                double dp = (T[i][j]     - 2 * T[i - 1][j] + T[i - 2][j]) / h;
                double psi = psiWENO(ap, bp, cp, dp);

                f[i][j] = ((T[i - 2][j] - 8 * T[i - 1][j] + 8 * T[i + 1][j] - T[i + 2][j]) / (12 * h) + psi) * a[i][j];
            }
        }
    }

    return f;
}

Array2D calc_bdTdy_WENO_2D(const Array2D& T, const Array2D& b, const SolverConfig& config) {
    const int M = static_cast<int>(T.size());
    const int N = static_cast<int>(T[0].size());
    double h = config.h;

    Array2D f(M, std::vector<double>(N, 0.0));

    for (int i = 3; i < M - 3; ++i) {
        for (int j = 3; j < N - 3; ++j) {
            if (b[i][j] > 0) {
                double ap = (T[i][j - 1] - 2 * T[i][j - 2] + T[i][j - 3]) / h;
                double bp = (T[i][j]     - 2 * T[i][j - 1] + T[i][j - 2]) / h;
                double cp = (T[i][j + 1] - 2 * T[i][j]     + T[i][j - 1]) / h;
                double dp = (T[i][j + 2] - 2 * T[i][j + 1] + T[i][j]) / h;
                double psi = psiWENO(ap, bp, cp, dp);

                f[i][j] = ((T[i][j - 2] - 8 * T[i][j - 1] + 8 * T[i][j + 1] - T[i][j + 2]) / (12 * h) - psi) * b[i][j];
            } else if (b[i][j] < 0) {
                double ap = (T[i][j + 3] - 2 * T[i][j + 2] + T[i][j + 1]) / h;
                double bp = (T[i][j + 2] - 2 * T[i][j + 1] + T[i][j]) / h;
                double cp = (T[i][j + 1] - 2 * T[i][j]     + T[i][j - 1]) / h;
                double dp = (T[i][j]     - 2 * T[i][j - 1] + T[i][j - 2]) / h;
                double psi = psiWENO(ap, bp, cp, dp);

                f[i][j] = ((T[i][j - 2] - 8 * T[i][j - 1] + 8 * T[i][j + 1] - T[i][j + 2]) / (12 * h) + psi) * b[i][j];
            }
        }
    }

    return f;
}


Array2D hyperbolic_T_WENO_2D(const Solution& sol, const Array2D& u_old, const Array2D& v_old, double t,const SolverConfig& config,const Mesh& mesh) {
    const Array2D& T = sol.T;
    auto dt = config.dt;
    int M = T.size();
    int N = T[0].size();

    int Mu = sol.u.size();
    int Nu = sol.u[0].size();
    int Mv = sol.v.size();
    int Nv = sol.v[0].size();

    // Cell-centered u velocity (ucen0,1,2)
    Array2D ucen0(Mu - 1, std::vector<double>(Nu));
    Array2D ucen1(Mu - 1, std::vector<double>(Nu));
    Array2D ucen2(Mu - 1, std::vector<double>(Nu));
    for (int i = 0; i < Mu - 1; ++i)
        for (int j = 0; j < Nu; ++j) {
            ucen0[i][j] = 0.5 * (u_old[i][j] + u_old[i + 1][j]);
            ucen1[i][j] = 0.5 * (sol.u[i][j] + sol.u[i + 1][j]);
            ucen2[i][j] = 0.5 * (ucen0[i][j] + ucen1[i][j]);
        }

    // Cell-centered v velocity (vcen0,1,2)
    Array2D vcen0(Mv, std::vector<double>(Nv - 1));
    Array2D vcen1(Mv, std::vector<double>(Nv - 1));
    Array2D vcen2(Mv, std::vector<double>(Nv - 1));
    for (int i = 0; i < Mv; ++i)
        for (int j = 0; j < Nv - 1; ++j) {
            vcen0[i][j] = 0.5 * (v_old[i][j] + v_old[i][j + 1]);
            vcen1[i][j] = 0.5 * (sol.v[i][j] + sol.v[i][j + 1]);
            vcen2[i][j] = 0.5 * (vcen0[i][j] + vcen1[i][j]);
        }

    // Allocate a*, b* velocity fields for WENO
    Array2D a0(M, std::vector<double>(N, 0.0)), a1 = a0, a2 = a0;
    Array2D b0 = a0, b1 = a0, b2 = a0;

    for (int i = 3; i < M - 3; ++i) {
        for (int j = 3; j < N - 3; ++j) {
            a0[i][j] = ucen0[i - 3][j - 2];
            a1[i][j] = ucen1[i - 3][j - 2];
            a2[i][j] = ucen2[i - 3][j - 2];

            b0[i][j] = vcen0[i - 2][j - 3];
            b1[i][j] = vcen1[i - 2][j - 3];
            b2[i][j] = vcen2[i - 2][j - 3];
        }
    }

    // Runge-Kutta stage 1
    Array2D f0x = calc_adTdx_WENO_2D(T, a0, config);
    Array2D f0y = calc_bdTdy_WENO_2D(T, b0, config);
    Array2D f0 = elementwiseAdd(f0x, f0y);
    Array2D T1 = elementwiseSubtract(T, scalarMultiply(f0, dt));
    bc_T(T1, mesh, config, t + dt);


    //saveMatrixToFile(T1, "T1.csv");
    //saveMatrixToFile(b0, "b0.csv");
    //saveMatrixToFile(f0x, "f0x.csv");
    //saveMatrixToFile(f0y, "f0y.csv");


    // Stage 2
    Array2D f1x = calc_adTdx_WENO_2D(T1, a1, config);
    Array2D f1y = calc_bdTdy_WENO_2D(T1, b1, config);
    Array2D f1 = elementwiseAdd(f1x, f1y);
    Array2D T2 = elementwiseAdd(T1,
                     elementwiseSubtract(
                         scalarMultiply(f0, 0.75 * dt),
                         scalarMultiply(f1, 0.25 * dt)));
    bc_T(T2, mesh, config, t + dt / 2.0);


    //saveMatrixToFile(T2, "T2.csv");

    // Stage 3
    Array2D f2x = calc_adTdx_WENO_2D(T2, a2, config);
    Array2D f2y = calc_bdTdy_WENO_2D(T2, b2, config);
    Array2D f2 = elementwiseAdd(f2x, f2y);
    Array2D T3 = elementwiseAdd(
                    T2,
                    elementwiseAdd(
                        scalarMultiply(f0, dt / 12.0),
                        elementwiseAdd(
                            scalarMultiply(f1, dt / 12.0),
                            scalarMultiply(f2, -2.0 * dt / 3.0))));
    bc_T(T3, mesh, config, t + dt);


    //saveMatrixToFile(T3, "Ts.csv");

    // Final output
    return scalarMultiply(elementwiseSubtract(T3, T), 1.0/config.dt);
}