#include "multi_grid.h"
#include "bc.h"
#include <string>
#include <vector>
#include "utils.h"
#include <cmath>
#include <limits>
#include <iostream>


void gaussSeidel(
    std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<double>>& f,
    double h,
    int niter)
{
    int M = phi.size();       // number of rows (with ghost cells)
    int N = phi[0].size();    // number of cols (with ghost cells)
    double h2 = h * h;

    for (int iter = 0; iter < niter; ++iter) {
        for (int i = 1; i < M - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                phi[i][j] = 0.25 * (
                    phi[i-1][j] + phi[i+1][j] +
                    phi[i][j-1] + phi[i][j+1]
                ) - 0.25 * h2 * f[i][j];
            }
        }

        // Apply zero-Neumann boundary conditions to ghost cells
        bcGS(phi);
    }
}

void correctOutletMassConservation(
    Solution& sol,
    const Mesh& mesh,
    const SolverConfig& config)
{
    double h = mesh.h;
    double qds = 0.0;
    double outletLength = 0.0;

    // === Sum normal fluxes across all boundaries (inlets and outlets) ===
    for (const auto& bc : config.boundaries) {
        const std::string& side = bc.side;
        const auto& type = bc.type;

        if (side == "left") {
            int i = 0;
            for (int j = 0; j < static_cast<int>(mesh.yc.size()); ++j)
                qds += sol.u[i][j] * h;

        } else if (side == "right") {
            int i = static_cast<int>(mesh.xf.size()) - 1;
            for (int j = 0; j < static_cast<int>(mesh.yc.size()); ++j)
                qds -= sol.u[i][j] * h;

        } else if (side == "bottom") {
            int j = 0;
            for (int i = 0; i < static_cast<int>(mesh.xc.size()); ++i)
                qds += sol.v[i][j] * h;

        } else if (side == "top") {
            int j = static_cast<int>(mesh.yf.size()) - 1;
            for (int i = 0; i < static_cast<int>(mesh.xc.size()); ++i)
                qds -= sol.v[i][j] * h;
        }
    }

    // === Apply correction to outlet velocities only in outlet ranges ===
    for (const auto& bc : config.boundaries) {
        if (bc.type != "outlet") continue;

        const std::string& side = bc.side;
        double a = bc.range.first;
        double b = bc.range.second;

        if (side == "left" || side == "right") {
            const auto& yvec = mesh.yc;
            int i = (side == "left") ? 0 : static_cast<int>(mesh.xf.size()) - 1;

            for (int j = 0; j < static_cast<int>(yvec.size()); ++j) {
                double y = yvec[j];
                if (y >= a && y <= b) {
                    outletLength += h;
                }
            }

        } else if (side == "bottom" || side == "top") {
            const auto& xvec = mesh.xc;
            int j = (side == "bottom") ? 0 : static_cast<int>(mesh.yf.size()) - 1;

            for (int i = 0; i < static_cast<int>(xvec.size()); ++i) {
                double x = xvec[i];
                if (x >= a && x <= b) {
                    outletLength += h;
                }
            }
        }
    }

    if (outletLength == 0.0) return;

    double ucorr = qds / outletLength;
    
    // === Reapply correction to outlet points in range ===
    for (const auto& bc : config.boundaries) {
        if (bc.type != "outlet") continue;

        const std::string& side = bc.side;
        double a = bc.range.first;
        double b = bc.range.second;

        if (side == "left" || side == "right") {
            const auto& yvec = mesh.yc;
            int i = (side == "left") ? 0 : static_cast<int>(mesh.xf.size()) - 1;

            for (int j = 0; j < static_cast<int>(yvec.size()); ++j) {
                double y = yvec[j];
                if (y >= a && y <= b) {
                    sol.u[i][j] += ucorr;
                }
            }

        } else if (side == "bottom" || side == "top") {
            const auto& xvec = mesh.xc;
            int j = (side == "bottom") ? 0 : static_cast<int>(mesh.yf.size()) - 1;

            for (int i = 0; i < static_cast<int>(xvec.size()); ++i) {
                double x = xvec[i];
                if (x >= a && x <= b) {
                    sol.v[i][j] += ucorr;
                }
            }
        }
    }
}

std::vector<std::vector<double>> restrict(const std::vector<std::vector<double>>& rh) {
    int M = rh.size();       // total rows (including ghost cells)
    int N = rh[0].size();    // total cols (including ghost cells)

    int Mc = (M - 2) / 2 + 2;  // coarse size with 1-layer ghost cells
    int Nc = (N - 2) / 2 + 2;

    std::vector<std::vector<double>> r2h(Mc, std::vector<double>(Nc, 0.0));

    for (int i = 1; i < Mc - 1; ++i) {
        for (int j = 1; j < Nc - 1; ++j) {
            int fi = 2 * i - 1;  // fine-grid index aligned with MATLAB
            int fj = 2 * j - 1;

            r2h[i][j] = 0.25 * (
                rh[fi][fj] +
                rh[fi + 1][fj] +
                rh[fi][fj + 1] +
                rh[fi + 1][fj + 1]
            );
        }
    }

    return r2h;
}

std::vector<std::vector<double>> prolong(const std::vector<std::vector<double>>& e2h) {
    int Mc = e2h.size();       // coarse grid size with ghost cells
    int Nc = e2h[0].size();

    int Mf = (Mc - 2) * 2 + 2; // fine grid size with ghost cells
    int Nf = (Nc - 2) * 2 + 2;

    std::vector<std::vector<double>> eh(Mf, std::vector<double>(Nf, 0.0));

    for (int i = 1; i < Mc - 1; ++i) {
        for (int j = 1; j < Nc - 1; ++j) {
            double val = e2h[i][j];

            int fi = 2 * i - 1;  // fine-grid corresponding indices
            int fj = 2 * j - 1;

            eh[fi][fj] = val;
            eh[fi + 1][fj] = val;
            eh[fi][fj + 1] = val;
            eh[fi + 1][fj + 1] = val;
        }
    }

    bcGS(eh);  // Apply BCs to ghost cells after prolongation
    return eh;
}

std::vector<std::vector<double>> calcDivV(const Solution& sol, const Mesh& mesh) {
    const auto& u = sol.u;
    const auto& v = sol.v;
    double h = mesh.h;

    int Mu = u.size();    // u: Mu x Nu
    int Nu = u[0].size();

    std::vector<std::vector<double>> divV(Mu + 1, std::vector<double>(Nu, 0.0));

    for (int i = 1; i < Mu; ++i) {
        for (int j = 1; j < Nu - 1; ++j) {
            double du = (u[i][j] - u[i - 1][j]) / h;
            double dv = (v[i][j] - v[i][j - 1]) / h;
            divV[i][j] = du + dv;
        }
    }

    return divV;
}



std::vector<std::vector<double>> multigrid(
    std::vector<std::vector<double>> phi,
    const std::vector<std::vector<double>>& f,
    double h)
{
    int M = static_cast<int>(phi.size()) - 2;     // interior rows
    int N = static_cast<int>(phi[0].size()) - 2;  // interior cols

    int n = 1;   // iterations per level
    int nc = 3;  // coarse mesh extra iterations

    // Pre-smoothing
    for (int i = 0; i < n; ++i)
        gaussSeidel(phi, f, h, 1);  // in-place Gauss-Seidel + BC inside

    // V-cycle condition: mesh divisible by 2
    if (M % 2 == 0 && N % 2 == 0) {
        auto rh = Residual(phi, f, h);            // residual on fine grid
        auto r2h = restrict(rh);                 // restrict to coarse grid

        int Mc = M / 2 + 2;
        int Nc = N / 2 + 2;
        std::vector<std::vector<double>> e2h(Mc, std::vector<double>(Nc, 0.0));

        // Recursive call on coarse mesh
        e2h = multigrid(e2h, r2h, 2.0 * h);

        // Prolong and correct fine-grid solution
        auto eh = prolong(e2h);
        for (size_t i = 0; i < phi.size(); ++i)
            for (size_t j = 0; j < phi[0].size(); ++j)
                phi[i][j] += eh[i][j];

        bcGS(phi);  // enforce BCs after correction

        // Post-smoothing
        for (int i = 0; i < n; ++i)
            gaussSeidel(phi, f, h, 1);

    } else {
        // On coarsest grid: perform extra Gauss-Seidel
        for (int i = 0; i < nc; ++i)
            gaussSeidel(phi, f, h, 1);
    }

    return phi;
}

std::vector<std::vector<double>> poissonSolver(
    std::vector<std::vector<double>> phi,
    const std::vector<std::vector<double>>& f,
    const Mesh& mesh,
    int nIterMax,
    double epsilon)
{
    double rr = std::numeric_limits<double>::max();
    int iter = 0;
    double h = mesh.h;

    while (iter < nIterMax && rr > epsilon) {
        iter += 1;
        phi = multigrid(phi, f, h);
        rr = RelativeResidualNorm(phi, f, h);
    }

    if (rr > epsilon) {
        std::cerr << "Warning: Poisson solver did not converge.\n";
        std::cerr << "Final residual = " << rr << " after " << iter << " iterations.\n";
    }

    return phi;
}


void projectV(Solution& sol, const Mesh& mesh, const SolverConfig& config, double t)
{
    double h = mesh.h;
    double dt = config.dt;
    const auto& phi = sol.phi;

    // Compute dphi/dx on u-staggered grid
    std::vector<std::vector<double>> dphi_dx(config.M + 1, std::vector<double>(config.N + 2, 0.0));

    for (int i = 0; i < config.M + 1; ++i) {
        for (int j = 0; j < config.N + 2; ++j) {
            dphi_dx[i][j] = (phi[i + 1][j] - phi[i][j]) / h;
        }
    }

    // Project u: interior
    for (int i = 1; i < config.M; ++i) {           // skip bottom and top ghost rows
        for (int j = 1; j < config.N + 1; ++j) {   // skip left and right ghost cols
            sol.u[i][j] -= dt * dphi_dx[i][j];
        }
    }

    bcGhost_u(sol.u, mesh, config, t + dt);


    std::vector<std::vector<double>> dphi_dy(config.M + 2, std::vector<double>(config.N, 0.0));

    for (int i = 0; i <= config.M + 1; ++i) {
        for (int j = 0; j < config.N; ++j) {
            dphi_dy[i][j] = (phi[i][j + 1] - phi[i][j]) / h;
        }
    }

    // Project v: interior
    for (int i = 1; i <= config.M; ++i) {
        for (int j = 1; j < config.N; ++j) {
            sol.v[i][j] -= dt * dphi_dy[i][j];
        }
    }

    bcGhost_v(sol.v, mesh, config, t + dt);
}