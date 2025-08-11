#include "MG_flat.h"
#include "bc.h"
#include "MG_flat_types.h"
#include <vector>
#include "config.h"
#include "solution.h"
#include "mesh.h"
#include <cmath>



void residual_flat(MGLevel& lv) {
    const int nx=lv.nx, ny=lv.ny; const double invh2=1.0/(lv.h*lv.h);
    for (int j=1;j<=ny;++j)
        for (int i=1;i<=nx;++i) {
            int p = idx(i,j,nx);
            double lap = (lv.phi[idx(i-1,j,nx)] + lv.phi[idx(i+1,j,nx)]
                        + lv.phi[idx(i,j-1,nx)] + lv.phi[idx(i,j+1,nx)]
                        - 4.0*lv.phi[p]) * invh2;
            lv.res[p] = lv.f[p] - lap;
        }
}

void mg_build_auto(MG& mg, int nx0, int ny0, double h0) {
    mg.L.clear();
    int nx = nx0, ny = ny0; double h = h0;
    while (true) {
        MGLevel lv;
        lv.nx = nx; lv.ny = ny; lv.h = h;
        const int sz = (nx + 2) * (ny + 2);
        lv.phi.assign(sz, 0.0);
        lv.f  .assign(sz, 0.0);
        lv.res.assign(sz, 0.0);
        mg.L.push_back(std::move(lv));

        if ((nx % 2) != 0 || (ny % 2) != 0) break;
        if (nx < 2 || ny < 2) break;
        nx /= 2; ny /= 2; h *= 2.0;
    }
}


double maxF_inf_flat(const std::vector<double>& f, int nx, int ny) {
    double m = 0.0;
    for (int j = 1; j <= ny; ++j)
        for (int i = 1; i <= nx; ++i)
            m = std::max(m, std::abs(f[idx(i,j,nx)]));
    return m;
}


double RelativeResidualNorm_from_buffer(const std::vector<double>& res,
                                        const std::vector<double>& f,
                                        int nx, int ny, double max_f /*precomputed or 0*/)
{
    double max_r = 0.0;
    if (max_f <= 0.0) {
        // compute both max_r and max_f here (slightly more work)
        for (int j = 1; j <= ny; ++j)
            for (int i = 1; i <= nx; ++i) {
                const int p = idx(i,j,nx);
                max_r = std::max(max_r, std::abs(res[p]));
                max_f = std::max(max_f, std::abs(f[p]));
            }
    } else {
        // only max_r; max_f provided
        for (int j = 1; j <= ny; ++j)
            for (int i = 1; i <= nx; ++i) {
                const int p = idx(i,j,nx);
                max_r = std::max(max_r, std::abs(res[p]));
            }
    }
    return (max_f > 0.0) ? (max_r / max_f) : 0.0;
}


void restrict_flat(const std::vector<double>& rh_fine,
                           int nx_f, int ny_f,
                           std::vector<double>& r2h_coarse,
                           int nx_c, int ny_c)
{
    // Coarse total size must already be (nx_c+2)*(ny_c+2)
    // Loop coarse interior: I=1..nx_c, J=1..ny_c
    for (int J = 1; J <= ny_c; ++J) {
        const int fj = 2*J - 1;
        for (int I = 1; I <= nx_c; ++I) {
            const int fi = 2*I - 1;
            const double v =
                0.25 * ( rh_fine[idx(fi,    fj,    nx_f)]
                       + rh_fine[idx(fi+1,  fj,    nx_f)]
                       + rh_fine[idx(fi,    fj+1,  nx_f)]
                       + rh_fine[idx(fi+1,  fj+1,  nx_f)] );
            r2h_coarse[idx(I, J, nx_c)] = v;
        }
    }
    // Ghosts of r2h (border) can be left as-is; they’re ignored by interior loops.
}


void prolong_replicate_add_flat(const std::vector<double>& e2h_coarse,
                                int nx_c, int ny_c,
                                std::vector<double>& phi_fine,
                                int nx_f, int ny_f)
{
    for (int J = 1; J <= ny_c; ++J) {
        for (int I = 1; I <= nx_c; ++I) {
            const double v = e2h_coarse[idx(I, J, nx_c)];
            const int i = 2*I - 1;
            const int j = 2*J - 1;
            phi_fine[idx(i,  j,  nx_f)] += v;
            phi_fine[idx(i+1,j,  nx_f)] += v;
            phi_fine[idx(i,  j+1,nx_f)] += v;
            phi_fine[idx(i+1,j+1,nx_f)] += v;
        }
    }
    bcGS_flat(phi_fine, nx_f, ny_f);
}


void smooth_gs_flat(std::vector<double>& phi,
                    const std::vector<double>& f,
                    int nx, int ny, double h,
                    int niter)
{
    const double h2 = h * h;
    for (int it = 0; it < niter; ++it) {
        for (int j = 1; j <= ny; ++j) {
            for (int i = 1; i <= nx; ++i) {
                const int p = idx(i,j,nx);
                phi[p] = 0.25 * (
                    phi[idx(i-1,j,nx)] + phi[idx(i+1,j,nx)] +
                    phi[idx(i,  j-1,nx)] + phi[idx(i,  j+1,nx)]
                ) - 0.25 * h2 * f[p];
            }
        }
        bcGS_flat(phi, nx, ny);  // keep ghosts identical to your 2D path
    }
}


void residual_flat_omp(MGLevel& lv)
{
#ifdef _OPENMP
    const int nx = lv.nx, ny = lv.ny;
    const double invh2 = 1.0 / (lv.h * lv.h);

    // For tiny grids, serial is faster than threading overhead
    if (1LL * nx * ny < 4096) { residual_flat(lv); return; }

    auto& phi = lv.phi;
    auto& f   = lv.f;
    auto& res = lv.res;

    #pragma omp parallel for schedule(static)
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            const int p = idx(i, j, nx);
            const double lap =
                (phi[idx(i-1,j,nx)] + phi[idx(i+1,j,nx)] +
                 phi[idx(i,  j-1,nx)] + phi[idx(i,  j+1,nx)]
                 - 4.0 * phi[p]) * invh2;
            res[p] = f[p] - lap;
        }
    }
#else
    residual_flat(lv);
#endif
}

// ------------------------------------------------------------
// restrict (2x2 mean): fine residual -> coarse RHS
// ------------------------------------------------------------
void restrict_flat_omp(const std::vector<double>& rh_fine, int nx_f, int ny_f,
                       std::vector<double>& r2h_coarse, int nx_c, int ny_c)
{
#ifdef _OPENMP
    if (1LL * nx_c * ny_c < 2048) { // tiny -> serial
        restrict_flat(rh_fine, nx_f, ny_f, r2h_coarse, nx_c, ny_c);
        return;
    }

    #pragma omp parallel for schedule(static)
    for (int J = 1; J <= ny_c; ++J) {
        const int fj = 2*J - 1;
        for (int I = 1; I <= nx_c; ++I) {
            const int fi = 2*I - 1;
            r2h_coarse[idx(I, J, nx_c)] = 0.25 * (
                rh_fine[idx(fi,   fj,   nx_f)] + rh_fine[idx(fi+1, fj,   nx_f)] +
                rh_fine[idx(fi,   fj+1, nx_f)] + rh_fine[idx(fi+1, fj+1, nx_f)]
            );
        }
    }
#else
    restrict_flat(rh_fine, nx_f, ny_f, r2h_coarse, nx_c, ny_c);
#endif
}

// ------------------------------------------------------------
// prolong (replicate-add) coarse error -> fine solution
// bcGS_flat is called *after* the parallel loop (not inside).
// ------------------------------------------------------------
void prolong_replicate_add_flat_omp(const std::vector<double>& e2h_coarse,
                                    int nx_c, int ny_c,
                                    std::vector<double>& phi_fine,
                                    int nx_f, int ny_f)
{
#ifdef _OPENMP
    if (1LL * nx_f * ny_f < 4096) {
        // If your existing prolong_flat already writes into a scratch buffer,
        // you can keep using your serial "add" variant here instead:
        // prolong_replicate_add_flat(e2h_coarse, nx_c, ny_c, phi_fine, nx_f, ny_f);
        // return;
    }

    #pragma omp parallel for schedule(static)
    for (int J = 1; J <= ny_c; ++J) {
        for (int I = 1; I <= nx_c; ++I) {
            const double v = e2h_coarse[idx(I, J, nx_c)];
            const int i = 2*I - 1;
            const int j = 2*J - 1;
            // Each (I,J) writes a unique 2x2 patch -> no overlap
            phi_fine[idx(i,  j,  nx_f)] += v;
            phi_fine[idx(i+1,j,  nx_f)] += v;
            phi_fine[idx(i,  j+1,nx_f)] += v;
            phi_fine[idx(i+1,j+1,nx_f)] += v;
        }
    }
    bcGS_flat(phi_fine, nx_f, ny_f); // outside the parallel region
#else
    prolong_replicate_add_flat(e2h_coarse, nx_c, ny_c, phi_fine, nx_f, ny_f);
#endif
}

double RelativeResidualNorm_from_buffer_omp(const std::vector<double>& res,
                                            const std::vector<double>& f,
                                            int nx,int ny)
{
#ifdef _OPENMP
    double max_r = 0.0;
    double max_f = 0.0;

    #pragma omp parallel for reduction(max:max_r, max_f) schedule(static)
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            const int p = idx(i,j,nx);
            const double ar = std::abs(res[p]);
            const double af = std::abs(f[p]);
            if (ar > max_r) max_r = ar;
            if (af > max_f) max_f = af;
        }
    }
    return (max_f > 0.0) ? (max_r / max_f) : 0.0;
#else
    return RelativeResidualNorm_from_buffer(res, f, nx, ny, /*max_f=*/0.0);
#endif
}

// ========= RBGS smoother =========
void smooth_rbgs_flat_omp(std::vector<double>& phi,
                          const std::vector<double>& f,
                          int nx,int ny,double h,int iters)
{
    const double h2 = h*h;

#ifdef _OPENMP
    for (int k = 0; k < iters; ++k) {
        // color = 0 (even), color = 1 (odd)
        for (int color = 0; color < 2; ++color) {
            #pragma omp parallel for schedule(static)
            for (int j = 1; j <= ny; ++j) {
                // start i so that (i + j) % 2 == color
                int i0 = 1 + ((j + color) & 1);   // 1 or 2
                for (int i = i0; i <= nx; i += 2) {
                    const int p = idx(i,j,nx);
                    phi[p] = 0.25 * (
                        phi[idx(i-1,j,nx)] + phi[idx(i+1,j,nx)] +
                        phi[idx(i,  j-1,nx)] + phi[idx(i,  j+1,nx)]
                    ) - 0.25 * h2 * f[p];
                }
            }
        }
        // keep ghosts consistent (outside parallel region)
        bcGS_flat(phi, nx, ny);
    }
#else
    // if not building with OpenMP, fall back to your serial GS
    smooth_gs_flat(phi, f, nx, ny, h, iters);
#endif
}


// ========= OpenMP V-cycle =========
void vcycle_openmp(MG& mg, int level, int nu1, int nu2, int nCoarse)
{
    auto& fine = mg.L[level];

    // Pre-smoothing
    if (nu1 > 0)
        smooth_rbgs_flat_omp(fine.phi, fine.f, fine.nx, fine.ny, fine.h, nu1);

    const bool hasCoarser = (level + 1 < (int)mg.L.size());
    if (hasCoarser) {
        auto& coarse = mg.L[level + 1];

        // relationships must hold (leave these asserts while testing)
        // assert(fine.nx % 2 == 0 && fine.ny % 2 == 0);
        // assert(coarse.nx == fine.nx/2 && coarse.ny == fine.ny/2);

        // Residual on fine
        residual_flat_omp(fine); // fills fine.res

        // Restrict fine.res -> coarse.f
        restrict_flat_omp(fine.res, fine.nx, fine.ny,
                          coarse.f,  coarse.nx, coarse.ny);

        // Zero coarse solution
        std::fill(coarse.phi.begin(), coarse.phi.end(), 0.0);

        // Recurse
        vcycle_openmp(mg, level + 1, nu1, nu2, nCoarse);

        // Prolong (replicate-add) and BC
        prolong_replicate_add_flat_omp(coarse.phi, coarse.nx, coarse.ny,
                                       fine.phi,  fine.nx,  fine.ny);
    } else {
        // Coarsest relax
        if (nCoarse > 0)
            smooth_rbgs_flat_omp(fine.phi, fine.f, fine.nx, fine.ny, fine.h, nCoarse);
    }

    // Post-smoothing
    if (nu2 > 0)
        smooth_rbgs_flat_omp(fine.phi, fine.f, fine.nx, fine.ny, fine.h, nu2);
}






void vcycle_serial(MG& mg, int level, int nu1, int nu2, int nCoarse)
{
    auto& fine = mg.L[level];



    // Pre-smoothing
    if (nu1 > 0)
        smooth_gs_flat(fine.phi, fine.f, fine.nx, fine.ny, fine.h, nu1);

    const bool hasCoarser = (level + 1 < (int)mg.L.size());
    if (hasCoarser) {
        auto& coarse = mg.L[level + 1];

        // Residual on fine
        residual_flat(fine);  // fills fine.res = f - A*phi

        // Restrict residual to coarse RHS
        restrict_flat(fine.res, fine.nx, fine.ny,
                      coarse.f,  coarse.nx, coarse.ny);

        // Zero coarse solution guess
        std::fill(coarse.phi.begin(), coarse.phi.end(), 0.0);

        // Recurse
        vcycle_serial(mg, level + 1, nu1, nu2, nCoarse);


        // Prolong coarse error to fine scratch (eh_fine) and add to phi
        prolong_replicate_add_flat(coarse.phi, coarse.nx, coarse.ny,
                            fine.phi,  fine.nx,  fine.ny);


    } else {
        // Coarsest "solve": a few extra smooths
        if (nCoarse > 0)
            smooth_gs_flat(fine.phi, fine.f, fine.nx, fine.ny, fine.h, nCoarse);
    }

    // Post-smoothing
    if (nu2 > 0)
        smooth_gs_flat(fine.phi, fine.f, fine.nx, fine.ny, fine.h, nu2);
}




void poisson_solve(Solution& sol,
                   const std::vector<std::vector<double>>& f2d,
                   const SolverConfig& config,
                   int maxVCycles, double eps,
                   int nu1, int nu2, int nCoarse)

{
    
    auto& phi2d = sol.phi;
    // sizes
    const int M  = (int)phi2d.size();        // rows (y)
    const int N  = (int)phi2d[0].size();     // cols (x)
    const int nx = N - 2;                    // interior cols
    const int ny = M - 2;                    // interior rows
    const double h = config.h;

    // Build hierarchy (once per call)
    MG mg; 
    mg_build_auto(mg, nx, ny, h);
    // Level-0 data
    mg.L[0].phi = flatten_2d(phi2d);
    mg.L[0].f   = flatten_2d(f2d);
    mg.L[0].res.assign((nx+2)*(ny+2), 0.0);

    // Loop V-cycles until converged
    for (int it = 0; it < maxVCycles; ++it) {
        // residual on finest level for stopping check
        residual_flat(mg.L[0]); // fills mg.L[0].res = f - A*phi

        double max_f = maxF_inf_flat(mg.L[0].f,nx,ny);

        double rr = RelativeResidualNorm_from_buffer(
                        mg.L[0].res, mg.L[0].f, nx, ny, max_f);
        // (Optional) std::cout << "V-cycle " << it << "  rr=" << rr << "\n";
        if (rr <= eps) break;

        // Backend dispatch
        switch (config.ellipticBackend) {
            case EllipticBackend::Serial:
                vcycle_serial(mg, /*level=*/0, nu1, nu2, nCoarse);
                break;

            case EllipticBackend::OpenMP:
                vcycle_openmp(mg, /*level=*/0, nu1, nu2, nCoarse);
                break;

            case EllipticBackend::CUDA:
                // TODO: CUDA backend (Jacobi/chebyshev smoother); fall back to serial for now.
                vcycle_serial(mg, /*level=*/0, nu1, nu2, nCoarse);
                break;
        }
    }

    // Return result to 2D container expected by the rest of your code
    unflatten_2d(mg.L[0].phi, M, N, sol.phi);

}




std::vector<double> flatten_2d(const std::vector<std::vector<double>>& a) {
    const int M = (int)a.size();        // rows
    const int N = (int)a[0].size();     // cols
    std::vector<double> out;
    out.reserve((size_t)M * N);
    for (int j = 0; j < M; ++j) {       // row-major copy
        out.insert(out.end(), a[j].begin(), a[j].end());
    }
    return out;
}

void unflatten_2d(const std::vector<double>& flat, int M, int N,
             std::vector<std::vector<double>>& a) {
    a.assign(M, std::vector<double>(N, 0.0));
    size_t k = 0;
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < N; ++i)
            a[j][i] = flat[k++];
}


