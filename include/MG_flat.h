#pragma once
#include "MG_flat_types.h"
#include "solution.h"
#include <vector>
#include "config.h"


void residual_flat(MGLevel& lv);


// Build hierarchy with as many coarsenings by 2 as possible
void mg_build_auto(MG& mg, int nx0, int ny0, double h0);

double maxF_inf_flat(const std::vector<double>& f, int nx, int ny);



void restrict_flat(const std::vector<double>& rh_fine,int nx_f, int ny_f,std::vector<double>& r2h_coarse,int nx_c, int ny_c);

void prolong_replicate_add_flat(const std::vector<double>& e2h_coarse,
                                int nx_c, int ny_c,
                                std::vector<double>& phi_fine,
                                int nx_f, int ny_f);

void smooth_gs_flat(std::vector<double>& phi,const std::vector<double>& f,int nx, int ny, double h,int niter);

void vcycle_serial(MG& mg, int level, int nu1, int nu2, int nCoarse);

void poisson_solve(Solution& sol,
                   const std::vector<std::vector<double>>& f2d,
                   const SolverConfig& config,
                   int maxVCycles, double eps,
                   int nu1 = 1, int nu2 = 1, int nCoarse = 3);




//open MP versions
void residual_flat_omp(MGLevel& lv);

void restrict_flat_omp(const std::vector<double>& rh_fine, int nx_f, int ny_f,
                       std::vector<double>& r2h_coarse, int nx_c, int ny_c);

void prolong_replicate_add_flat_omp(const std::vector<double>& e2h_coarse,
                                    int nx_c, int ny_c,
                                    std::vector<double>& phi_fine,
                                    int nx_f, int ny_f);
// --- OpenMP smoother (RBGS) ---
void smooth_rbgs_flat_omp(std::vector<double>& phi,
                          const std::vector<double>& f,
                          int nx,int ny,double h,int iters);

// --- OpenMP V-cycle ---
void vcycle_openmp(MG& mg, int level, int nu1, int nu2, int nCoarse);

// --- Optional: parallel RR norm (∞-norm over interior) ---
double RelativeResidualNorm_from_buffer_omp(const std::vector<double>& res,
                                            const std::vector<double>& f,
                                            int nx,int ny);                                    




// Flatten 2D vector-with-ghosts (size (nx+2)x(ny+2)) to row-major 1D
std::vector<double> flatten_2d(const std::vector<std::vector<double>>& a);

void unflatten_2d(const std::vector<double>& flat, int nx,int ny,std::vector<std::vector<double>>& a);


double RelativeResidualNorm_from_buffer(const std::vector<double>& res,const std::vector<double>& f,int nx, int ny, double max_f /*precomputed or 0*/);
