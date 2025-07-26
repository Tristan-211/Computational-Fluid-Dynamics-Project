#pragma once
#include <vector>
#include "config.h"
#include "solution.h"
#include "mesh.h"

void gaussSeidel(std::vector<std::vector<double>>& phi,const std::vector<std::vector<double>>& f,double h,int niter); // simple Gauss Seidel solver

void correctOutletMassConservation(Solution& sol,const Mesh& mesh,const SolverConfig& config); // corrects outlet velocities to conserve mass based on inlet velocityies (generalized)


//mesh restriction and prolongation functions 
std::vector<std::vector<double>> restrict(const std::vector<std::vector<double>>& rh);
std::vector<std::vector<double>> prolong(const std::vector<std::vector<double>>& e2h);

//calculates cell center velocity divergence 
std::vector<std::vector<double>> calcDivV(const Solution& sol, const Mesh& mesh) ;

// v cycle muiltigrid solver
std::vector<std::vector<double>> multigrid(std::vector<std::vector<double>> phi,const std::vector<std::vector<double>>& f,double h);


//poisson solver
std::vector<std::vector<double>> poissonSolver(std::vector<std::vector<double>> phi,const std::vector<std::vector<double>>& f,const Mesh& mesh,int nIterMax,double epsilon);

//porjects velocity vector fields using lagrange multiplier
void projectV(Solution& sol,const Mesh& mesh,const SolverConfig& config,double t);