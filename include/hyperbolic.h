#pragma once
#include <utility> // for std::pair
#include <vector>
#include <cmath>
#include "solution.h"
#include "config.h"
#include "mesh.h"
using Array2D = std::vector<std::vector<double>>;

std::pair<Array2D, Array2D> hyperbolic_uv_2D(Array2D& u,Array2D& v,const SolverConfig& config); //u-v hyperbolic forcing term calculator 

double psiWENO(double a, double b, double c, double d) noexcept; //psiWENO value calculator 

Array2D calc_adTdx_WENO_2D(const Array2D& T, const Array2D& a, const SolverConfig& config);
Array2D calc_bdTdy_WENO_2D(const Array2D& T, const Array2D& b, const SolverConfig& config);

Array2D hyperbolic_T_WENO_2D(const Solution& sol, const Array2D& u_old, const Array2D& v_old, double t,const SolverConfig& config,const Mesh& mesh);