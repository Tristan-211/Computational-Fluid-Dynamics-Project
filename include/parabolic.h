#pragma once
#include <vector>
#include "mesh.h"
#include "config.h"
#include "solution.h"

// Each parabolic Crank-Nicolson ADI step takes:
// - Solution& sol: the full solution struct (u, v, T, Qu, Qv, QT)
// - const Mesh& mesh: mesh information (xf, xc, etc.)
// - const SolverConfig& config: simulation parameters (Re, Pr, dt, etc.)
// - double time: current simulation time

void parabolic_CN1_u(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time);

void parabolic_CN2_u(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time);

void parabolic_CN1_v(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time);

void parabolic_CN2_v(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time);

void parabolic_CN1_T(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time);

void parabolic_CN2_T(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time);

