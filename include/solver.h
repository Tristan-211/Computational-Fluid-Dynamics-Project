#pragma once

#include "solution.h"
#include "mesh.h"
#include "config.h"


void runSolver(Solution& sol, const Mesh& mesh, SolverConfig& config);

Solution initializeSolution(const Mesh& mesh, const SolverConfig& config);

