#pragma once
#include <vector>
#include "mesh.h"
#include "config.h"

void bc_u(std::vector<std::vector<double>>& u, const Mesh& mesh, const SolverConfig& config, double time);
void bc_v(std::vector<std::vector<double>>& v, const Mesh& mesh, const SolverConfig& config, double time);
void bc_T(std::vector<std::vector<double>>& T, const Mesh& mesh, const SolverConfig& config, double time);

void bcCN1_u(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time);

void bcCN2_u(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time);

void bcCN1_v(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time);

void bcCN2_v(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time);

void bcCN1_T(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time);

void bcCN2_T(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time);


void bcGS(std::vector<std::vector<double>>& phi);


void bcGhost_u(std::vector<std::vector<double>>& u, const Mesh& mesh, const SolverConfig& config, double time);
void bcGhost_v(std::vector<std::vector<double>>& v, const Mesh& mesh, const SolverConfig& config, double time);