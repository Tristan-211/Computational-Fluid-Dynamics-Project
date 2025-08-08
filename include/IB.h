#pragma once
#include <vector>
#include <string>
#include "config.h"
#include "mesh.h"
#include "solution.h"



//Mask Functions 
std::vector<std::vector<bool>> makeMask(const std::vector<std::vector<double>>& XR,const std::vector<std::vector<double>>& YR,double IBlh,double IBwh);
void applyMaskSet(std::vector<std::vector<double>>& target,const std::vector<std::vector<bool>>& mask,double value);

//checks that rotors fit in bounds 
void checkRotorBounds(const std::vector<Rotor>& rotors, const Mesh& mesh);

//reshapes a 1D vector into a 2d array of size rows x cols 
std::vector<std::vector<double>> reshape(const std::vector<double>& flat, size_t rows, size_t cols);

void initIBMeshU(const Mesh& mesh,SolverConfig& config,const  Solution& sol);

void calcSourceIB(Solution& sol,const Mesh& mesh,const SolverConfig& config,double t);