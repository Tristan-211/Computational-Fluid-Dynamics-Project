#pragma once
#include <vector>
#include <string>
#include "config.h"
#include "mesh.h"
#include "solution.h"


std::vector<double> linspace(double start, double end, int num); //linspace function works the same as MATLAB 

std::vector<double> solveTriDiag( //Tri diagonal system solver
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d);

void saveMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename); //saves a matrix to a .csv file
bool calcDt(double t, double outputTime, SolverConfig& config, const Solution& sol); //calculates dt for parabolic solver 

void meshgrid(const std::vector<double>& x,const std::vector<double>& y,std::vector<std::vector<double>>& X,std::vector<std::vector<double>>& Y); //same as meshgrid function from MATLAB


//equivlent to MATLAB element wise matrix operations 
std::vector<std::vector<double>> elementwiseAdd(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> elementwiseSubtract(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> elementwiseMultiply(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> elementwiseDivide(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> scalarMultiply(const std::vector<std::vector<double>>& A,double scalar);

// residual functions 
double RelativeResidualNorm(const std::vector<std::vector<double>>& phi,const std::vector<std::vector<double>>& f,double h);
std::vector<std::vector<double>> Residual(const std::vector<std::vector<double>>& phi,const std::vector<std::vector<double>>& f,double h);

double maxAbs(const std::vector<std::vector<double>>& mat);

std::vector<std::vector<double>> zeros(size_t rows, size_t cols);

void printTimeProgress(double t, double totalTime);

//matrix multiplier
std::vector<std::vector<double>> matMul(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B);