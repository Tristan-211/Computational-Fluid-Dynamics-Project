#include "parabolic.h"
#include "utils.h"
#include "bc.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <type_traits>

void parabolic_CN1_u(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time){

    double nu = 1.0 / config.Re;
    size_t M = sol.u.size() - 1;
    size_t N = sol.u[0].size() - 2;

    //const auto& Qu = sol.Qu;
    std::vector<std::vector<double>> a(M+1, std::vector<double>(N+2,0.0));
    std::vector<std::vector<double>> b(M+1, std::vector<double>(N+2,0.0));
    std::vector<std::vector<double>> c(M+1, std::vector<double>(N+2,0.0));
    std::vector<std::vector<double>> d(M+1, std::vector<double>(N+2,0.0));


    double delx = config.h;
    double dely = delx;

    double dx = nu * config.dt / (delx * delx);
    double dy = nu * config.dt / (dely * dely);
    double d1 = dx / 2.0;
    double d2 = dy / 2.0;

    // Fill a, b, c with constants
    for (size_t i = 0; i <= M; ++i)
        for (size_t j = 0; j <= N+1; ++j) {
            a[i][j] = -d1;
            b[i][j] = 1.0 + 2.0 * d1;
            c[i][j] = -d1;
        }
    

    // Fill d with internal values

    for (size_t i = 1; i < M; ++i)
        for (size_t j = 1; j <= N; ++j)
            d[i][j] = d2 * sol.u[i][j+1] + (1.0 - 2.0*d2) * sol.u[i][j] + d2 * sol.u[i][j-1] + sol.Qu[i][j] * config.dt * 0.5;
        

    //saveMatrixToFile(d, "d.csv");
    bcCN1_u(a, b, c, d, mesh, config, time);
    

    
    for (size_t j = 1; j <= N; ++j) {  // MATLAB j=2:N+1 corresponds to C++ j=1 to N
        std::vector<double> a_col(M - 1);
        std::vector<double> b_col(M - 1);
        std::vector<double> c_col(M - 1);
        std::vector<double> d_col(M - 1);

        for (size_t i = 1; i < M; ++i) {
            a_col[i - 1] = a[i][j];
            b_col[i - 1] = b[i][j];
            c_col[i - 1] = c[i][j];
            d_col[i - 1] = d[i][j];
        }

        std::vector<double> result = solveTriDiag(a_col, b_col, c_col, d_col);

        for (size_t i = 1; i < M; ++i) {
            sol.u[i][j] = result[i - 1];
        }




    }

    bc_u(sol.u,mesh,config, time + config.dt / 2.0);
    
}

void parabolic_CN2_u(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time){

    double nu = 1.0 / config.Re;
    size_t M = sol.u.size() - 1;
    size_t N = sol.u[0].size() - 2;

    const auto& Qu = sol.Qu;
    std::vector<std::vector<double>> a(M+1, std::vector<double>(N+2,0.0));
    std::vector<std::vector<double>> b(M+1, std::vector<double>(N+2,0.0));
    std::vector<std::vector<double>> c(M+1, std::vector<double>(N+2,0.0));
    std::vector<std::vector<double>> d(M+1, std::vector<double>(N+2,0.0));


    double delx = config.h;
    double dely = delx;


    double dy = nu * config.dt / (dely * dely);

    double d2 = dy / 2.0;

    // Fill a, b, c with constants
    for (size_t i = 0; i <= M; ++i)
        for (size_t j = 0; j <= N+1; ++j) {
            a[i][j] = -d2;
            b[i][j] = 1.0 + 2.0 * d2;
            c[i][j] = -d2;
        }

    // Fill d with internal values
    for (size_t i = 1; i < M; ++i) {             // MATLAB i = 2:M
        for (size_t j = 1; j <= N; ++j)         // MATLAB j = 2:N+1
            d[i][j] = d2 * sol.u[i+1][j] + (1.0 - 2.0 * d2) * sol.u[i][j] + d2 * sol.u[i-1][j] + Qu[i][j] * config.dt * 0.5;
        
}

    bcCN2_u(a, b, c, d, mesh, config, time);

    for (size_t i = 1; i <= M; ++i) {  // MATLAB i=2:M+1
        std::vector<double> a_row(N);
        std::vector<double> b_row(N);
        std::vector<double> c_row(N);
        std::vector<double> d_row(N);

        for (size_t j = 1; j <= N; ++j) {  // extract row from columns 2:N+1 → j=1:N
            a_row[j - 1] = a[i][j];
            b_row[j - 1] = b[i][j];
            c_row[j - 1] = c[i][j];
            d_row[j - 1] = d[i][j];
        }

        std::vector<double> result = solveTriDiag(a_row, b_row, c_row, d_row);

        for (size_t j = 1; j <= N; ++j) {
            sol.u[i][j] = result[j - 1];
        }

        }

    bc_u(sol.u, mesh, config, time + config.dt);


}



void parabolic_CN1_v(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time){

    double nu = 1.0 / config.Re;
    size_t M = sol.v.size() - 2;      
    size_t N = sol.v[0].size() - 1;   

    // Allocate matrices
    const auto& Qv = sol.Qv;
    std::vector<std::vector<double>> a(M + 2, std::vector<double>(N + 1, 0.0));
    std::vector<std::vector<double>> b(M + 2, std::vector<double>(N + 1, 0.0));
    std::vector<std::vector<double>> c(M + 2, std::vector<double>(N + 1, 0.0));
    std::vector<std::vector<double>> d(M + 2, std::vector<double>(N + 1, 0.0));

    // Parameters
    double delx = config.h;
    double dely = delx;

    double dx = nu * config.dt / (delx * delx);
    double dy = nu * config.dt / (dely * dely);
    double d1 = dx / 2.0;
    double d2 = dy / 2.0;

    // Fill a, b, c with constants
    for (size_t i = 0; i < M + 2; ++i) {
        for (size_t j = 0; j < N + 1; ++j) {
            a[i][j] = -d1;
            b[i][j] = 1.0 + 2.0 * d1;
            c[i][j] = -d1;
        }
    }

    // Fill d with internal values (2:M+1, 2:N in MATLAB => 1 to Mv, 1 to Nv-1)

    //saveMatrixToFile(Qv, "Qv.csv");
    for (size_t i = 1; i <= M; ++i) {
        for (size_t j = 1; j < N; ++j) 
            d[i][j] = d2 * sol.v[i][j+1] + (1.0 - 2.0 * d2) * sol.v[i][j] + d2 * sol.v[i][j-1] + Qv[i][j] * config.dt * 0.5;
        
    }

    
    bcCN1_v(a, b, c, d, mesh, config, time);


    for (size_t j = 1; j < N; ++j) {  // MATLAB j = 2:N
        std::vector<double> a_col(M);
        std::vector<double> b_col(M);
        std::vector<double> c_col(M);
        std::vector<double> d_col(M);

        for (size_t i = 1; i <= M; ++i) {
            a_col[i - 1] = a[i][j];
            b_col[i - 1] = b[i][j];
            c_col[i - 1] = c[i][j];
            d_col[i - 1] = d[i][j];
        }

        std::vector<double> result = solveTriDiag(a_col, b_col, c_col, d_col);


        for (size_t i = 1; i <= M; ++i) {
            sol.v[i][j] = result[i - 1];
    }
    }
    bc_v(sol.v, mesh, config, time + config.dt / 2.0);

}

void parabolic_CN2_v(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time){

    double nu = 1.0 / config.Re;
    size_t M = sol.v.size() - 2;      
    size_t N = sol.v[0].size() - 1;   

    // Allocate matrices
    const auto& Qv = sol.Qv;
    std::vector<std::vector<double>> a(M + 2, std::vector<double>(N + 1, 0.0));
    std::vector<std::vector<double>> b(M + 2, std::vector<double>(N + 1, 0.0));
    std::vector<std::vector<double>> c(M + 2, std::vector<double>(N + 1, 0.0));
    std::vector<std::vector<double>> d(M + 2, std::vector<double>(N + 1, 0.0));

    // Parameters
    double delx = config.h;
    double dely = delx;

    double dx = nu * config.dt / (delx * delx);
    double dy = nu * config.dt / (dely * dely);
    double d1 = dx / 2.0;
    double d2 = dy / 2.0;

    // Fill a, b, c with constants
    for (size_t i = 0; i < M + 2; ++i) {
        for (size_t j = 0; j < N + 1; ++j) {
            a[i][j] = -d1;
            b[i][j] = 1.0 + 2.0 * d1;
            c[i][j] = -d1;
        }
    }

    // Fill d with internal values (2:M+1, 2:N in MATLAB => 1 to Mv, 1 to Nv-1)
    for (size_t i = 1; i <= M; ++i) {
        for (size_t j = 1; j < N; ++j) 
            d[i][j] = d2 * sol.v[i+1][j] + (1.0 - 2.0 * d2) * sol.v[i][j] + d2 * sol.v[i-1][j] + Qv[i][j] * config.dt * 0.5;
        
    }

    bcCN2_v(a, b, c, d, mesh, config, time);

    for (size_t i = 1; i <= M; ++i) {  // MATLAB i = 2:M+1
        std::vector<double> a_row(N - 1);
        std::vector<double> b_row(N - 1);
        std::vector<double> c_row(N - 1);
        std::vector<double> d_row(N - 1);

        for (size_t j = 1; j < N; ++j) {  // MATLAB j = 2:N
            a_row[j - 1] = a[i][j];
            b_row[j - 1] = b[i][j];
            c_row[j - 1] = c[i][j];
            d_row[j - 1] = d[i][j];
        }

        std::vector<double> result = solveTriDiag(a_row, b_row, c_row, d_row);


        for (size_t j = 1; j < N; ++j) {
            sol.v[i][j] = result[j - 1];
        }
    }
    bc_v(sol.v, mesh, config, time + config.dt);
}


void parabolic_CN1_T(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time){
    size_t M = sol.T.size();       // total size with ghost cells
    size_t N = sol.T[0].size();

    const auto& QT = sol.QT;
    std::vector<std::vector<double>> a(M, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> b(M, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> c(M, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> d(M, std::vector<double>(N,0.0));

    size_t Mint = M - 6;
    size_t Nint = N - 6;

    double nu = 1.0 / (config.Re * config.Pr);
    double h = config.Lx / config.M;  // uniform spacing in x and y

    double dx = nu * config.dt / (h * h);
    double dy = dx;                   // same spacing in both directions
    double d1 = dx / 2.0;
    double d2 = dy / 2.0;

    // Fill a, b, c
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j) {
            a[i][j] = -d1;
            b[i][j] = 1.0 + 2.0 * d1;
            c[i][j] = -d1;
        }

    // Fill d using internal cell indices (skip 3 ghost layers on each side)
    for (size_t i = 3; i < Mint + 3; ++i)
        for (size_t j = 3; j < Nint + 3; ++j)
            d[i][j] = d2 * sol.T[i][j+1] + (1.0 - 2.0 * d2) * sol.T[i][j] + d2 * sol.T[i][j-1] + QT[i][j] * config.dt * 0.5;


    bcCN1_T(a, b, c, d, mesh, config, time);

    for (size_t j = 3; j < Nint + 3; ++j) {
        std::vector<double> a_col(Mint);
        std::vector<double> b_col(Mint);
        std::vector<double> c_col(Mint);
        std::vector<double> d_col(Mint);

        for (size_t i = 3; i < Mint + 3; ++i) {
            size_t idx = i - 3;
            a_col[idx] = a[i][j];
            b_col[idx] = b[i][j];
            c_col[idx] = c[i][j];
            d_col[idx] = d[i][j];
        }

        std::vector<double> result = solveTriDiag(a_col, b_col, c_col, d_col);
        

        for (size_t i = 3; i < Mint + 3; ++i) {
            sol.T[i][j] = result[i - 3];
        }
    } 

    bc_T(sol.T, mesh, config, time + config.dt / 2.0);

}



void parabolic_CN2_T(Solution& sol, const Mesh& mesh, const SolverConfig& config, double time){
    size_t M = sol.T.size();       // total size with ghost cells
    size_t N = sol.T[0].size();

    const auto& QT = sol.QT;
    std::vector<std::vector<double>> a(M, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> b(M, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> c(M, std::vector<double>(N,0.0));
    std::vector<std::vector<double>> d(M, std::vector<double>(N,0.0));

    size_t Mint = M - 6;
    size_t Nint = N - 6;

    double nu = 1.0 / (config.Re * config.Pr);
    double h = config.Lx / config.M;  // uniform spacing in x and y

    double dx = nu * config.dt / (h * h);
    double dy = dx;                   // same spacing in both directions
    double d1 = dx / 2.0;
    double d2 = dy / 2.0;

    // Fill a, b, c
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j) {
            a[i][j] = -d1;
            b[i][j] = 1.0 + 2.0 * d1;
            c[i][j] = -d1;
        }

    // Fill d using internal cell indices (skip 3 ghost layers on each side)
    for (size_t i = 3; i < Mint + 3; ++i)
        for (size_t j = 3; j < Nint + 3; ++j)
            d[i][j] = d2 * sol.T[i+1][j] + (1.0 - 2.0 * d2) * sol.T[i][j] + d2 * sol.T[i-1][j] + QT[i][j] * config.dt * 0.5;


    bcCN2_T(a, b, c, d, mesh, config, time);

    for (size_t i = 3; i < Mint + 3; ++i) {
        std::vector<double> a_row(Nint);
        std::vector<double> b_row(Nint);
        std::vector<double> c_row(Nint);
        std::vector<double> d_row(Nint);

        for (size_t j = 3; j < Nint + 3; ++j) {
            size_t idx = j - 3;
            a_row[idx] = a[i][j];
            b_row[idx] = b[i][j];
            c_row[idx] = c[i][j];
            d_row[idx] = d[i][j];
        }

        std::vector<double> result = solveTriDiag(a_row, b_row, c_row, d_row);

        for (size_t j = 3; j < Nint + 3; ++j) {
            sol.T[i][j] = result[j - 3];
        }
    }


    bc_T(sol.T, mesh, config, time + config.dt);

}