#include "bc.h"
#include <iostream>
#include "MG_flat_types.h"

void bc_u(std::vector<std::vector<double>>& u, const Mesh& mesh, const SolverConfig& config, double time) {
    size_t Nx = mesh.xf.size();
    size_t Ny = mesh.yc.size();

    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second)
                        u[i][0] = 2 * bc.uFunction(mesh.xf[i], 0.0, time) - u[i][1];
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second)
                        u[i][Ny - 1] = 2 * bc.uFunction(mesh.xf[i], config.Ly, time) - u[i][Ny - 2];
            }
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second)
                        u[0][j] = bc.uFunction(0.0, mesh.yc[j], time);
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second)
                        u[Nx - 1][j] = bc.uFunction(config.Lx, mesh.yc[j], time);
            }
        } else if (bc.type == "outlet") {
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second)
                        u[0][j] = 4.0 / 3.0 * u[1][j] - 1.0 / 3.0 * u[2][j]; // 2nd-order forward
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second)
                        u[Nx - 1][j] = 4.0 / 3.0 * u[Nx - 2][j] - 1.0 / 3.0 * u[Nx - 3][j]; // 2nd-order backward
            }
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second)
                        u[i][0] = u[i][1]; // central
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second)
                        u[i][Ny - 1] = u[i][Ny - 2]; // central
            }
        }
    }

    // Wall fallback
    for (size_t i = 0; i < Nx; ++i) {
        bool bottom = false, top = false;
        for (const auto& bc : config.boundaries) {
            if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second) {
                if (bc.side == "bottom") bottom = true;
                if (bc.side == "top") top = true;
            }
        }
        if (!bottom) u[i][0] = -u[i][1];
        if (!top) u[i][Ny - 1] = -u[i][Ny - 2];
    }

    for (size_t j = 0; j < Ny; ++j) {
        bool left = false, right = false;
        for (const auto& bc : config.boundaries) {
            if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second) {
                if (bc.side == "left") left = true;
                if (bc.side == "right") right = true;
            }
        }
        if (!left) u[0][j] = 0.0;
        if (!right) u[Nx - 1][j] = 0.0;
    }
}

void bc_v(std::vector<std::vector<double>>& v, const Mesh& mesh, const SolverConfig& config, double time) {
    size_t Nx = mesh.xc.size();
    size_t Ny = mesh.yf.size();

    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second)
                        v[i][0] = bc.vFunction(mesh.xc[i], 0.0, time);
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second)
                        v[i][Ny - 1] = bc.vFunction(mesh.xc[i], mesh.yf[Ny - 1], time);
            }
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second)
                        v[0][j] = 2 * bc.vFunction(0.0, mesh.yf[j], time) - v[1][j];
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second)
                        v[Nx - 1][j] = 2* bc.vFunction(mesh.xc[Nx - 1], mesh.yf[j], time) - v[Nx-2][j];
            }
        } else if (bc.type == "outlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second)
                        v[i][0] = 4.0 / 3.0 * v[i][1] - 1.0 / 3.0 * v[i][2]; // forward
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second)
                        v[i][Ny - 1] = 4.0 / 3.0 * v[i][Ny - 2] - 1.0 / 3.0 * v[i][Ny - 3]; // backward
            }
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second)
                        v[0][j] = v[1][j]; // central
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second)
                        v[Nx - 1][j] = v[Nx - 2][j]; // central
            }
        }
    }

    for (size_t i = 0; i < Nx; ++i) {
        bool bottom = false, top = false;
        for (const auto& bc : config.boundaries) {
            if (bc.side == "bottom" && mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second) bottom = true;
            if (bc.side == "top" && mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second) top = true;
        }
        if (!bottom) v[i][0] = 0.0;
        if (!top) v[i][Ny - 1] = 0.0;
    }

    for (size_t j = 0; j < Ny; ++j) {
        bool left = false, right = false;
        for (const auto& bc : config.boundaries) {
            if (bc.side == "left" && mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second) left = true;
            if (bc.side == "right" && mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second) right = true;
        }
        if (!left) v[0][j] = -v[1][j];
        if (!right) v[Nx - 1][j] = -v[Nx - 2][j];
    }
}

void bc_T(std::vector<std::vector<double>>& T, const Mesh& mesh, const SolverConfig& config, double time) {
    size_t Nx = mesh.xc3.size();
    size_t Ny = mesh.yc3.size();

    // Apply only Dirichlet inlet BCs from boundary struct
    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) {
                        double val = bc.TFunction(mesh.xc3[i], 0.0, time);
                        T[i][0] = 2 * val - T[i][5];
                        T[i][1] = 2 * val - T[i][4];
                        T[i][2] = 2 * val - T[i][3];
                    }
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) {
                        double val = bc.TFunction(mesh.xc3[i], config.Ly, time);
                        T[i][Ny - 1] = 2 * val - T[i][Ny - 6];
                        T[i][Ny - 2] = 2 * val - T[i][Ny - 5];
                        T[i][Ny - 3] = 2 * val - T[i][Ny - 4];
                    }
            }
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) {
                        double val = bc.TFunction(0.0, mesh.yc3[j], time);
                        T[0][j] = 2 * val - T[5][j];
                        T[1][j] = 2 * val - T[4][j];
                        T[2][j] = 2 * val - T[3][j];
                    }
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) {
                        double val = bc.TFunction(config.Lx, mesh.yc3[j], time);
                        T[Nx - 1][j] = 2 * val - T[Nx - 6][j];
                        T[Nx - 2][j] = 2 * val - T[Nx - 5][j];
                        T[Nx - 3][j] = 2 * val - T[Nx - 4][j];
                    }
            }
        }
    }

    // Walls and outlets use mirrored fallback condition
    for (size_t i = 0; i < Nx; ++i) {
        bool bottom = false, top = false;
        for (const auto& bc : config.boundaries) {
            if (bc.type == "inlet") {
                if (bc.side == "bottom" && mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) bottom = true;
                if (bc.side == "top" && mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) top = true;
            }
        }
        if (!bottom) {
            T[i][0] = T[i][5];
            T[i][1] = T[i][4];
            T[i][2] = T[i][3];
        }
        if (!top) {
            T[i][Ny - 1] = T[i][Ny - 6];
            T[i][Ny - 2] = T[i][Ny - 5];
            T[i][Ny - 3] = T[i][Ny - 4];
        }
    }
    for (size_t j = 0; j < Ny; ++j) {
        bool left = false, right = false;
        for (const auto& bc : config.boundaries) {
            if (bc.type == "inlet") {
                if (bc.side == "left" && mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) left = true;
                if (bc.side == "right" && mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) right = true;
            }
        }
        if (!left) {
            T[0][j] = T[5][j];
            T[1][j] = T[4][j];
            T[2][j] = T[3][j];
        }
        if (!right) {
            T[Nx - 1][j] = T[Nx - 6][j];
            T[Nx - 2][j] = T[Nx - 5][j];
            T[Nx - 3][j] = T[Nx - 4][j];
        }
    }
}



void bcCN1_u(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time) {

    size_t Nx = mesh.xf.size();
    size_t Ny = mesh.yc.size();


    
    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "left")
                for (size_t j = 0; j < Ny; ++j){


                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second){
                        
                        d[1][j] = d[1][j] - a[1][j] * bc.uFunction(0.0, mesh.yc[j], time);
                        a[1][j] = 0.0;
                    }
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second){
                        
                        d[Nx-2][j] = d[Nx-2][j] - c[Nx-2][j] * bc.uFunction(config.Lx, mesh.yc[j], time);
                        c[Nx-2][j] = 0.0;
                    }    
            }
        } else if (bc.type == "outlet") {
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second){
                        b[1][j] = b[1][j] + 4.0/3.0 * a[1][j];
                        c[1][j] = c[1][j] - 1.0/3.0 * a[1][j];
                        a[1][j] = 0.0;
                    }
                        
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second){
                        b[Nx-2][j] = b[Nx-2][j] + 4.0/3.0 * c[Nx-2][j];
                        a[Nx-2][j] = a[Nx-2][j] - 1.0/3.0 * c[Nx-2][j];
                        c[Nx-2][j] = 0.0;

                    }
                        
            }
        }
    }

    // Wall fallback
    for (size_t j = 0; j < Ny; ++j) {
        bool left = false, right = false;
        for (const auto& bc : config.boundaries) {
            if (mesh.yc[j] >= bc.range.first && mesh.yc[j] <= bc.range.second) {
                if (bc.side == "left") left = true;
                if (bc.side == "right") right = true;
            }
        }
        if (!left) a[1][j] = 0.0;
        if (!right) c[Nx-2][j] = 0.0;
    }   

    }


void bcCN2_u(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time) { 

    size_t Nx = mesh.xf.size();
    size_t Ny = mesh.yc.size();

    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second){
                        b[i][1] = b[i][1] - a[i][1];
                        d[i][1] = d[i][1] - a[i][1] * 2 * bc.uFunction(mesh.xf[i], 0.0, time);
                        a[i][1] = 0.0;
                    }
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second){
                        b[i][Ny-2] = b[i][Ny-2] - c[i][Ny-2];
                        d[i][Ny-2] = d[i][Ny-2] - c[i][Ny-2] * 2 * bc.uFunction(mesh.xf[i],config.Ly, time);
                        c[i][Ny-2] = 0.0;
                    }
            }
        } else if (bc.type == "outlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second){
                        b[i][1] = b[i][1] + a[i][1];
                        a[i][1] = 0.0;
                    }
                        
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second){
                        b[i][Ny-2] = b[i][Ny-2] + c[i][Ny-2];
                        c[i][Ny-2] = 0.0;
                    }
                    
            }
        }
    }

    // Wall fallback
    for (size_t i = 0; i < Nx; ++i) {
        bool bottom = false, top = false;
        for (const auto& bc : config.boundaries) {
            if (mesh.xf[i] >= bc.range.first && mesh.xf[i] <= bc.range.second) {
                if (bc.side == "bottom") bottom = true;
                if (bc.side == "top") top = true;
            }
        }
        if (!bottom) {
            b[i][1] = b[i][1] - a[i][1];
            a[i][1] = 0.0;
         }
        if (!top){
            b[i][Ny-2] = b[i][Ny-2] - c[i][Ny-2];
            c[i][Ny-2] = 0.0; 
        }
    }


    }


void bcCN1_v(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time){
    size_t Nx = mesh.xc.size();
    size_t Ny = mesh.yf.size();

    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second){
                        b[1][j] = b[1][j] - a[1][j];
                        d[1][j] = d[1][j] - a[1][j] * 2 * bc.vFunction(0.0, mesh.yf[j], time);
                        a[1][j] = 0.0;
                    }
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second){
                        b[j][Nx - 2] = b[j][Nx - 2] - c[j][Nx - 2];
                        d[j][Nx - 2] = d[j][Nx - 2] - c[j][Nx - 2] * 2 * bc.vFunction(config.Lx, mesh.yf[j], time);
                        c[j][Nx - 2] = 0.0;
                    }

            }
        } else if (bc.type == "outlet") {
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second){
                        b[1][j] = b[1][j] + a[1][j];
                        a[1][j] = 0.0;
                    }

            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second){
                        b[Nx-2][j] = b[Nx-2][j] + c[Nx-2][j];
                        c[Nx-2][j] = 0.0; 
                    }
            }
        }
    }

    for (size_t j = 0; j < Ny; ++j) {
        bool left = false, right = false;
        for (const auto& bc : config.boundaries) {
            if (bc.side == "left" && mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second) left = true;
            if (bc.side == "right" && mesh.yf[j] >= bc.range.first && mesh.yf[j] <= bc.range.second) right = true;
        }
        if (!left) {
            b[1][j] = b[1][j] - a[1][j];
            a[1][j] = 0.0;
        }
        if (!right){
            b[Nx-2][j] = b[Nx-2][j] - c[Nx-2][j];
            c[Nx-2][j] = 0.0;  
        }
    }
        
    }


void bcCN2_v(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time){

    size_t Nx = mesh.xc.size();
    size_t Ny = mesh.yf.size();

    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second){
                        d[i][1] = d[i][1] - a[i][1] * bc.vFunction(mesh.xc[i], 0.0, time);
                        a[i][1] = 0.0; 
                    }

            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second){
                        d[i][Ny-2] = d[i][Ny-2] - c[i][Ny-2] * bc.vFunction(mesh.xc[i], config.Ly, time);
                        c[i][Ny-2] = 0.0;
                    } 
            }
        } else if (bc.type == "outlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second){
                        b[i][1] = b[i][1] + (4.0 / 3.0) * a[i][1];
                        c[i][1] = c[i][1] - (1.0 / 3.0) * a[i][1];
                        a[i][1] = 0.0;
                    }

            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second){
                        b[i][Ny - 2] = b[i][Ny - 2] + (4.0 / 3.0) * c[i][Ny - 2];
                        a[i][Ny - 2] = a[i][Ny - 2] - (1.0 / 3.0) * c[i][Ny - 2];
                        c[i][Ny - 2] = 0.0;
                    }

            }
        }
    }

    for (size_t i = 0; i < Nx; ++i) {
        bool bottom = false, top = false;
        for (const auto& bc : config.boundaries) {
            if (bc.side == "bottom" && mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second) bottom = true;
            if (bc.side == "top" && mesh.xc[i] >= bc.range.first && mesh.xc[i] <= bc.range.second) top = true;
        }
        if (!bottom) a[i][1] = 0.0;
        if (!top) c[i][Ny-2] = 0.0;
    }



    }


void bcCN1_T(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time){
    size_t Nx = mesh.xc3.size();
    size_t Ny = mesh.yc3.size();

    // Apply only Dirichlet inlet BCs from boundary struct
    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "left") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) {
                        double val = bc.TFunction(0.0, mesh.yc3[j], time);
                        for(size_t k = 1; k<=3; ++k){
                            b[k][j] = b[k][j] - a[k][j];
                            d[k][j] = d[j][j] - 2*val*a[k][j];
                            a[k][j] = 0.0;
                        }

                    }
            }
            if (bc.side == "right") {
                for (size_t j = 0; j < Ny; ++j)
                    if (mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) {
                        double val = bc.TFunction(config.Lx, mesh.yc3[j], time);
                        for(size_t k = Nx-4; k<=Nx-2; ++k){
                            b[k][j] = b[k][j] - c[k][j];
                            d[k][j] = d[j][j] - 2*val*c[k][j];
                            c[k][j] = 0.0;
                        }
                    }
            }
        }
    }

    // Walls and outlets use mirrored fallback condition
    for (size_t j = 0; j < Ny; ++j) {
        bool left = false, right = false;
        for (const auto& bc : config.boundaries) {
            if (bc.type == "inlet") {
                if (bc.side == "left" && mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) left = true;
                if (bc.side == "right" && mesh.yc3[j] >= bc.range.first && mesh.yc3[j] <= bc.range.second) right = true;
            }
        }
        if (!left) {
            for(size_t k = 1; k<=3; ++k){
                b[k][j] = b[k][j] + a[k][j];
                a[k][j] = 0.0;
            }
        }
        if (!right) {
            for(size_t k = Nx-4; k<=Nx-2; ++k){
                b[k][j] = b[k][j] + c[k][j];
                c[k][j] = 0.0;
            }
        }
    }
    


    }


void bcCN2_T(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
             std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
             const Mesh& mesh, const SolverConfig& config, double time){

    size_t Nx = mesh.xc3.size();
    size_t Ny = mesh.yc3.size();

    // Apply only Dirichlet inlet BCs from boundary struct
    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            if (bc.side == "bottom") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) {
                        double val = bc.TFunction(mesh.xc3[i], 0.0, time);
                        for(size_t k = 1; k<=3; ++k){
                            b[i][k] = b[i][k] - a[i][k];
                            d[i][k] = d[i][k] - 2 * val * a[i][k];
                            a[i][k] = 0.0;
                        }
                    }
            }
            if (bc.side == "top") {
                for (size_t i = 0; i < Nx; ++i)
                    if (mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) {
                        double val = bc.TFunction(mesh.xc3[i], config.Ly, time);
                        for(size_t k = Ny-4; k<=Ny-2; ++k){
                            b[i][k] = b[i][k] - c[i][k];
                            d[i][k] = d[i][k] - 2 * val * c[i][k];
                            c[i][k] = 0.0;
                        }
                    }
            }

        }
    }

    // Walls and outlets use mirrored fallback condition
    for (size_t i = 0; i < Nx; ++i) {
        bool bottom = false, top = false;
        for (const auto& bc : config.boundaries) {
            if (bc.type == "inlet") {
                if (bc.side == "bottom" && mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) bottom = true;
                if (bc.side == "top" && mesh.xc3[i] >= bc.range.first && mesh.xc3[i] <= bc.range.second) top = true;
            }
        }
        if (!bottom) {
            for(size_t k = 1; k<=3; ++k){
                b[i][k] = b[i][k] + a[i][k];
                a[i][k] = 0.0;
            }
        }
        if (!top) {
            for(size_t k = Ny-4; k<=Ny-2; ++k){
                b[i][k] = b[i][k] + c[i][k];
                c[i][k] = 0.0;
            }
        }
    }

    }


void bcGS(std::vector<std::vector<double>>& phi){
    int M = phi.size();       // total rows including ghost cells
    int N = phi[0].size();    // total cols including ghost cells

    // Left boundary (i = 0)
    for (int j = 0; j < N; ++j) {
        phi[0][j] = phi[1][j];
    }

    // Right boundary (i = M-1)
    for (int j = 0; j < N; ++j) {
        phi[M-1][j] = phi[M-2][j];
    }

    // Bottom boundary (j = 0)
    for (int i = 0; i < M; ++i) {
        phi[i][0] = phi[i][1];
    }

    // Top boundary (j = N-1)
    for (int i = 0; i < M; ++i) {
        phi[i][N-1] = phi[i][N-2];
    }
}

// Zero-Neumann ghost update (flat, row-major, with ghosts)
void bcGS_flat(std::vector<double>& phi, int nx, int ny) {
    // top and bottom ghost rows
    for (int i = 0; i <= nx + 1; ++i) {
        phi[idx(i, 0,      nx)] = phi[idx(i, 1,      nx)];   // j = 0   <- 1
        phi[idx(i, ny + 1, nx)] = phi[idx(i, ny,     nx)];   // j = ny+1<- ny
    }
    // left and right ghost columns
    for (int j = 0; j <= ny + 1; ++j) {
        phi[idx(0,      j, nx)] = phi[idx(1,      j, nx)];   // i = 0   <- 1
        phi[idx(nx + 1, j, nx)] = phi[idx(nx,     j, nx)];   // i = nx+1<- nx
    }
    // corners are set consistently by the two passes (redundant writes are fine)
}



void bcGhost_u(std::vector<std::vector<double>>& u, const Mesh& mesh, const SolverConfig& config, double time) {
    size_t Nx = mesh.xf.size();
    size_t Ny = mesh.yc.size();

    for (const auto& bc : config.boundaries) {
        const auto& side = bc.side;
        const auto& range = bc.range;

        for (size_t i = 0; i < Nx; ++i) {
            double x = mesh.xf[i];

            if (side == "bottom" && x >= range.first && x <= range.second) {
                if (bc.type == "inlet") {
                    double u_bc = bc.uFunction(x, 0.0, time);
                    u[i][0] = 2.0 * u_bc - u[i][1];
                } else if (bc.type == "outlet") {
                    u[i][0] = u[i][1];  // simple central copy for ghost cell
                }
            }

            if (side == "top" && x >= range.first && x <= range.second) {
                if (bc.type == "inlet") {
                    double u_bc = bc.uFunction(x, config.Ly, time);
                    u[i][Ny - 1] = 2.0 * u_bc - u[i][Ny - 2];
                } else if (bc.type == "outlet") {
                    u[i][Ny - 1] = u[i][Ny - 2];
                }
            }
        }
    }

    // Wall fallback
    for (size_t i = 0; i < Nx; ++i) {
        double x = mesh.xf[i];
        bool bottomFound = false;
        bool topFound = false;

        for (const auto& bc : config.boundaries) {
            if (x >= bc.range.first && x <= bc.range.second) {
                if (bc.side == "bottom") bottomFound = true;
                if (bc.side == "top") topFound = true;
            }
        }

        if (!bottomFound) u[i][0] = -u[i][1];
        if (!topFound) u[i][Ny - 1] = -u[i][Ny - 2];
    }
}


void bcGhost_v(std::vector<std::vector<double>>& v, const Mesh& mesh, const SolverConfig& config, double time) {
    size_t Nx = mesh.xc.size();
    size_t Ny = mesh.yf.size();

    for (const auto& bc : config.boundaries) {
        const auto& side = bc.side;
        const auto& range = bc.range;

        for (size_t j = 0; j < Ny; ++j) {
            double y = mesh.yf[j];

            if (side == "left" && y >= range.first && y <= range.second) {
                if (bc.type == "inlet") {
                    double v_bc = bc.vFunction(0.0, y, time);
                    v[0][j] = 2.0 * v_bc - v[1][j];
                } else if (bc.type == "outlet") {
                    v[0][j] = v[1][j];
                }
            }

            if (side == "right" && y >= range.first && y <= range.second) {
                if (bc.type == "inlet") {
                    double v_bc = bc.vFunction(config.Lx, y, time);
                    v[Nx - 1][j] = 2.0 * v_bc - v[Nx - 2][j];
                } else if (bc.type == "outlet") {
                    v[Nx - 1][j] = v[Nx - 2][j];
                }
            }
        }
    }

    // Wall fallback
    for (size_t j = 0; j < Ny; ++j) {
        double y = mesh.yf[j];
        bool leftFound = false;
        bool rightFound = false;

        for (const auto& bc : config.boundaries) {
            if (y >= bc.range.first && y <= bc.range.second) {
                if (bc.side == "left") leftFound = true;
                if (bc.side == "right") rightFound = true;
            }
        }

        if (!leftFound)  v[0][j] = -v[1][j];
        if (!rightFound) v[Nx - 1][j] = -v[Nx - 2][j];
    }
}