#include "solver.h"
#include "parabolic.h"
#include "utils.h"
#include <iostream>
#include "multi_grid.h"
#include "solution.h"
#include "bc.h"
#include "hyperbolic.h"
#include "IB.h"
#include "MG_flat.h"

using Array2D = std::vector<std::vector<double>>;

void runSolver(Solution& sol, const Mesh& mesh, SolverConfig& config) {
    double t = 0.0;
    auto u_old = sol.u;
    auto v_old = sol.v;
    Array2D Hu, Hv, Hu_nm1, Hv_nm1, HT;
    double finalTime = config.totalTime.back();  

    for (double tf : config.totalTime) {
        bool flag = false;
        //int test = 1;

        while (t < tf) {
        //while (test <= 1) {

            flag = calcDt(t, tf, config, sol);  // sets config.dt and flag
            

            sol.Qu = zeros(config.M + 1,config.N + 2);
            sol.Qv = zeros(config.M + 2,config.N + 1);
            sol.QT = zeros(config.M + 6,config.N + 6);

            if (config.enableRotors)
                calcSourceIB(sol,mesh,config,t);

            auto QT1 = sol.QT; //store for temp ADI
            auto QT2 = sol.QT; //store for temp ADI 

            if (config.enableHyperbolicSolver) {


                //solve hyperbolic forcing terms for u and v
                std::tie(Hu_nm1, Hv_nm1) = hyperbolic_uv_2D(u_old,v_old,config);
                std::tie(Hu, Hv) = hyperbolic_uv_2D(sol.u,sol.v,config);

                // save old velocities
                u_old = sol.u;
                v_old = sol.v;

                // add to existing forcing terms if present
                auto Term1 = scalarMultiply(Hu,3.0/2.0);
                auto Term2 = scalarMultiply(Hu_nm1,1.0/2.0);
                sol.Qu = elementwiseAdd(sol.Qu,Term1);
                sol.Qu = elementwiseSubtract(sol.Qu,Term2);

                Term1 = scalarMultiply(Hv,3.0/2.0);
                Term2 = scalarMultiply(Hv_nm1,1.0/2.0);
                sol.Qv = elementwiseAdd(sol.Qv,Term1);
                sol.Qv = elementwiseSubtract(sol.Qv,Term2);
            } 

            // Apply step 1 of ADI
            parabolic_CN1_u(sol, mesh, config, t);
            parabolic_CN1_v(sol, mesh, config, t);

            sol.Qu = zeros(config.M + 1,config.N + 2);
            sol.Qv = zeros(config.M + 2,config.N + 1);
            sol.QT = zeros(config.M + 6,config.N + 6);
            
            if (config.enableRotors)
                calcSourceIB(sol,mesh,config,t+config.dt/2.0);

            

            if (config.enableHyperbolicSolver) {

                // add to existing forcing terms if present
                auto Term1 = scalarMultiply(Hu,3.0/2.0);
                auto Term2 = scalarMultiply(Hu_nm1,1.0/2.0);
                sol.Qu = elementwiseAdd(sol.Qu,Term1);
                sol.Qu = elementwiseSubtract(sol.Qu,Term2);

                Term1 = scalarMultiply(Hv,3.0/2.0);
                Term2 = scalarMultiply(Hv_nm1,1.0/2.0);
                sol.Qv = elementwiseAdd(sol.Qv,Term1);
                sol.Qv = elementwiseSubtract(sol.Qv,Term2);
            }


            // Apply step 2 of ADI
            parabolic_CN2_u(sol, mesh, config, t);
            parabolic_CN2_v(sol, mesh, config, t);

            if (config.enableHyperbolicSolver) {  //calculate T hyperbolic term using solve velocity
                HT = hyperbolic_T_WENO_2D(sol, u_old, v_old, t, config, mesh);
                sol.QT = elementwiseAdd(QT1,HT);
            }


            //perform ADI for T 
            parabolic_CN1_T(sol, mesh, config, t);

            if (config.enableRotors)
                calcSourceIB(sol,mesh,config,t+config.dt/2.0);
            
            sol.QT = zeros(config.M + 6,config.N + 6); 
            QT2 = sol.QT; //store for temp ADI


            if (config.enableHyperbolicSolver) {  //calculate T hyperbolic term using solve velocity
                sol.QT = elementwiseAdd(QT2,HT);
            }

            parabolic_CN2_T(sol, mesh, config, t);



            if (config.enableEllipticSolver) { // elliptic solver section

                correctOutletMassConservation(sol, mesh, config);

                auto divV = calcDivV(sol, mesh);
                auto f = scalarMultiply(divV, 1.0 / config.dt);

                int nIterMax = 100;
                double epsilon = 1e-6;
                int nu1 = 1;
                int nu2 = 1;
                int nCoarse = 3;

                poisson_solve(sol,f,config,nIterMax,epsilon,nu1,nu2,nCoarse);
                //sol.phi = poissonSolver(sol.phi, f, mesh, nIterMax, epsilon);

                bcGS(sol.phi);
                projectV(sol, mesh, config, t);
                auto divV_post = calcDivV(sol, mesh);


            }


            t += config.dt;
            printTimeProgress(t, finalTime);

            //test += 1;

            if (flag){
                std::cout << std::endl << "Flag time reached t =  "<< t << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

Solution initializeSolution(const Mesh& mesh, const SolverConfig& config) {
    Solution sol;

    // Initialize Solution Fields
    sol.u.resize(config.M + 1, std::vector<double>(config.N + 2, 0.0));    // xf × yc
    sol.v.resize(config.M + 2, std::vector<double>(config.N + 1, 0.0));    // xc × yf
    sol.T.resize(config.M + 6, std::vector<double>(config.N + 6, 1.0));    // xc3 × yc3
    sol.phi.resize(config.M + 2, std::vector<double>(config.N + 2, 0.0));  // xc × yc
    sol.Qu.resize(config.M + 1, std::vector<double>(config.N + 2, 0.0));
    sol.Qv.resize(config.M + 2, std::vector<double>(config.N + 1, 0.0));
    sol.QT.resize(config.M + 6, std::vector<double>(config.N + 6, 0.0));

    double t = 0.0;

    // Apply Velocity and Temperature Boundary Conditions
    bc_u(sol.u, mesh, config, t);
    bc_v(sol.v, mesh, config, t);
    bc_T(sol.T, mesh, config, t);

    // Elliptic Solver Intialization (mass conservation)
    if (config.enableEllipticSolver) {
        correctOutletMassConservation(sol, mesh, config);
        
        
        auto divV = calcDivV(sol, mesh);
        auto f = scalarMultiply(divV, 1.0 / config.dt);

        int nIterMax = 100;
        double epsilon = 1e-8;
        int nu1 = 1;
        int nu2 = 1;
        int nCoarse = 3;

        //sol.phi = poissonSolver(sol.phi, f, mesh, nIterMax, epsilon);

        poisson_solve(sol,f,config,nIterMax,epsilon,nu1,nu2,nCoarse);

        bcGS(sol.phi);
        projectV(sol, mesh, config, t);
        
        
    }

    return sol;
}