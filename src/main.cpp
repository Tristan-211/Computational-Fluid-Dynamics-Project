#include <iostream>
#include <vector>
#include <iomanip>
#include "mesh.h"
#include "utils.h"
#include "config.h"
#include <cmath>
#include "solution.h"
#include "bc.h"
#include "parabolic.h"
#include "solver.h"
#include "multi_grid.h"
#include "hyperbolic.h"
#include "IB.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {


    // Solver/Domain parameters
    SolverConfig config;
    config.M = 64;
    config.N = 3 * config.M / 4;
    config.Lx = 4.0;
    config.Ly = 3.0;
    config.Re = 200.0;
    config.Pr = 3.0;
    config.totalTime = {3};
    config.CFL = 0.4;
    config.h = config.Lx/config.M;
    config.dt = 0.01; //value for intialization only calc dt at start of solver loop
    config.enableEllipticSolver = true;
    config.enableHyperbolicSolver = true;
    config.enableRotors = true;


    // Boundary Condition Structs 
    BoundaryCondition inlet1;
    inlet1.side = "bottom"; 
    inlet1.range = {0.75, 1.25};
    inlet1.type = "inlet";
    inlet1.uFunction = [](double x, double, double) {
        double z = (x - 0.75) / 0.5;
        return 2.0 * cos(2 * M_PI / 3) * (6 * z * (1 - z));
    };
    inlet1.vFunction = [](double x, double, double) {
        double z = (x - 0.75) / 0.5;
        return 2.0 * sin(2 * M_PI / 3) * (6 * z * (1 - z));  
    };
    inlet1.TFunction = [](double, double, double t) {
        return 1.5 + 0.5 * sin(2 * M_PI * t); 
    };
    config.boundaries.push_back(inlet1);

    BoundaryCondition inlet2;
    inlet2.side = "top"; 
    inlet2.range = {0.5, 1.25};
    inlet2.type = "inlet";
    inlet2.uFunction = [](double x, double, double) {
        double z = (x - 0.5) / 0.75;
        return 1.5 * cos(4 * M_PI / 3) * (6 * z * (1 - z));
    };
    inlet2.vFunction = [](double x, double, double) {
        double z = (x - 0.5) / 0.75;
        return 1.5 * sin(4 * M_PI / 3) * (6 * z * (1 - z));  
    };
    inlet2.TFunction = [](double, double, double t) {
        return 0.5 + 0.25 * sin(5 * M_PI / 2 * t); 
    };
    config.boundaries.push_back(inlet2);


    BoundaryCondition outlet1;
    outlet1.side = "right";
    outlet1.range = {0.75, 2.25};
    outlet1.type = "outlet";
    outlet1.uFunction = [](double x, double, double) {
        return 0.0;
    };
    outlet1.vFunction = [](double x, double, double) {
        return 0.0;
    };
    outlet1.TFunction = [](double x, double, double) {
        return 1.0;
    };
    config.boundaries.push_back(outlet1);


    bool hasInlet = false;
    bool hasOutlet = false;

    for (const auto& bc : config.boundaries) {
        if (bc.type == "inlet") {
            hasInlet = true;
        }
        if (bc.type == "outlet") {
            hasOutlet = true;
        }
    }

    if (hasInlet && !hasOutlet) {
        std::cerr << "ERROR: You have at least one inlet but no outlet defined. Please add an outlet boundary condition." << std::endl;
        return 1;  // exit with error code
    }

    std::cout << "Solver configuration set up with "
    << config.boundaries.size() << " boundary conditions." << std::endl;

    //Rotors 
    // Inside main.cpp or a separate setup function for now
    Rotor rotor1, rotor2;

    // Set center, angular velocity, initial angle, and temperature
    rotor1.center = {1.0, 1.0};
    rotor2.center = {1.25, 2.0};
    rotor1.omega = 0.5 * 2.0 * M_PI;
    rotor2.omega = -0.5 * 2.0 * M_PI;
    rotor1.theta0 = M_PI / 6.0;
    rotor2.theta0 = -M_PI / 4.0;
    rotor1.temp = 3.0;
    rotor2.temp = 0.0;

    // Set size
    rotor1.width = rotor2.width = 0.25;
    rotor1.length = rotor2.length = 0.75;

    // Cache half-dimensions
    rotor1.halfWidth = rotor1.width / 2.0;
    rotor1.halfLength = rotor1.length / 2.0;
    rotor2.halfWidth = rotor2.width / 2.0;
    rotor2.halfLength = rotor2.length / 2.0;

    config.rotors.push_back(rotor1);
    config.rotors.push_back(rotor2);


    // Create Mesh
    Mesh mesh = generateMesh(config.Lx, config.Ly, config.M, config.N);



    //check rotors 
    checkRotorBounds(config.rotors,mesh);



    // Initialize Solution Vars
    Solution sol = initializeSolution(mesh,config);
    initIBMeshU(mesh,config,sol);



    //saveMatrixToFile(sol.u, "u.csv");
    //saveMatrixToFile(sol.v, "v.csv");
    //saveMatrixToFile(sol.T, "T.csv");




    runSolver(sol, mesh, config);

    


    //calcSourceIB(sol,mesh,config,t);

    //saveMatrixToFile(sol.Qu, "Qu.csv");
    //saveMatrixToFile(sol.Qv, "Qv.csv");
    //saveMatrixToFile(sol.QT, "QT.csv");

    //saveMatrixToFile(sol.u, "u.csv");
    //saveMatrixToFile(sol.v, "v.csv");
    //saveMatrixToFile(sol.T, "T.csv");





    
    /*
    double t = 0.0;
 
    bool flag = calcDt(t,1.0, config, sol);
    auto u_old = sol.u;
    auto v_old = sol.v;


    auto [Hu_nm1, Hv_nm1] = hyperbolic_uv_2D(u_old,v_old,config);
    auto [Hu, Hv] = hyperbolic_uv_2D(sol.u,sol.v,config);


    auto Term1 = scalarMultiply(Hu,3.0/2.0);
    auto Term2 = scalarMultiply(Hu_nm1,1.0/2.0);
    sol.Qu = elementwiseAdd(sol.Qu,Term1);
    sol.Qu = elementwiseSubtract(sol.Qu,Term2);

    Term1 = scalarMultiply(Hv,3.0/2.0);
    Term2 = scalarMultiply(Hv_nm1,1.0/2.0);
    sol.Qv = elementwiseAdd(sol.Qv,Term1);
    sol.Qv = elementwiseSubtract(sol.Qv,Term2);

    // Apply step 1 of ADI
    parabolic_CN1_u(sol, mesh, config, t);
    parabolic_CN1_v(sol, mesh, config, t);

    // Apply step 2 of ADI
    parabolic_CN2_u(sol, mesh, config, t);
    parabolic_CN2_v(sol, mesh, config, t);

    auto HT = hyperbolic_T_WENO_2D(sol, u_old, v_old, t, config, mesh);

    sol.QT = elementwiseAdd(sol.QT,HT);
    parabolic_CN1_T(sol, mesh, config, t);
    parabolic_CN2_T(sol, mesh, config, t);


    //saveMatrixToFile(sol.phi, "phi.csv");
    //saveMatrixToFile(sol.Qu, "Qu.csv");
    //saveMatrixToFile(sol.Qv, "Qv.csv");
    //saveMatrixToFile(HT, "HT.csv");
    */

    //saveMatrixToFile(sol.u, "u.csv");
    //saveMatrixToFile(sol.v, "v.csv");
    //saveMatrixToFile(sol.T, "T.csv");







    std::cout << std::endl;


    return 0;
}


