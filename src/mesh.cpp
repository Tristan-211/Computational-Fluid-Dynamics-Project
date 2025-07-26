#include "mesh.h"
#include "utils.h"
#include <cmath>

Mesh generateMesh(double Lx, double Ly, int M, int N) {
    Mesh mesh;

    // creating u staggered mesh
    mesh.xf = linspace(0, Lx, M + 1);
    mesh.h = std::fabs(mesh.xf[1] - mesh.xf[0]);

    // Build yc mesh: from -h/2 to Ly + h/2, step h
    int num_yc = N + 2;  // N intervals + 1, plus extra for staggered
    double start_yc = -mesh.h / 2;
    double end_yc = Ly + mesh.h / 2;

    mesh.yc = linspace(start_yc, end_yc, num_yc);


    // creating v staggered mesh
    mesh.yf = linspace(0,Ly,N+1);

    // Build xc mesh: from -h/2 to Lx + h/2, step h
    int num_xc = M + 2;  // M intervals + 1, plus extra for staggered
    double start_xc = -mesh.h / 2;
    double end_xc = Lx + mesh.h / 2;

    mesh.xc = linspace(start_xc, end_xc, num_xc);



    // creating cell centred mesh with 3 ghost cell layers
    mesh.xc3 = linspace(0-2.5*mesh.h,Lx+2.5*mesh.h,M+6);
    mesh.yc3 = linspace(0-2.5*mesh.h,Ly+2.5*mesh.h,N+6);


    return mesh;
}