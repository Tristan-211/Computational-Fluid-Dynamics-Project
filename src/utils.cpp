#include "utils.h"
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <algorithm> 
#include "config.h"
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include "solution.h"

// linspace function similar to MATLAB vec = linspace(start,end,num_points)
std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;

    if (num == 0) {
        return result;
    }
    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (num - 1);
    result.reserve(num);

    for (int i = 0; i < num; ++i) {
        result.push_back(start + step * i);
    }

    return result;
}

std::vector<double> solveTriDiag(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d) {

    size_t n = a.size();
    if (b.size() != n || c.size() != n || d.size() != n) {
        throw std::runtime_error("Vectors a, b, c, d must have the same length");
    }

    std::vector<double> cp(c);  // upper diagonal
    std::vector<double> bp(b);  // main diagonal
    std::vector<double> dp(d);  // RHS

    // Forward elimination
    for (size_t i = 1; i < n; ++i) {
        double m = a[i] / bp[i - 1];
        bp[i] -= m * cp[i - 1];
        dp[i] -= m * dp[i - 1];
    }

    // Back substitution
    std::vector<double> x(n);
    x[n - 1] = dp[n - 1] / bp[n - 1];
    for (size_t i = n - 1; i-- > 0;) {
        x[i] = (dp[i] - cp[i] * x[i + 1]) / bp[i];
    }

    return x;
}




void saveMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    file << std::setprecision(16) << std::fixed;
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            double value = row[j];
            if (value == 0.0) value = 0.0;  // clean -0.0 to +0.0
            file << value;
            if (j < row.size() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Saved matrix to " << filename << std::endl;
}


bool calcDt(double t, double outputTime, SolverConfig& config, const Solution& sol) {
    // Calculate maximum velocity magnitudes (used only if hyperbolic solver is enabled)
    double dtConv = std::numeric_limits<double>::max();

    if (config.enableHyperbolicSolver) {
        double umax = maxAbs(sol.u);
        double vmax = maxAbs(sol.v);

        double ueq = 2.0 * umax + vmax;
        double veq = 2.0 * vmax + umax;

        double dtu = config.h / (ueq + 1e-14);  // avoid divide-by-zero
        double dtv = config.h / (veq + 1e-14);

        dtConv = std::min(dtu, dtv);
    }

    // Elliptic constraints (viscous terms)
    double dtDiff1 = 0.25 * config.h * config.h / (1.0 / config.Re);
    double dtDiff2 = 0.25 * config.h * config.h / (1.0 / (config.Re * config.Pr));

    // Minimum of all
    double dtMin = std::min({dtConv, dtDiff1, dtDiff2});
    config.dt = config.CFL * dtMin;

    // Adjust timestep if it would overshoot output time
    bool outputFlag = false;
    if (t < outputTime && (t + config.dt) >= outputTime) {
        config.dt = outputTime - t;
        outputFlag = true;
    }

    return outputFlag;
}


void meshgrid(const std::vector<double>& x,
              const std::vector<double>& y,
              std::vector<std::vector<double>>& X,
              std::vector<std::vector<double>>& Y)
{
    size_t Ny = y.size();
    size_t Nx = x.size();

    X.resize(Ny, std::vector<double>(Nx));
    Y.resize(Ny, std::vector<double>(Nx));

    for (size_t i = 0; i < Ny; ++i) {
        for (size_t j = 0; j < Nx; ++j) {
            X[i][j] = x[j];  // replicate x across rows
            Y[i][j] = y[i];  // replicate y across columns
        }
    }
}

std::vector<std::vector<double>> reshape1DTo2D(const std::vector<double>& flat, size_t rows, size_t cols) {
    if (flat.size() != rows * cols) throw std::runtime_error("reshape1DTo2D: size mismatch");
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            result[i][j] = flat[i * cols + j];  // row-major layout
    return result;
}


std::vector<std::vector<double>> elementwiseAdd(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B) {
    
    size_t rows = A.size(), cols = A[0].size();
    std::vector<std::vector<double>> C(rows, std::vector<double>(cols));
    
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

std::vector<std::vector<double>> elementwiseSubtract(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B) {

    size_t rows = A.size(), cols = A[0].size();
    std::vector<std::vector<double>> C(rows, std::vector<double>(cols));
    
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

std::vector<std::vector<double>> elementwiseMultiply(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B) {
    
    size_t rows = A.size(), cols = A[0].size();
    std::vector<std::vector<double>> C(rows, std::vector<double>(cols));
    
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            C[i][j] = A[i][j] * B[i][j];
    return C;
}

std::vector<std::vector<double>> elementwiseDivide(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B) {
    
    size_t rows = A.size(), cols = A[0].size();
    std::vector<std::vector<double>> C(rows, std::vector<double>(cols));
    
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            C[i][j] = A[i][j] / B[i][j];
    return C;
}

std::vector<std::vector<double>> scalarMultiply(const std::vector<std::vector<double>>& A,double scalar) {
    
    size_t rows = A.size(), cols = A[0].size();
    std::vector<std::vector<double>> C(rows, std::vector<double>(cols));
    
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            C[i][j] = A[i][j] * scalar;
    return C;
}


std::vector<std::vector<bool>> makeMask(const std::vector<std::vector<double>>& XR,
                                        const std::vector<std::vector<double>>& YR,
                                        double IBlh,
                                        double IBwh) {
    size_t M = XR.size();
    size_t N = XR[0].size();
    std::vector<std::vector<bool>> mask(M, std::vector<bool>(N, false));

    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            if (std::fabs(XR[i][j]) <= IBlh && std::fabs(YR[i][j]) < IBwh)
                mask[i][j] = true;

    return mask;
}


void applyMaskSet(std::vector<std::vector<double>>& target,
                  const std::vector<std::vector<bool>>& mask,
                  double value) {
    size_t M = target.size();
    size_t N = target[0].size();
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            if (mask[i][j])
                target[i][j] = value;
}


void checkRotorBounds(const std::vector<Rotor>& rotors, const Mesh& mesh) {
    double xmin = mesh.xf.front();
    double xmax = mesh.xf.back();
    double ymin = mesh.yf.front();
    double ymax = mesh.yf.back();

    for (size_t i = 0; i < rotors.size(); ++i) {
        const Rotor& r = rotors[i];

        double x0 = r.center[0];
        double y0 = r.center[1];

        double left   = x0 - r.halfLength;
        double right  = x0 + r.halfLength;
        double bottom = y0 - r.halfWidth;
        double top    = y0 + r.halfWidth;

        if (left < xmin || right > xmax || bottom < ymin || top > ymax) {
            std::ostringstream msg;
            msg << "Rotor " << i << " exceeds mesh bounds:\n"
                << "  Rotor bounds: ["
                << left << ", " << right << "] x ["
                << bottom << ", " << top << "]\n"
                << "  Mesh bounds: ["
                << xmin << ", " << xmax << "] x ["
                << ymin << ", " << ymax << "]";
            throw std::runtime_error(msg.str());
        }
    }
}

std::vector<std::vector<double>> Residual(
    const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<double>>& f,
    double h)
{
    int M = phi.size();
    int N = phi[0].size();
    double h2 = h * h;

    std::vector<std::vector<double>> r(M, std::vector<double>(N, 0.0));

    for (int i = 1; i < M - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            r[i][j] = f[i][j] - (
                (phi[i+1][j] - 2.0 * phi[i][j] + phi[i-1][j]) / h2 +
                (phi[i][j+1] - 2.0 * phi[i][j] + phi[i][j-1]) / h2
            );
        }
    }

    return r;
}

double RelativeResidualNorm(
    const std::vector<std::vector<double>>& phi,
    const std::vector<std::vector<double>>& f,
    double h)
{
    int M = phi.size();
    int N = phi[0].size();
    
    std::vector<std::vector<double>> r = Residual(phi, f, h);

    double max_r = 0.0;
    double max_f = 0.0;

    for (int i = 1; i < M - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            max_r = std::max(max_r, std::abs(r[i][j]));
            max_f = std::max(max_f, std::abs(f[i][j]));
        }
    }

    return (max_f > 0.0) ? max_r / max_f : 0.0; // avoid division by zero
}


double maxAbs(const std::vector<std::vector<double>>& mat) {
    double max_val = 0.0;
    for (const auto& row : mat) {
        for (double val : row) {
            max_val = std::max(max_val, std::abs(val));
        }
    }
    return max_val;
}