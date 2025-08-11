#pragma once
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// Internal smoothers matched to backends
enum class InternalSmoother { GS, RBGS, Jacobi };

inline InternalSmoother smoother_for_backend(EllipticBackend be) {
    switch (be) {
        case EllipticBackend::Serial: return InternalSmoother::GS;
        case EllipticBackend::OpenMP: return InternalSmoother::RBGS;
        case EllipticBackend::CUDA:   return InternalSmoother::Jacobi;
    }
    return InternalSmoother::GS; // default fallback
}

// Thread selection helper
inline int choose_omp_threads(int requested) {
#ifdef _OPENMP
    int maxT = omp_get_max_threads();
    if (requested <= 0) return maxT;
    if (requested > maxT) return maxT;  // clamp
    return requested;
#else
    (void)requested;
    return 1;
#endif
}