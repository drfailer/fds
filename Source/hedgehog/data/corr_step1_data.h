#ifndef CORR_STEP1_DATA_H
#define CORR_STEP1_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel corrector step 1 kernel execution.
/// Bundles COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES + DENSITY kernels.
struct CorrStep1Work {
    int nm;      ///< Mesh index
    double t;    ///< Simulation time
    double dt;   ///< Time step

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    CorrStep1Work(int nm_, double t_, double dt_,
                  std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // CORR_STEP1_DATA_H
