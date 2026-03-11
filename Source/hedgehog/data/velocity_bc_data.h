#ifndef VELOCITY_BC_DATA_H
#define VELOCITY_BC_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel VELOCITY_BC_PROCESS_EDGES_KERNEL execution.
struct VelocityBCWork {
    int nm;               ///< Mesh index
    double t;             ///< Simulation time
    double dt;            ///< Time step
    int applyToEstimated; ///< 1 = estimated (predictor), 0 = actual (corrector)

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    VelocityBCWork(int nm_, double t_, double dt_, int estimated_,
                   std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), applyToEstimated(estimated_),
          originalMeshData(md) {}
};

#endif // VELOCITY_BC_DATA_H
