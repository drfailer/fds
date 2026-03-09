#ifndef VELOCITY_PREDICTOR_DATA_H
#define VELOCITY_PREDICTOR_DATA_H

#include "mesh_data.h"
#include <memory>
#include <ostream>

/// Work token for parallel velocity predictor kernel execution.
/// This token flows through the velocity predictor sub-graph, carrying
/// mesh-specific parameters for kernel execution while preserving the
/// original MeshData for downstream graph routing.
struct VelocityPredictorWork {
    int nm;          ///< Mesh index (1-based, Fortran convention)
    double t;        ///< Simulation time
    double dt;       ///< Time step

    /// Original MeshData token to be emitted downstream after kernel completion
    std::shared_ptr<MeshData> originalMeshData;

    VelocityPredictorWork(int nm_, double t_, double dt_,
                          std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}

    friend std::ostream &operator<<(std::ostream &os, const VelocityPredictorWork &w) {
        os << "VelocityPredictorWork{nm=" << w.nm
           << ", t=" << w.t << ", dt=" << w.dt << "}";
        return os;
    }
};

#endif // VELOCITY_PREDICTOR_DATA_H
