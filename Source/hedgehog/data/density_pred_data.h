#ifndef DENSITY_PRED_DATA_H
#define DENSITY_PRED_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel density predictor kernel execution.
struct DensityPredWork {
    int nm;      ///< Mesh index
    double t;    ///< Simulation time
    double dt;   ///< Time step

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    DensityPredWork(int nm_, double t_, double dt_,
                    std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // DENSITY_PRED_DATA_H
