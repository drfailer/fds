#ifndef CORR_DIV_PART1_DATA_H
#define CORR_DIV_PART1_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel corrector divergence part 1 kernel execution.
struct CorrDivPart1Work {
    int nm;      ///< Mesh index
    double t;    ///< Simulation time
    double dt;   ///< Time step

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    CorrDivPart1Work(int nm_, double t_, double dt_,
                     std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // CORR_DIV_PART1_DATA_H
