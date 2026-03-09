#ifndef DIVERGENCE_PART2_DATA_H
#define DIVERGENCE_PART2_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel divergence part 2 kernel execution.
struct DivergencePart2Work {
    int nm;      ///< Mesh index
    double dt;   ///< Time step

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    DivergencePart2Work(int nm_, double dt_,
                        std::shared_ptr<MeshData> md)
        : nm(nm_), dt(dt_), originalMeshData(md) {}
};

#endif // DIVERGENCE_PART2_DATA_H
