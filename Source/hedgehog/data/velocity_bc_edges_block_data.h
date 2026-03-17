#ifndef VELOCITY_BC_EDGES_BLOCK_DATA_H
#define VELOCITY_BC_EDGES_BLOCK_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for block-decomposed VelocityBC edge processing.
/// Each token represents a K-range sub-block of a single mesh,
/// carrying the edge processing parameters and DRAG_UVWMAX accumulator.
struct VelocityBCEdgesBlockWork {
    int nm;                ///< Mesh index (1-based)
    int k1;                ///< K-range start (1-based inclusive)
    int k2;                ///< K-range end (1-based inclusive)
    double t;              ///< Simulation time
    int applyToEstimated;  ///< 1 for predictor, 0 for corrector
    int totalBlocks;       ///< Total blocks for this mesh (for reassembly)
    double dragUvwMax;     ///< Per-block DRAG_UVWMAX (output from kernel, reduced in collector)

    std::shared_ptr<MeshData> originalMeshData;  ///< Parent token for downstream routing

    VelocityBCEdgesBlockWork(int nm_, int k1_, int k2_, double t_,
                              int est_, int total_,
                              std::shared_ptr<MeshData> md)
        : nm(nm_), k1(k1_), k2(k2_), t(t_),
          applyToEstimated(est_), totalBlocks(total_),
          dragUvwMax(0.0), originalMeshData(std::move(md)) {}
};

#endif // VELOCITY_BC_EDGES_BLOCK_DATA_H
