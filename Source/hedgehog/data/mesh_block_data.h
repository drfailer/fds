// UNUSED — Block decomposition disabled. Kept for reference.
#ifndef MESH_BLOCK_DATA_H
#define MESH_BLOCK_DATA_H

#include <memory>
#include <ostream>
#include "mesh_data.h"

/// Token representing a sub-range (block) of a mesh for intra-mesh parallelism.
/// The block covers cells [k1, k2] along the K dimension (1-based inclusive).
/// All I and J cells are included (full I,J range).
struct MeshBlockData {
    int nm;              ///< Mesh index (1-based, Fortran convention)
    int k1;              ///< K-range start (1-based inclusive)
    int k2;              ///< K-range end (1-based inclusive)
    double t;            ///< Current simulation time
    double dt;           ///< Current time step
    int phase;           ///< 0 = predictor, 1 = corrector
    int totalBlocks;     ///< Total number of blocks for this mesh (for reassembly)

    std::shared_ptr<MeshData<>> originalMeshData;  ///< Parent token for reassembly

    MeshBlockData() = default;
    MeshBlockData(int nm_, int k1_, int k2_, double t_, double dt_, int phase_,
                  int totalBlocks_, std::shared_ptr<MeshData<>> md)
        : nm(nm_), k1(k1_), k2(k2_), t(t_), dt(dt_), phase(phase_),
          totalBlocks(totalBlocks_), originalMeshData(std::move(md)) {}

    friend std::ostream &operator<<(std::ostream &os, const MeshBlockData &mb) {
        os << "MeshBlockData{nm=" << mb.nm << ", k=[" << mb.k1 << ":" << mb.k2
           << "], t=" << mb.t << ", dt=" << mb.dt << "}";
        return os;
    }
};

#endif // MESH_BLOCK_DATA_H
