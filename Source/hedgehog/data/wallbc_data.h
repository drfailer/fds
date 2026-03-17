#ifndef WALLBC_DATA_H
#define WALLBC_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel WALL_BC cell processing kernel execution.
/// Bundles WALL_BC_PROCESS_CELLS_KERNEL (processes ~90% of wall cells
/// without cross-mesh dependencies).
struct WallBCWork {
    int nm;             ///< Mesh index
    double t;           ///< Simulation time
    double dt;          ///< Time step
    double dt_bc;       ///< Boundary condition time step (Fortran-computed)
    int call_ht_1d;     ///< Flag to call 1-D heat transfer (0=false, 1=true)

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    WallBCWork(int nm_, double t_, double dt_, double dt_bc_,
               int call_ht_1d_, std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), dt_bc(dt_bc_),
          call_ht_1d(call_ht_1d_), originalMeshData(md) {}
};

/// Work token for block-decomposed WALL_BC kernel execution.
/// Extends WallBCWork with K-range block information for intra-mesh parallelism.
struct WallBCBlockWork {
    int nm;             ///< Mesh index
    int k1;             ///< K-range start (1-based inclusive)
    int k2;             ///< K-range end (1-based inclusive)
    double t;           ///< Simulation time
    double dt;          ///< Time step
    double dt_bc;       ///< Boundary condition time step
    int call_ht_1d;     ///< Flag to call 1-D heat transfer (0=false, 1=true)
    int totalBlocks;    ///< Total number of blocks for this mesh (for reassembly)

    /// Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    WallBCBlockWork(int nm_, int k1_, int k2_, double t_, double dt_,
                    double dt_bc_, int call_ht_1d_, int totalBlocks_,
                    std::shared_ptr<MeshData> md)
        : nm(nm_), k1(k1_), k2(k2_), t(t_), dt(dt_), dt_bc(dt_bc_),
          call_ht_1d(call_ht_1d_), totalBlocks(totalBlocks_),
          originalMeshData(std::move(md)) {}
};

#endif // WALLBC_DATA_H
