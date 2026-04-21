#ifndef MESH_DIM_H
#define MESH_DIM_H

/// Target mesh dimensions for automatic mesh re-decomposition.
/// When set (all > 0), each user-configured mesh is split into
/// sub-meshes of approximately this many cells per dimension.
struct MeshDim {
    int i = 0;  ///< Target cells in I direction (0 = no splitting)
    int j = 0;  ///< Target cells in J direction (0 = no splitting)
    int k = 0;  ///< Target cells in K direction (0 = no splitting)

    /// Returns true if mesh re-decomposition is requested.
    [[nodiscard]] bool enabled() const { return i > 0 && j > 0 && k > 0; }
};

#endif // MESH_DIM_H
