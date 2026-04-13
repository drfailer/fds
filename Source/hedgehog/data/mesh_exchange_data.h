#ifndef MESH_EXCHANGE_DATA_H
#define MESH_EXCHANGE_DATA_H

#include <memory>
#include <vector>
#include <ostream>
#include "mesh_data.h"

/// Token emitted by MeshDependenciesManagerState to FluxExchangeTask.
///
/// Carries the mesh to be exchanged and its full list of same-rank neighbors.
/// The exchange task uses pull-only copies: each mesh pulls from all its
/// neighbors, writing only to its own OMESH buffers.
struct MeshExchangeData {
    std::shared_ptr<MeshData<>> mesh;   ///< The mesh to exchange
    std::vector<int> neighbors;       ///< Same-rank neighbors to pull from (1-based)

    MeshExchangeData() = default;
    MeshExchangeData(std::shared_ptr<MeshData<>> m, std::vector<int> nbrs)
        : mesh(std::move(m)), neighbors(std::move(nbrs)) {}

    friend std::ostream &operator<<(std::ostream &os, const MeshExchangeData &d) {
        os << "MeshExchangeData{nm=" << d.mesh->nm << ", neighbors=[";
        for (size_t i = 0; i < d.neighbors.size(); ++i) {
            if (i > 0) os << ",";
            os << d.neighbors[i];
        }
        os << "]}";
        return os;
    }
};

#endif // MESH_EXCHANGE_DATA_H
