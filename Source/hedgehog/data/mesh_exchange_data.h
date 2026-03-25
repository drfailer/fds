#ifndef MESH_EXCHANGE_DATA_H
#define MESH_EXCHANGE_DATA_H

#include <memory>
#include <vector>
#include <ostream>
#include "mesh_data.h"

/// Token emitted by MeshDependenciesManagerState to MeshExchangeTask.
///
/// Carries the mesh to be exchanged and the list of neighbors that still
/// need bidirectional copies.  Neighbors already exchanged by a previous
/// mesh are excluded from this list.
struct MeshExchangeData {
    std::shared_ptr<MeshData> mesh;   ///< The mesh to exchange
    std::vector<int> neighbors;       ///< Same-rank neighbors needing exchange (1-based)

    MeshExchangeData() = default;
    MeshExchangeData(std::shared_ptr<MeshData> m, std::vector<int> nbrs)
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
