#ifndef EXCHANGE_WORK_H
#define EXCHANGE_WORK_H

#include <memory>
#include <vector>
#include "mesh_data.h"

/// Work token emitted when all same-rank dependencies of destNM are satisfied.
///
/// The parallel ExchangeTask copies data from each sender into destNM's OMESH
/// arrays (via fds_flux_copy_neighbor_ts), then emits destNM's MeshData.
///
/// For future MPI support: crossRankSenders would trigger MPI_Wait/unpack.
struct ExchangePullWork {
    std::shared_ptr<MeshData> mesh;            ///< The destination mesh
    std::vector<int> sameRankSenders;           ///< Same-rank meshes that send TO destNM
    std::vector<int> crossRankSenders;          ///< Cross-rank meshes (future)

    ExchangePullWork(std::shared_ptr<MeshData> m,
                     std::vector<int> sameSenders,
                     std::vector<int> crossSenders)
        : mesh(std::move(m)),
          sameRankSenders(std::move(sameSenders)),
          crossRankSenders(std::move(crossSenders)) {}
};

#endif // EXCHANGE_WORK_H
