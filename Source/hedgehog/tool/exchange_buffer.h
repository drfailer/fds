#ifndef EXCHANGE_BUFFER_H
#define EXCHANGE_BUFFER_H

#include <vector>
#include <cassert>
#include "mesh_dependency_graph.h"
#include "../fds_fortran_interface.h"

/// Pre-allocated flat storage for exchange data, indexed by (source NM, dest NOM).
///
/// Each MeshExchangeGraph instance owns its own ExchangeBuffer, so concurrent
/// exchanges (e.g. pre-solve and post-solve) never share buffers.
///
/// Thread safety: different (NM, NOM) pairs occupy non-overlapping regions in
/// the flat storage.  Pipeline ordering (deps gate emits only after all pushes
/// complete) ensures that push and pull never access the same entry concurrently.
///
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename Strategy>
class ExchangeBuffer {
public:
    /// Build the buffer from the pre-computed dependency graph.
    /// Allocates storage for all same-rank (source, dest) pairs.
    ExchangeBuffer(const MeshDependencyGraph &depGraph)
        : lower_(depGraph.lowerMesh()),
          upper_(depGraph.upperMesh()),
          totalMeshes_(depGraph.totalMeshes()) {

        int n = totalMeshes_ + 1;  // 1-based indexing
        offsets_.resize(static_cast<size_t>(n * n), -1);
        sizes_.resize(static_cast<size_t>(n * n), 0);
        recvSources_.resize(static_cast<size_t>(n));

        int myRank = fds_mesh_process(lower_);
        int totalSize = 0;

        // Forward pass: compute offsets and sizes for each (source, dest) pair
        for (int nm = lower_; nm <= upper_; ++nm) {
            for (int nom : depGraph.sendTargets(nm)) {
                if (fds_mesh_process(nom) != myRank) continue;
                int sz = Strategy::bufferSize(nm, nom);
                if (sz <= 0) continue;
                offsets_[index(nm, nom)] = totalSize;
                sizes_[index(nm, nom)] = sz;
                totalSize += sz;
                // Build reverse mapping: nom receives from nm
                recvSources_[static_cast<size_t>(nom)].push_back(nm);
            }
        }

        storage_.resize(static_cast<size_t>(totalSize), 0.0);
    }

    /// Get pointer to the buffer for the (source nm, dest nom) pair.
    [[nodiscard]] double *buffer(int nm, int nom) {
        int off = offsets_[index(nm, nom)];
        assert(off >= 0 && "No buffer allocated for this (nm, nom) pair");
        return storage_.data() + off;
    }

    /// Get buffer size (doubles) for the (source nm, dest nom) pair.
    [[nodiscard]] int bufferSize(int nm, int nom) const {
        return sizes_[index(nm, nom)];
    }

    /// Get the list of same-rank source meshes that send TO mesh nm.
    /// Used by the pull task to know which buffers to read from.
    [[nodiscard]] const std::vector<int> &recvSources(int nm) const {
        return recvSources_[static_cast<size_t>(nm)];
    }

    /// Total storage size in doubles.
    [[nodiscard]] size_t totalStorageSize() const { return storage_.size(); }

private:
    [[nodiscard]] size_t index(int nm, int nom) const {
        return static_cast<size_t>(nm) * static_cast<size_t>(totalMeshes_ + 1)
             + static_cast<size_t>(nom);
    }

    int lower_;
    int upper_;
    int totalMeshes_;
    std::vector<double> storage_;
    std::vector<int> offsets_;
    std::vector<int> sizes_;
    std::vector<std::vector<int>> recvSources_;
};

#endif // EXCHANGE_BUFFER_H
