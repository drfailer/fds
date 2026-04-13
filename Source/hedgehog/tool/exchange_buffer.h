#ifndef EXCHANGE_BUFFER_H
#define EXCHANGE_BUFFER_H

#include <vector>
#include <cassert>
#include "mesh_dependency_graph.h"
#include "../fds_fortran_interface.h"

/// Pre-allocated flat storage for exchange data, indexed by (source NM, dest NOM).
///
/// Allocates buffers for ALL (source, dest) pairs where dest is a local mesh,
/// including cross-rank sources.  This enables both same-rank and cross-rank
/// exchange data to flow through the same ExchangeBuffer.
///
/// Each ExchangeBuffer is specific to one exchange code (e.g. CODE 5 for flux).
/// Buffer sizes are determined at construction time by calling the unified
/// dispatch function fds_exchange_slab_size_recv(code, nom, nm).
///
/// Thread safety: different (NM, NOM) pairs occupy non-overlapping regions in
/// the flat storage.  Pipeline ordering (deps gate emits only after all writes
/// complete) ensures that write and pull never access the same entry concurrently.
class ExchangeBuffer {
public:
    /// Build the buffer from the pre-computed dependency graph for a given
    /// exchange code. Allocates storage for all (source, dest) pairs where
    /// dest is local.
    ExchangeBuffer(const MeshDependencyGraph &depGraph, int exchangeCode)
        : lower_(depGraph.lowerMesh()),
          upper_(depGraph.upperMesh()),
          totalMeshes_(depGraph.totalMeshes()) {

        int n = totalMeshes_ + 1;  // 1-based indexing
        offsets_.resize(static_cast<size_t>(n * n), -1);
        sizes_.resize(static_cast<size_t>(n * n), 0);
        recvSources_.resize(static_cast<size_t>(n));

        int totalSize = 0;

        // Iterate over local dest meshes and their recv deps (all ranks).
        // This allocates buffers for same-rank AND cross-rank source meshes.
        for (int nom = lower_; nom <= upper_; ++nom) {
            const auto &deps = depGraph.recvDeps(nom);
            for (int i = 0; i < totalMeshes_; ++i) {
                if (!deps.contains(static_cast<size_t>(i))) continue;
                int nm = i + 1;  // 1-based mesh index
                int sz = fds_exchange_slab_size_recv(exchangeCode, nom, nm);
                if (sz <= 0) continue;
                offsets_[index(nm, nom)] = totalSize;
                sizes_[index(nm, nom)] = sz;
                totalSize += sz;
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

    /// Get the list of all source meshes that send TO mesh nm.
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
