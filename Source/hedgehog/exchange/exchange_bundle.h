#ifndef EXCHANGE_BUNDLE_H
#define EXCHANGE_BUNDLE_H

#include <memory>
#include <vector>
#include <ostream>
#include "../data/mesh_data.h"
#include "pack_data.h"

/// Bundle emitted by NeighborWaitState when mesh NM is released.
/// Carries the local MeshData plus all PackData buffers from cross-rank
/// neighbors (destNom == NM). CopyTask consumes this to do all data movement.
template<MeshState K>
struct ExchangeBundle {
    std::shared_ptr<MeshData<K>> mesh;
    std::vector<std::shared_ptr<PackData<K>>> remoteHalos;

    ExchangeBundle() = default;
    ExchangeBundle(std::shared_ptr<MeshData<K>> m)
        : mesh(std::move(m)) {}

    friend std::ostream &operator<<(std::ostream &os, const ExchangeBundle &b) {
        os << "ExchangeBundle<" << PackData<K>::exchangeCode()
           << ">{nm=" << b.mesh->nm
           << ", remoteHalos=" << b.remoteHalos.size() << "}";
        return os;
    }
};

#endif // EXCHANGE_BUNDLE_H
