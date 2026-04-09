#ifndef EXCHANGE_FLUX_MESH_DATA_H
#define EXCHANGE_FLUX_MESH_DATA_H

#include <vector>
#include <memory>
#include <hedgehog_comm.h>

/// Lightweight dep-satisfaction signal from WriteBufferTask to Gate.
struct ExchangeDepSignal {
    int sourceNm;
    int destNom;
};

/// Slab exchange data carrying flux arrays for one (source, dest) pair.
///
/// Used by the MeshExchangeGraph pipeline:
///   FanOutTask -> CommunicatorTask -> WriteBufferTask -> Gate -> PullTask
///
/// For same-rank pairs, the CommunicatorTask loops back (addResult, no MPI).
/// For cross-rank pairs, pack/unpack/package handle MPI serialization.
///
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename Strategy>
struct ExchangeFluxMeshData {
    struct Header {
        int sourceNm = 0;
        int destNom  = 0;
        int slabSize = 0;
    } header_;

    std::vector<double> buffer_;

    ExchangeFluxMeshData() = default;

    /// Reconfigure for a new (source, dest) pair.
    void reset(int sourceNm, int destNom) {
        header_.sourceNm = sourceNm;
        header_.destNom  = destNom;
        header_.slabSize = Strategy::bufferSize(sourceNm, destNom);
        buffer_.resize(static_cast<size_t>(header_.slabSize));
    }

    /// Copy slab data FROM source Fortran arrays INTO buffer_.
    void push() {
        Strategy::pushToBuffer(header_.sourceNm, header_.destNom,
                               buffer_.data(), header_.slabSize);
    }

    // --- MPI serialization (2 buffers: header + slab data) ---

    /// Pack for sending: returns pointers to header and slab (actual size).
    hh::comm::Package pack() {
        return {{
            {reinterpret_cast<char *>(&header_), sizeof(Header)},
            {reinterpret_cast<char *>(buffer_.data()),
             static_cast<size_t>(header_.slabSize) * sizeof(double)}
        }};
    }

    /// Package for receiving: returns pointers to header and slab (max size).
    hh::comm::Package package() {
        return {{
            {reinterpret_cast<char *>(&header_), sizeof(Header)},
            {reinterpret_cast<char *>(buffer_.data()),
             buffer_.size() * sizeof(double)}
        }};
    }

    /// Called before receive: allocate max buffer to receive any slab size.
    void preRecv() {
        buffer_.resize(static_cast<size_t>(Strategy::maxBufferSize()));
    }

    /// Called after receive: trim buffer to actual slab size from header.
    void unpack(hh::comm::Package &&) {
        buffer_.resize(static_cast<size_t>(header_.slabSize));
    }

    /// Called by pool on release: reset header for reuse.
    void cleanMemory() { header_ = {}; }
};

#endif // EXCHANGE_FLUX_MESH_DATA_H
