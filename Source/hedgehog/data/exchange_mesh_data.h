#ifndef EXCHANGE_MESH_DATA_H
#define EXCHANGE_MESH_DATA_H

#include <vector>
#include <memory>
#include <hedgehog_comm.h>
#include "../fds_fortran_interface.h"

/// Lightweight dep-satisfaction signal from WriteBufferTask to Gate.
struct ExchangeDepSignal {
    int sourceNm;
    int destNom;
    int roundId;  ///< Exchange round for double-buffered Gate slot selection
};

/// Slab exchange data carrying arrays for one (source, dest) pair.
///
/// Used by the MeshExchangeGraph pipeline:
///   FanOutTask -> CommunicatorTask -> WriteBufferTask -> Gate -> PullTask
///
/// The exchange code (which data to push/pull) is carried in the header
/// and dispatched at runtime through unified Fortran functions.
///
/// For same-rank pairs, the CommunicatorTask loops back (addResult, no MPI).
/// For cross-rank pairs, pack/unpack/package handle MPI serialization.
struct ExchangeMeshData {
    struct Header {
        int sourceNm = 0;
        int destNom  = 0;
        int slabSize = 0;
        int exchangeCode = 0;
        int roundId = 0;  ///< Exchange round for double-buffered state selection
    } header_;

    std::vector<double> buffer_;

    /// Global max slab size across all supported codes.
    /// Set once at startup from MeshExchangeGraph constructor.
    static inline int globalMaxSlabSize_ = 0;

    ExchangeMeshData() = default;

    /// Reconfigure for a new (source, dest, code) triple.
    /// Uses sender-side size function (sourceNm must be local).
    void reset(int sourceNm, int destNom, int exchangeCode) {
        header_.sourceNm = sourceNm;
        header_.destNom = destNom;
        header_.exchangeCode = exchangeCode;
        header_.slabSize = fds_exchange_slab_size(exchangeCode, sourceNm, destNom);
        buffer_.resize(static_cast<size_t>(header_.slabSize));
    }

    /// Copy slab data FROM source Fortran arrays INTO buffer_.
    void push() {
        fds_exchange_push_slab(header_.exchangeCode, header_.sourceNm,
                               header_.destNom, buffer_.data(), header_.slabSize);
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
        buffer_.resize(static_cast<size_t>(globalMaxSlabSize_));
    }

    /// Called after receive: trim buffer to actual slab size from header.
    void unpack(hh::comm::Package &&) {
        buffer_.resize(static_cast<size_t>(header_.slabSize));
    }

    /// Called by pool on release: reset header for reuse.
    void cleanMemory() { header_ = {}; }
};

#endif // EXCHANGE_MESH_DATA_H
