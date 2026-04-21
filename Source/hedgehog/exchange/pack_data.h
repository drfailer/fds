#ifndef PACK_DATA_H
#define PACK_DATA_H

#include <vector>
#include <ostream>
#include <hedgehog_comm.h>
#include "../fds_fortran_interface.h"

/// Packed halo data for one (source mesh, dest mesh) pair.
///
/// Template parameter K selects the exchange code at compile time for
/// type-based Hedgehog routing. The Fortran pack/unpack dispatch uses
/// exchangeCode() which maps K → FDS CODE integer.
///
/// Lifecycle:
///   PackTask creates one PackData per remote neighbor:
///     reset(sourceNm, destNom) → fillBuffer() → emit to CommunicatorTask
///   CommunicatorTask sends/receives via MPI:
///     pack() → send; package() → receive → unpack(Package&&)
///   CopyTask downstream:
///     unpackIntoOMesh() → emit downstream
template<MeshState K>
struct PackData {
    int sourceNm = 0;
    int destNom  = 0;
    int slabSize = 0;
    std::vector<double> buffer;

    PackData() = default;

    static constexpr int exchangeCode() {
        if constexpr (K == MeshState::MeshExch1) return 1;
        else if constexpr (K == MeshState::MeshExch2) return 2;
        else if constexpr (K == MeshState::MeshExch3) return 3;
        else if constexpr (K == MeshState::MeshExch4) return 4;
        else if constexpr (K == MeshState::MeshExch5) return 5;
        else if constexpr (K == MeshState::MeshExch6) return 6;
        else if constexpr (K == MeshState::MeshExch7) return 7;
        else if constexpr (K == MeshState::PreSolveExch) return 5;
        else if constexpr (K == MeshState::PostSolveExch) return 5;
        else return 0;
    }

    /// Configure for a (source, dest) pair. Source must be local.
    void reset(int src, int dst) {
        sourceNm = src;
        destNom  = dst;
        slabSize = fds_exchange_slab_size(exchangeCode(), src, dst);
        buffer.resize(static_cast<size_t>(slabSize));
    }

    /// Pack source mesh halo data into buffer (sender side).
    void fillBuffer() {
        fds_exchange_push_slab(exchangeCode(), sourceNm, destNom,
                               buffer.data(), slabSize);
    }

    /// Unpack received buffer into dest mesh OMESH arrays (receiver side).
    void unpackIntoOMesh() {
        fds_exchange_pull_slab_recv(exchangeCode(), destNom, sourceNm,
                                    buffer.data(), slabSize);
    }

    // --- MPI serialization for CommunicatorTask ---

    struct Header {
        int sourceNm = 0;
        int destNom  = 0;
        int slabSize = 0;
    };

    hh::comm::Package pack() {
        header_ = Header{sourceNm, destNom, slabSize};
        hh::comm::Package pkg;
        pkg.data.emplace_back(
            reinterpret_cast<char *>(&header_), sizeof(Header));
        pkg.data.emplace_back(
            reinterpret_cast<char *>(buffer.data()),
            static_cast<size_t>(slabSize) * sizeof(double));
        return pkg;
    }

    hh::comm::Package package() {
        if (buffer.size() < static_cast<size_t>(maxSlabSize_)) {
            buffer.resize(static_cast<size_t>(maxSlabSize_));
        }
        hh::comm::Package pkg;
        pkg.data.emplace_back(
            reinterpret_cast<char *>(&header_), sizeof(Header));
        pkg.data.emplace_back(
            reinterpret_cast<char *>(buffer.data()),
            buffer.size() * sizeof(double));
        return pkg;
    }

    void unpack(hh::comm::Package &&) {
        sourceNm = header_.sourceNm;
        destNom  = header_.destNom;
        slabSize = header_.slabSize;
        buffer.resize(static_cast<size_t>(slabSize));
    }

    void cleanMemory() {
        sourceNm = 0;
        destNom = 0;
        slabSize = 0;
    }

    friend std::ostream &operator<<(std::ostream &os, const PackData &d) {
        os << "PackData<" << exchangeCode() << ">{src=" << d.sourceNm
           << ", dst=" << d.destNom << ", size=" << d.slabSize << "}";
        return os;
    }

    static inline int maxSlabSize_ = 0;

private:
    Header header_{};
};

#endif // PACK_DATA_H
