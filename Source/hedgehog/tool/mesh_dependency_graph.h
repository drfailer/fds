#ifndef MESH_DEPENDENCY_GRAPH_H
#define MESH_DEPENDENCY_GRAPH_H

#include <vector>
#include <string>
#include <sstream>
#include "dyn_bitset.h"
#include "../fds_fortran_interface.h"

/// Pre-computed mesh exchange dependency graph.
///
/// For each local mesh NM, stores:
///   - recvDeps_[nm]: bitset of meshes that send TO nm (nm's prerequisites)
///   - sendTargets_[nm]: list of meshes that nm sends TO
///   - sameRankRecvDeps_[nm]: bitset of same-rank senders (subset of recvDeps_)
///
/// Built once at init time from Fortran queries.  Used by
/// ExchangeOrchestratorState to track when a mesh's dependencies are satisfied.
///
/// Mesh indices: 1-based (Fortran convention).  Internal arrays are sized
/// [0..totalMeshes] with index 0 unused.
class MeshDependencyGraph {
public:
    /// Build the dependency graph for locally-owned meshes.
    /// Must be called after FDS initialization (INITIALIZE_MESH_EXCHANGE_1).
    MeshDependencyGraph(int lowerMesh, int upperMesh)
        : lower_(lowerMesh), upper_(upperMesh) {
        int totalMeshes = fds_get_total_meshes();
        totalMeshes_ = totalMeshes;
        int myRank = fds_mesh_process(lowerMesh);

        recvDeps_.resize(totalMeshes + 1, DynBitset(static_cast<size_t>(totalMeshes)));
        sameRankRecvDeps_.resize(totalMeshes + 1, DynBitset(static_cast<size_t>(totalMeshes)));
        sendTargets_.resize(totalMeshes + 1);

        for (int nm = lowerMesh; nm <= upperMesh; ++nm) {
            // Who sends TO nm?
            int nRecv = fds_exchange_recv_dep_count(nm);
            for (int i = 1; i <= nRecv; ++i) {
                int nom = fds_exchange_recv_dep_mesh(nm, i);
                recvDeps_[nm].set(static_cast<size_t>(nom - 1)); // 0-based bit index
                if (fds_mesh_process(nom) == myRank) {
                    sameRankRecvDeps_[nm].set(static_cast<size_t>(nom - 1));
                }
            }

            // Who does nm send TO?
            int nSend = fds_exchange_send_dep_count(nm);
            sendTargets_[nm].reserve(static_cast<size_t>(nSend));
            for (int i = 1; i <= nSend; ++i) {
                int nom = fds_exchange_send_dep_mesh(nm, i);
                sendTargets_[nm].push_back(nom);
            }
        }
    }

    /// Bitset of all meshes that send data TO nm.
    [[nodiscard]] const DynBitset &recvDeps(int nm) const { return recvDeps_[nm]; }

    /// Bitset of same-rank meshes that send data TO nm.
    [[nodiscard]] const DynBitset &sameRankRecvDeps(int nm) const { return sameRankRecvDeps_[nm]; }

    /// List of meshes that nm sends data TO.
    [[nodiscard]] const std::vector<int> &sendTargets(int nm) const { return sendTargets_[nm]; }

    /// Total number of meshes (global).
    [[nodiscard]] int totalMeshes() const { return totalMeshes_; }

    [[nodiscard]] int lowerMesh() const { return lower_; }
    [[nodiscard]] int upperMesh() const { return upper_; }

    /// Debug: print the dependency graph.
    [[nodiscard]] std::string dump() const {
        std::ostringstream oss;
        oss << "MeshDependencyGraph (meshes " << lower_ << ".." << upper_ << "):\n";
        for (int nm = lower_; nm <= upper_; ++nm) {
            oss << "  M" << nm << " recv_deps={";
            bool first = true;
            for (int i = 0; i < totalMeshes_; ++i) {
                if (recvDeps_[nm].contains(static_cast<size_t>(i))) {
                    if (!first) oss << ",";
                    oss << (i + 1);
                    first = false;
                }
            }
            oss << "} send_to={";
            first = true;
            for (int t : sendTargets_[nm]) {
                if (!first) oss << ",";
                oss << t;
                first = false;
            }
            oss << "}\n";
        }
        return oss.str();
    }

private:
    int lower_;
    int upper_;
    int totalMeshes_;
    std::vector<DynBitset> recvDeps_;
    std::vector<DynBitset> sameRankRecvDeps_;
    std::vector<std::vector<int>> sendTargets_;
};

#endif // MESH_DEPENDENCY_GRAPH_H
