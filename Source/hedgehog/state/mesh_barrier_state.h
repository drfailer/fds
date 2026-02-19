#ifndef MESH_BARRIER_STATE_H
#define MESH_BARRIER_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Barrier state that collects N mesh tokens, calls fds_mesh_exchange(code),
/// then re-emits all N tokens. This replaces the synchronization points in the
/// original FDS main loop where all meshes must complete before exchanging data.
class MeshBarrierState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    MeshBarrierState(int nmeshes, int exchangeCode)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes), exchangeCode_(exchangeCode) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // All meshes arrived at the barrier - perform the exchange
            fds_mesh_exchange(exchangeCode_);
            // Re-emit all tokens
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
        }
    }

private:
    int nmeshes_;
    int exchangeCode_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // MESH_BARRIER_STATE_H
