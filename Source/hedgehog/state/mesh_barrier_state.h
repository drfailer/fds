#ifndef MESH_BARRIER_STATE_H
#define MESH_BARRIER_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"

/// Passthrough barrier: collects N mesh tokens, then re-emits all without
/// calling any global routine. Used to prevent concurrent Fortran calls
/// between adjacent tasks that share no barrier, since Hedgehog runs each
/// task on its own thread and module-level Fortran pointers (via POINT_TO_MESH)
/// are not thread-safe.
///
/// This is pure data-flow control — no computation.
class PassthroughBarrierState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PassthroughBarrierState(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // MESH_BARRIER_STATE_H
