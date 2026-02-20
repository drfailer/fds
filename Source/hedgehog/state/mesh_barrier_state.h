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

/// Barrier state for COMBUSTION_LOAD_BALANCED. Collects N mesh tokens,
/// calls fds_combustion(t, dt), then re-emits all tokens.
/// The original code (main.f90:837) calls COMBUSTION_LOAD_BALANCED(T,DT)
/// as a global routine, not a mesh exchange.
class CombustionBarrierState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    CombustionBarrierState(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            fds_combustion(collected_[0]->t, collected_[0]->dt);
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

/// Barrier state for HVAC_CALC. Collects N mesh tokens,
/// calls fds_hvac_calc(t, dt, first), then re-emits all tokens.
class HvacBarrierState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    HvacBarrierState(int nmeshes, int first)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes), first_(first) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            fds_hvac_calc(collected_[0]->t, collected_[0]->dt, first_);
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
        }
    }

private:
    int nmeshes_;
    int first_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Passthrough barrier: collects N mesh tokens, then re-emits all without
/// calling any global routine. Used to prevent concurrent Fortran calls
/// between adjacent tasks that share no barrier, since Hedgehog runs each
/// task on its own thread and module-level Fortran pointers (via POINT_TO_MESH)
/// are not thread-safe.
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
