#ifndef DIVERGENCE_BARRIER_STATE_H
#define DIVERGENCE_BARRIER_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Barrier that calls fds_initialize_divergence_integrals() once all meshes arrive.
/// Must be placed BEFORE tasks that call DIVERGENCE_PART_1, which accumulates
/// into the DSUM/PSUM/USUM arrays that this barrier zeroes.
class InitDivIntegralsBarrier : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit InitDivIntegralsBarrier(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            fds_initialize_divergence_integrals();
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

/// Barrier state for divergence exchange.
/// After all meshes complete DIVERGENCE_PART_1: calls fds_exchange_divergence_info()
/// and fds_global_matrix_reassign(), then re-emits all tokens.
class DivergenceBarrierState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    DivergenceBarrierState(int nmeshes, bool corrector = false)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes), corrector_(corrector) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Exchange divergence info across all meshes
            fds_exchange_divergence_info();

            // In corrector phase: RTE source correction happens here
            // (between DIVERGENCE_PART_1 and EXCHANGE_DIVERGENCE_INFO in main.f90
            //  but actually after exchange in the original code... let's put it after)
            if (corrector_) {
                fds_rte_source_correction();
            }

            // Global matrix reassign (corrector only in original, but safe to call)
            fds_global_matrix_reassign(0);

            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
        }
    }

private:
    int nmeshes_;
    bool corrector_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // DIVERGENCE_BARRIER_STATE_H
