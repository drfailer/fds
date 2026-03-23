#ifndef PIPELINE_FORK2_STATE_H
#define PIPELINE_FORK2_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Fork state for Corrector Fork 2: RADIATION || DIV_P1.
/// Collects N MeshData tokens, runs InitDivIntegrals (zero DSUM/PSUM/USUM),
/// then emits MeshData. Hedgehog multicasts to both branches.
class PipelineFork2State
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PipelineFork2State(int nmeshes)
        : nmeshes_(nmeshes) {
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
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // PIPELINE_FORK2_STATE_H
