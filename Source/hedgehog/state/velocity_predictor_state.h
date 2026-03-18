#ifndef VELOCITY_PREDICTOR_STATE_H
#define VELOCITY_PREDICTOR_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Collector for CC_IBM post-processing after velocity predictor kernel.
///
/// Gathers all N kernel results, then for each mesh runs the sequential
/// operations that must follow VELOCITY_PREDICTOR_KERNEL for CC_IBM:
///   1. CC_PROJECT_VELOCITY(STORE=FALSE) — project velocities onto cut-cells
///   2. WALL_VELOCITY_NO_GRADH(STORE=FALSE) — fix wall velocities for sparse solvers
///   3. CHECK_STABILITY_KERNEL — compute CFL-limited DT_NEW
///
/// This matches the ordering in velo.f90 VELOCITY_PREDICTOR (lines 574-612).
class VelocityPredictorCCCollector
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit VelocityPredictorCCCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            for (auto &md : collected_) {
                fds_cc_project_velocity(md->nm, md->dt, 0);  // STORE=.FALSE.
                fds_wall_velocity_no_gradh_kernel(md->nm, md->dt, 0, 1);  // store=0, predictor=1
                fds_check_stability_kernel_only(md->nm, md->t + md->dt, md->dt);
            }

            for (auto &md : collected_) {
                this->addResult(md);
            }

            std::fill(collected_.begin(), collected_.end(), nullptr);
            count_ = 0;
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // VELOCITY_PREDICTOR_STATE_H
