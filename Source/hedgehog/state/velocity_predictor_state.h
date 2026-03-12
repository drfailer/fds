#ifndef VELOCITY_PREDICTOR_STATE_H
#define VELOCITY_PREDICTOR_STATE_H

#include <algorithm>
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
        : nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        results_.push_back(data);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) { return a->nm < b->nm; });

            for (auto &md : results_) {
                fds_cc_project_velocity(md->nm, md->dt, 0);  // STORE=.FALSE.
                fds_wall_velocity_no_gradh(md->nm, md->dt, 0);  // STORE=.FALSE.
                fds_check_stability_kernel_only(md->nm, md->t + md->dt, md->dt);
            }

            for (auto &md : results_) {
                this->addResult(md);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> results_;
};

#endif // VELOCITY_PREDICTOR_STATE_H
