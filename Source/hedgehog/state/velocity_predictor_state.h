#ifndef VELOCITY_PREDICTOR_STATE_H
#define VELOCITY_PREDICTOR_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Collector for CC_PROJECT_VELOCITY after velocity predictor kernel (CC_IBM only).
///
/// Gathers all N kernel results, runs sequential CC_PROJECT_VELOCITY(STORE=FALSE),
/// then emits N MeshData tokens.
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
