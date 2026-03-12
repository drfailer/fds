#ifndef COLLECTOR_STATE_H
#define COLLECTOR_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Generic barrier state that collects N MeshData tokens and emits a single
/// BarrierData containing all of them. This is pure data-flow control with
/// no computation — the actual work is done by a downstream barrier task.
///
/// Meshes are placed directly at their correct position using NM as the index,
/// avoiding any sorting overhead.
class CollectorState : public hh::AbstractState<1, MeshData, BarrierData> {
public:
    explicit CollectorState(int nmeshes)
        : hh::AbstractState<1, MeshData, BarrierData>(),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;
        if (count_ == nmeshes_) {
            auto bd = std::make_shared<BarrierData>();
            bd->meshes = std::move(collected_);
            collected_.resize(nmeshes_, nullptr);
            count_ = 0;
            this->addResult(bd);
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // COLLECTOR_STATE_H
