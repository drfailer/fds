#ifndef COLLECTOR_STATE_H
#define COLLECTOR_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"

/// Generic barrier state that collects N MeshData tokens and emits a single
/// BarrierData containing all of them. This is pure data-flow control with
/// no computation — the actual work is done by a downstream barrier task.
class CollectorState : public hh::AbstractState<1, MeshData, BarrierData> {
public:
    explicit CollectorState(int nmeshes)
        : hh::AbstractState<1, MeshData, BarrierData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sort by mesh index to match original FDS ordering (ascending NM)
            std::sort(collected_.begin(), collected_.end(),
                      [](const auto &a, const auto &b) { return a->nm < b->nm; });
            auto bd = std::make_shared<BarrierData>();
            bd->meshes = std::move(collected_);
            collected_ = {};
            collected_.reserve(nmeshes_);
            this->addResult(bd);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // COLLECTOR_STATE_H
