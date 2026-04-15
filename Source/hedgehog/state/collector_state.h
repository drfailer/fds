#ifndef COLLECTOR_STATE_H
#define COLLECTOR_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Task that collects N MeshData<> tokens and emits a single BarrierData
/// containing all of them. Pure data-flow control with no computation —
/// the actual work is done by a downstream barrier task.
///
/// Runs on a single thread. Meshes are placed directly at their correct
/// position using NM as the index, avoiding any sorting overhead.
class CollectorTask : public hh::AbstractTask<1, MeshData<>, BarrierData> {
public:
    /// @param nmeshes    Number of distinct meshes (array size).
    /// @param name       Task name for graph visualization.
    /// @param totalExpected  Number of tokens to collect before firing.
    ///                       Defaults to nmeshes.  Set to 2*nmeshes for
    ///                       fork-join collection (second token per mesh
    ///                       overwrites first; count gates the fire).
    explicit CollectorTask(int nmeshes, std::string name = "Collector",
                           int totalExpected = 0)
        : hh::AbstractTask<1, MeshData<>, BarrierData>(std::move(name), 1),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()),
          totalExpected_(totalExpected > 0 ? totalExpected : nmeshes) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;
        if (count_ == totalExpected_) {
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
    int totalExpected_;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

#endif // COLLECTOR_STATE_H
