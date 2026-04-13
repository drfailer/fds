#ifndef FORK_JOIN_STATE_H
#define FORK_JOIN_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Join task for fork-join patterns with MeshData<>.
/// Counts arrivals per mesh from multiple branches and emits after all arrive.
///
/// Runs on a single thread.
class ForkJoinTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit ForkJoinTask(int nmeshes, int numBranches = 2,
                          int numDownstreamThreads = 1,
                          std::string name = "ForkJoin")
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(std::move(name), 1),
          nmeshes_(nmeshes), numBranches_(numBranches),
          numDownstreamThreads_(numDownstreamThreads),
          nmOffset_(fds_get_lower_mesh_index()) {
        counts_.resize(nmeshes, 0);
        readyList_.reserve(numDownstreamThreads);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        int idx = data->nm - nmOffset_;
        counts_[idx]++;
        if (counts_[idx] == numBranches_) {
            counts_[idx] = 0;
            completedCount_++;
            readyList_.push_back(data);

            if (static_cast<int>(readyList_.size()) == numDownstreamThreads_
                || completedCount_ == nmeshes_) {
                this->batchAddResult(readyList_);
                readyList_.clear();
                if (completedCount_ == nmeshes_) {
                    completedCount_ = 0;
                }
            }
        }
    }

private:
    int nmeshes_;
    int numBranches_;
    int numDownstreamThreads_;
    int nmOffset_;
    int completedCount_ = 0;
    std::vector<int> counts_;
    std::vector<std::shared_ptr<MeshData<>>> readyList_;
};

/// Join task for fork-join patterns with BarrierData.
/// Counts arrivals from multiple branches and emits after all arrive.
///
/// Runs on a single thread.
class BarrierJoinTask : public hh::AbstractTask<1, BarrierData, BarrierData> {
public:
    explicit BarrierJoinTask(int numBranches = 2, std::string name = "BarrierJoin")
        : hh::AbstractTask<1, BarrierData, BarrierData>(std::move(name), 1),
          numBranches_(numBranches) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        count_++;
        if (!lastData_) lastData_ = data;
        if (count_ == numBranches_) {
            count_ = 0;
            this->addResult(lastData_);
            lastData_ = nullptr;
        }
    }

private:
    int numBranches_;
    int count_ = 0;
    std::shared_ptr<BarrierData> lastData_;
};

#endif // FORK_JOIN_STATE_H
