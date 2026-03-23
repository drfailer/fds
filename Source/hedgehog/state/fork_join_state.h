#ifndef FORK_JOIN_STATE_H
#define FORK_JOIN_STATE_H

#include <hedgehog/hedgehog.h>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"

/// Join state for fork-join patterns with MeshData.
/// Counts arrivals per mesh from multiple branches and emits after all arrive.
class ForkJoinState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit ForkJoinState(int numBranches = 2) : numBranches_(numBranches) {}

    void execute(std::shared_ptr<MeshData> data) override {
        int nm = data->nm;
        counts_[nm]++;
        if (counts_[nm] == numBranches_) {
            counts_.erase(nm);
            this->addResult(data);
        }
    }

private:
    int numBranches_;
    std::unordered_map<int, int> counts_;
};

/// Join state for fork-join patterns with BarrierData.
/// Counts arrivals from multiple branches and emits after all arrive.
class BarrierJoinState : public hh::AbstractState<1, BarrierData, BarrierData> {
public:
    explicit BarrierJoinState(int numBranches = 2) : numBranches_(numBranches) {}

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
