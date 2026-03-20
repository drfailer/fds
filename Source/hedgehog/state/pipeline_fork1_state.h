#ifndef PIPELINE_FORK1_STATE_H
#define PIPELINE_FORK1_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/pipeline_fork1_data.h"

/// Fork state for Corrector Fork 1: VFLUX || COMBUSTION.
/// Receives MeshData from MeshExchange(4) and dispatches both branches
/// concurrently for each mesh. No collection needed — each MeshData
/// immediately spawns one VFluxWork + one CombWork.
class PipelineFork1State
    : public hh::AbstractState<1, MeshData, Fork1VFluxWork, Fork1CombWork> {
public:
    PipelineFork1State() = default;

    void execute(std::shared_ptr<MeshData> data) override {
        this->addResult(std::make_shared<Fork1VFluxWork>(data));
        this->addResult(std::make_shared<Fork1CombWork>(data));
    }
};

/// Join state for Corrector Fork 1.
/// Collects VFluxResult and CombResult per mesh. When both arrive for
/// the same mesh, emits the original MeshData downstream.
class PipelineJoin1State
    : public hh::AbstractState<2, Fork1VFluxResult, Fork1CombResult, MeshData> {
public:
    PipelineJoin1State() = default;

    void execute(std::shared_ptr<Fork1VFluxResult> result) override {
        int nm = result->meshData->nm;
        vfluxDone_[nm] = result->meshData;
        tryEmit(nm);
    }

    void execute(std::shared_ptr<Fork1CombResult> result) override {
        int nm = result->meshData->nm;
        combDone_[nm] = result->meshData;
        tryEmit(nm);
    }

private:
    void tryEmit(int nm) {
        auto vIt = vfluxDone_.find(nm);
        auto cIt = combDone_.find(nm);
        if (vIt != vfluxDone_.end() && cIt != combDone_.end()) {
            this->addResult(vIt->second);
            vfluxDone_.erase(vIt);
            combDone_.erase(cIt);
        }
    }

    std::unordered_map<int, std::shared_ptr<MeshData>> vfluxDone_;
    std::unordered_map<int, std::shared_ptr<MeshData>> combDone_;
};

#endif // PIPELINE_FORK1_STATE_H
