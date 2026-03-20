#ifndef PRED_FORK_STATE_H
#define PRED_FORK_STATE_H

#include <hedgehog/hedgehog.h>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../data/pred_fork_data.h"

/// Fork state for Predictor: (VFLUX+PMOM) || (WALLBC+DIV_P1_early).
/// Receives MeshData from DIV_P1_prefork and dispatches both branches
/// concurrently for each mesh.
class PredForkState
    : public hh::AbstractState<1, MeshData, PredForkVFluxWork, PredForkDivWork> {
public:
    PredForkState() = default;

    void execute(std::shared_ptr<MeshData> data) override {
        this->addResult(std::make_shared<PredForkVFluxWork>(data));
        this->addResult(std::make_shared<PredForkDivWork>(data));
    }
};

/// Join state for Predictor Fork.
/// Collects VFluxResult and DivResult per mesh. When both arrive for
/// the same mesh, emits the original MeshData downstream.
class PredJoinState
    : public hh::AbstractState<2, PredForkVFluxResult, PredForkDivResult, MeshData> {
public:
    PredJoinState() = default;

    void execute(std::shared_ptr<PredForkVFluxResult> result) override {
        int nm = result->meshData->nm;
        vfluxDone_[nm] = result->meshData;
        tryEmit(nm);
    }

    void execute(std::shared_ptr<PredForkDivResult> result) override {
        int nm = result->meshData->nm;
        divDone_[nm] = result->meshData;
        tryEmit(nm);
    }

private:
    void tryEmit(int nm) {
        auto vIt = vfluxDone_.find(nm);
        auto dIt = divDone_.find(nm);
        if (vIt != vfluxDone_.end() && dIt != divDone_.end()) {
            this->addResult(vIt->second);
            vfluxDone_.erase(vIt);
            divDone_.erase(dIt);
        }
    }

    std::unordered_map<int, std::shared_ptr<MeshData>> vfluxDone_;
    std::unordered_map<int, std::shared_ptr<MeshData>> divDone_;
};

#endif // PRED_FORK_STATE_H
