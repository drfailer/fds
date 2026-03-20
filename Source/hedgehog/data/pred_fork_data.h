#ifndef PRED_FORK_DATA_H
#define PRED_FORK_DATA_H

#include <memory>
#include "mesh_data.h"

/// Branch A work token for Predictor Fork: VFLUX + PARTICLE_MOMENTUM
struct PredForkVFluxWork {
    std::shared_ptr<MeshData> meshData;
    PredForkVFluxWork(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch A result token for Predictor Fork
struct PredForkVFluxResult {
    std::shared_ptr<MeshData> meshData;
    PredForkVFluxResult(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch B work token for Predictor Fork: WALL_BC + DIV_P1_early
struct PredForkDivWork {
    std::shared_ptr<MeshData> meshData;
    PredForkDivWork(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch B result token for Predictor Fork
struct PredForkDivResult {
    std::shared_ptr<MeshData> meshData;
    PredForkDivResult(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

#endif // PRED_FORK_DATA_H
