#ifndef PIPELINE_FORK1_DATA_H
#define PIPELINE_FORK1_DATA_H

#include <memory>
#include "mesh_data.h"

/// Branch A work token for Fork 1: VELOCITY_FLUX
struct Fork1VFluxWork {
    std::shared_ptr<MeshData> meshData;
    Fork1VFluxWork(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch A result token for Fork 1
struct Fork1VFluxResult {
    std::shared_ptr<MeshData> meshData;
    Fork1VFluxResult(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch B work token for Fork 1: COMBUSTION
struct Fork1CombWork {
    std::shared_ptr<MeshData> meshData;
    Fork1CombWork(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch B result token for Fork 1
struct Fork1CombResult {
    std::shared_ptr<MeshData> meshData;
    Fork1CombResult(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

#endif // PIPELINE_FORK1_DATA_H
