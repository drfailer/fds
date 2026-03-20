#ifndef PIPELINE_FORK2_DATA_H
#define PIPELINE_FORK2_DATA_H

#include <memory>
#include <vector>
#include "mesh_data.h"
#include "barrier_data.h"

/// Branch C work token for Fork 2: RADIATION (per mesh)
struct Fork2RadWork {
    std::shared_ptr<MeshData> meshData;
    Fork2RadWork(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch C barrier token: all N meshes completed radiation
struct Fork2RadBarrier {
    std::shared_ptr<BarrierData> barrierData;
    Fork2RadBarrier(std::shared_ptr<BarrierData> bd) : barrierData(std::move(bd)) {}
};

/// Branch D work token for Fork 2: DIV_P1 (per mesh)
struct Fork2DivP1Work {
    std::shared_ptr<MeshData> meshData;
    Fork2DivP1Work(std::shared_ptr<MeshData> md) : meshData(std::move(md)) {}
};

/// Branch D barrier token: all N meshes completed DIV_P1
struct Fork2DivP1Barrier {
    std::shared_ptr<BarrierData> barrierData;
    Fork2DivP1Barrier(std::shared_ptr<BarrierData> bd) : barrierData(std::move(bd)) {}
};

#endif // PIPELINE_FORK2_DATA_H
