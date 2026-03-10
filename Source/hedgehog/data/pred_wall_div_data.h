#ifndef PRED_WALL_DIV_DATA_H
#define PRED_WALL_DIV_DATA_H

#include <memory>
#include "mesh_data.h"

/// Work token for the predictor wall+div sub-graph.
/// Carries parameters needed by PARTICLE_MOMENTUM_KERNEL + DIVERGENCE_PART_1_KERNEL.
struct PredWallDivWork {
    int nm;
    double t;
    double dt;
    std::shared_ptr<MeshData> originalMeshData;

    PredWallDivWork(int nm_, double t_, double dt_, std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // PRED_WALL_DIV_DATA_H
