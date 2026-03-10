#ifndef PRED_STEP1_DATA_H
#define PRED_STEP1_DATA_H

#include <memory>
#include "mesh_data.h"

/// Work token for the predictor step 1 sub-graph.
/// Carries parameters needed by COMPUTE_VISCOSITY_KERNEL + MASS_FINITE_DIFFERENCES_NEW_KERNEL.
struct PredStep1Work {
    int nm;
    double t;
    double dt;
    std::shared_ptr<MeshData> originalMeshData;

    PredStep1Work(int nm_, double t_, double dt_, std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // PRED_STEP1_DATA_H
