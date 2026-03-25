#ifndef PRESSURE_ITERATION_DATA_H
#define PRESSURE_ITERATION_DATA_H

#include <memory>
#include <vector>
#include "mesh_data.h"

/// MeshData wrapper for the pressure iteration cycle.
struct PressureIterMeshData {
    std::shared_ptr<MeshData> mesh;
};

#endif // PRESSURE_ITERATION_DATA_H
