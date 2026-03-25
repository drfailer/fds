#ifndef BAROCLINIC_KERNEL_TASK_H
#define BAROCLINIC_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel baroclinic correction kernel task.
///
/// Calls fds_baroclinic_correction(t, nm) per mesh. The Fortran routine
/// internally checks BAROCLINIC/SOLID_PHASE_ONLY/FREEZE_VELOCITY and is
/// a no-op when correction is not applicable.
///
/// Multi-threaded: each clone processes one mesh independently.
class BaroclinicKernelTask
    : public hh::AbstractTask<2, PressureIterMeshData, MeshData, MeshData> {
public:
    explicit BaroclinicKernelTask(size_t kernelThreads)
        : hh::AbstractTask<2, PressureIterMeshData, MeshData, MeshData>(
              "BaroclinicKernel", kernelThreads) {}

    // simply unwrap the mesh when starting a new iteration
    void execute(std::shared_ptr<PressureIterMeshData> pimd) override {
        execute(pimd->mesh);
    }

    void execute(std::shared_ptr<MeshData> md) override {
        // Only apply baroclinic correction when ITERATE_BAROCLINIC_TERM is true.
        // This flag starts as true (set by init) and is cleared by the
        // convergence check when pressure error drops below tolerance.
        // The flag is stable during parallel execution (set/cleared only
        // in barrier contexts between iterations).
        if (fds_pressure_iteration_needs_baroclinic()) {
            fds_baroclinic_correction(md->t, md->nm);
        }
        this->addResult(md);
    }

    std::shared_ptr<hh::AbstractTask<2, PressureIterMeshData, MeshData, MeshData>> copy() override {
        return std::make_shared<BaroclinicKernelTask>(this->numberThreads());
    }
};

#endif // BAROCLINIC_KERNEL_TASK_H
