# FDS Hedgehog Parallelization Progress

Tracking the systematic conversion of sequential Hedgehog tasks into parallel sub-graphs using thread-safe kernels.

## Pattern

Each sub-graph replaces a sequential task with:
```
[Orchestrator State] collects N mesh tokens
    |
[Parallel Kernel Task] processes N meshes concurrently (kernelThreads)
    |
[Collector State] gathers results, sorts by NM, emits N tokens
```

## Completed Sub-Graphs

### 1. Velocity Corrector (corrector phase)
- **Task replaced**: CorrVelocityTask
- **Kernels**: VELOCITY_CORRECTOR_KERNEL, CHECK_DIVERGENCE_KERNEL
- **Files**: data/velocity_corrector_data.h, state/velocity_corrector_state.h, task/velocity_corrector_kernel_task.h
- **Status**: Implemented and verified byte-identical

### 2. Velocity Predictor (predictor phase)
- **Task replaced**: VelPredictorTask
- **Kernels**: VELOCITY_PREDICTOR_KERNEL, CHECK_STABILITY_KERNEL
- **Files**: data/velocity_predictor_data.h, state/velocity_predictor_state.h, task/velocity_predictor_kernel_task.h
- **Status**: Implemented and verified byte-identical

### 3. Divergence Part 2 (predictor + corrector)
- **Tasks replaced**: DivPart2PredTask, CorrDivPart2Task
- **Kernel**: DIVERGENCE_PART_2_KERNEL(M, DT, NM)
- **Files**: data/divergence_part2_data.h, state/divergence_part2_state.h, task/divergence_part2_kernel_task.h
- **Status**: COMPLETE - verified byte-identical (DEVC) across all 5 test cases

## Planned Sub-Graphs

### 4. Corrector Step 1 (viscosity + mass FD + density)
- **Task replaced**: CorrStep1Task
- **Kernels**: COMPUTE_VISCOSITY_KERNEL, MASS_FINITE_DIFFERENCES_NEW_KERNEL, DENSITY_KERNEL
- **Files**: data/corr_step1_data.h, state/corr_step1_state.h, task/corr_step1_kernel_task.h
- **Status**: COMPLETE - verified byte-identical (DEVC) across all 5 test cases

### 5. Density Predictor
- **Task replaced**: DensityPredTask
- **Kernel**: DENSITY_KERNEL(M, T, DT, NM)
- **Files**: data/density_pred_data.h, state/density_pred_state.h, task/density_pred_kernel_task.h
- **Status**: COMPLETE - verified byte-identical (DEVC) across all 5 test cases

### 6. Corrector Divergence Part 1
- **Task to replace**: CorrDivPart1Task
- **Kernel**: DIVERGENCE_PART_1_KERNEL(M, T, DT, NM)
- **Pattern**: Pre-processing — sequential COMBUSTION_BC (OMESH access) in orchestrator, parallel kernel
- **Difficulty**: Medium — COMBUSTION_BC reads OMESH(NOM)%Q for ghost cells
- **Status**: NOT STARTED

### 7. Predictor/Corrector Div Setup (velocity flux)
- **Tasks to replace**: PredDivSetupTask, CorrDivSetupTask
- **Kernel**: VELOCITY_FLUX_KERNEL(M, T, DT, NM, ...)
- **Pattern**: Pre-processing — sequential VISCOSITY_BC (OMESH access), parallel VELOCITY_FLUX_KERNEL
- **Difficulty**: Medium — VISCOSITY_BC reads OMESH; CC_IBM has pre/post kernel paths
- **Status**: NOT STARTED

## Not Parallelizable (no kernel or cross-mesh dependency)

| Task | Reason |
|------|--------|
| CorrWallBCTask | WALL_BC: 239-line orchestration, no kernel, global state |
| CorrRadiationTask | COMPUTE_RADIATION: no kernel, complex iterative solver |
| CorrCondensTask | CONDENSATION_EVAPORATION: no kernel |
| CorrParticleTask | Particle routines: no kernels, cross-mesh transfer |
| PredFinalTask / CorrFinalTask | MATCH_VELOCITY, VELOCITY_BC use OMESH |
| PredStep1Task | INSERT_ALL_PARTICLES has no kernel (2/3 routines have kernels) |
| Barrier tasks | Inherently sequential (MESH_EXCHANGE, PRESSURE_ITERATION, etc.) |

## Available Thread-Safe Kernels (not yet used in sub-graphs)

| Kernel | File | Used By |
|--------|------|---------|
| BAROCLINIC_CORRECTION_KERNEL | velo_kernels.f90 | (internal to other kernels) |
| COMPUTE_VISCOSITY_KERNEL | velo_kernels.f90 | Planned: CorrStep1 |
| VELOCITY_FLUX_KERNEL | velo_kernels.f90 | Planned: DivSetup |
| MASS_FINITE_DIFFERENCES_NEW_KERNEL | mass_kernels.f90 | Planned: CorrStep1 |
| DENSITY_KERNEL | mass_kernels.f90 | Planned: CorrStep1, DensityPred |
| DIVERGENCE_PART_1_KERNEL | divg_kernels.f90 | Planned: CorrDivPart1 |
| Turb kernels (EX2G3D, FILL_EDGES, etc.) | turb_kernels.f90 | (called from COMPUTE_VISCOSITY_KERNEL) |
| Wall kernels (PYROLYSIS, etc.) | wall_kernels.f90 | (called from WALL_BC orchestration) |
| CCIB kernels | ccib_*_kernels.f90 | (called from CC_IBM paths) |
