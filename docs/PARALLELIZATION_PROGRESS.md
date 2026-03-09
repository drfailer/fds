# FDS Hedgehog Parallelization Progress

Tracking the systematic conversion of sequential Hedgehog tasks into parallel sub-graphs using thread-safe kernels.

## Methodology

See these docs for the step-by-step procedures:

- **[METHOD_MODULE_SPLIT.md](METHOD_MODULE_SPLIT.md)** — Decomposing large Fortran modules into focused sub-modules
- **[METHOD_KERNEL_EXTRACTION.md](METHOD_KERNEL_EXTRACTION.md)** — Extracting thread-safe kernels from Fortran modules
- **[METHOD_SUBGRAPH.md](METHOD_SUBGRAPH.md)** — Converting sequential Hedgehog tasks into parallel sub-graphs

## Pipeline Overview

The full pipeline to parallelize an FDS routine:

```
1. MODULE SPLIT (if module > 5K lines)
   Large module → functional sub-modules
   See: METHOD_MODULE_SPLIT.md

2. KERNEL EXTRACTION
   Computation routine → *_kernels.f90 with TYPE(MESH_TYPE) argument
   See: METHOD_KERNEL_EXTRACTION.md

3. SUB-GRAPH CREATION
   Sequential Hedgehog task → Orchestrator → Parallel Kernel → Collector
   See: METHOD_SUBGRAPH.md
```

## Completed Sub-Graphs (7 sub-graphs, 9 graph node replacements)

### 1. Velocity Corrector (corrector phase) — Pattern A
- **Task replaced**: CorrVelocityTask
- **Kernels**: VELOCITY_CORRECTOR_KERNEL, CHECK_DIVERGENCE_KERNEL
- **Files**: data/velocity_corrector_data.h, state/velocity_corrector_state.h, task/velocity_corrector_kernel_task.h

### 2. Velocity Predictor (predictor phase) — Pattern A
- **Task replaced**: VelPredictorTask
- **Kernels**: VELOCITY_PREDICTOR_KERNEL, CHECK_STABILITY_KERNEL
- **Files**: data/velocity_predictor_data.h, state/velocity_predictor_state.h, task/velocity_predictor_kernel_task.h

### 3. Divergence Part 2 (predictor + corrector) — Pattern A
- **Tasks replaced**: DivPart2PredTask, CorrDivPart2Task (2 instances)
- **Kernel**: DIVERGENCE_PART_2_KERNEL

### 4. Corrector Step 1 (viscosity + mass FD + density) — Pattern A
- **Task replaced**: CorrStep1Task
- **Kernels**: COMPUTE_VISCOSITY_KERNEL, MASS_FINITE_DIFFERENCES_NEW_KERNEL, DENSITY_KERNEL

### 5. Density Predictor — Pattern A
- **Task replaced**: DensityPredTask
- **Kernel**: DENSITY_KERNEL

### 6. Corrector Divergence Part 1 — Pattern B
- **Task replaced**: CorrDivPart1Task
- **Sequential pre-processing**: COMBUSTION_BC (reads OMESH%Q)
- **Kernel**: DIVERGENCE_PART_1_KERNEL

### 7. Predictor/Corrector Div Setup (velocity flux) — Pattern B
- **Tasks replaced**: PredDivSetupTask, CorrDivSetupTask (2 instances)
- **Sequential pre-processing**: VISCOSITY_BC (reads OMESH%MU/D/DS) + AGGLOMERATION (corr only)
- **Kernel**: VELOCITY_FLUX_KERNEL

All verified byte-identical (DEVC) across 5 test cases (1/3/4/5-mesh configurations).

## Remaining Sequential Tasks (not parallelizable with current kernels)

| Task | Routines | Blocker |
|------|----------|---------|
| PredStep1Task | INSERT_ALL_PARTICLES, COMPUTE_VISCOSITY, MASS_FINITE_DIFFERENCES, DENSITY | INSERT_ALL_PARTICLES has no kernel |
| PredWallDivTask | WALL_BC, PARTICLE_MOMENTUM, DIVERGENCE_PART_1 | WALL_BC: 239-line orchestration, OMESH access |
| PredFinalTask | MATCH_VELOCITY, VELOCITY_BC, CC_END_STEP | MATCH_VELOCITY/VELOCITY_BC use OMESH |
| CorrCondensTask | CONDENSATION_EVAPORATION | No kernel, uses POINT_TO_MESH |
| CorrParticleTask | MOVE_PARTICLES, PARTICLE_MASS_ENERGY, REMOVE_PARTICLES | No kernels, cross-mesh particle transfer |
| CorrWallBCTask | WALL_BC | 239-line orchestration, global state, OMESH |
| CorrRadiationTask | COMPUTE_RADIATION | Complex iterative solver, no kernel |
| CorrFinalTask | MATCH_VELOCITY, VELOCITY_BC, CC_END_STEP, outputs | OMESH access in velocity routines |
| All barrier tasks | MESH_EXCHANGE, PRESSURE_ITERATION, HVAC_CALC, etc. | Inherently global/sequential |

## Next Steps

### Phase 1: Extract More Kernels

The remaining sequential tasks contain routines that **could** become kernels but haven't been extracted yet. These are the candidates for new kernel extraction work (see METHOD_KERNEL_EXTRACTION.md):

#### 1a. PredStep1Task — partial parallelization

PredStep1Task calls 4 routines. Three already have kernels (COMPUTE_VISCOSITY, MASS_FINITE_DIFFERENCES, DENSITY) but INSERT_ALL_PARTICLES does not. Options:

- **Extract INSERT_ALL_PARTICLES kernel**: Requires analysis of particle insertion logic for thread safety. If it only accesses `MESHES(NM)`, it's a candidate.
- **Split the task**: Move the 3 kernel-ready routines into a sub-graph, keep INSERT_ALL_PARTICLES sequential before it.

#### 1b. WALL_BC decomposition

WALL_BC (wall.f90) is a 239-line orchestration routine that calls multiple sub-routines. Some of its callees (PYROLYSIS, CALCULATE_RHO_D_F) already have kernel extractions in `wall_kernels.f90`. A deeper decomposition could:

- Extract the per-wall-cell computation loops as kernels
- Keep the ghost-cell exchange (OMESH access) sequential
- Apply Pattern B (sequential BC + parallel kernel)

#### 1c. CONDENSATION_EVAPORATION kernel

CONDENSATION_EVAPORATION (fire.f90) uses `POINT_TO_MESH(NM)` but does not access `OMESH`. It could be a kernel extraction candidate if the mesh-pointer variables are replaced with `M%` prefixes.

#### 1d. Particle routines

MOVE_PARTICLES and PARTICLE_MASS_ENERGY_TRANSFER use `POINT_TO_MESH` but the cross-mesh particle transfer (REMOVE_PARTICLES) makes the overall particle step hard to parallelize. Individual routines could still be kernelized.

#### 1e. COMPUTE_RADIATION

The radiation solver is iterative and complex. It may benefit from a dedicated analysis to identify parallelizable inner loops. Low priority.

### Phase 2: CC_IBM Integration

The current DivSetup sub-graph bypasses CC_IBM pre/post kernel processing (calls `fds_velocity_flux_kernel` directly instead of the full `fds_velocity_flux`). For CC_IBM cases:

- `CC_VELOCITY_BC` reads OMESH (115 occurrences) → must stay sequential
- `CUTFACE_VELOCITIES` and `CC_VELOCITY_FLUX` are in `ccib_velocity_kernels.f90` → could run in parallel
- The orchestrator would need to call CC_VELOCITY_BC sequentially, then dispatch both the main kernel and CC_IBM kernels

### Phase 3: Performance Profiling

With all possible sub-graphs implemented, profile the application to identify:

1. **Amdahl's law bottleneck**: What fraction of total time is in sequential tasks vs parallel kernels?
2. **Kernel dominance**: Which kernels consume the most wall-clock time?
3. **Scaling**: How does wall-clock time scale with `kernelThreads = 1, 2, 4, N`?
4. **Overhead**: Is the orchestrator/collector overhead significant for lightweight kernels?

### Phase 4: Multi-Process + Multi-Thread

Currently FDS uses MPI for multi-mesh (one process per mesh group). The Hedgehog integration adds intra-process parallelism (multiple threads for meshes within one process). The next frontier:

1. **Hybrid MPI+threads**: Each MPI rank runs a Hedgehog graph with `kernelThreads > 1`
2. **Load balancing**: Assign meshes to MPI ranks considering both mesh count and kernel thread availability
3. **MESH_EXCHANGE optimization**: Overlap MPI communication with kernel computation using Hedgehog's asynchronous task execution

## Unused Thread-Safe Kernels

These kernels exist in `*_kernels.f90` but are not directly used in sub-graphs (they're called internally by other kernels or from orchestration code):

| Kernel | File | Called By |
|--------|------|-----------|
| BAROCLINIC_CORRECTION_KERNEL | velo_kernels.f90 | Internal to VELOCITY_*_KERNEL |
| Turb kernels (WALE_VISCOSITY, WALL_MODEL, etc.) | turb_kernels.f90 | Called from COMPUTE_VISCOSITY_KERNEL |
| Wall kernels (PYROLYSIS, CALCULATE_RHO_D_F, etc.) | wall_kernels.f90 | Called from WALL_BC orchestration |
| CCIB kernels (21 routines) | ccib_*_kernels.f90 | Called from CC_IBM orchestration paths |
| Fire kernels (COMBUSTION_MODEL, etc.) | fire_kernels.f90 | Called from COMBUSTION_LOAD_BALANCED |
| Pressure kernels (FFT, RHS, residuals) | pres_kernels.f90 | Called from PRESSURE_ITERATION |
