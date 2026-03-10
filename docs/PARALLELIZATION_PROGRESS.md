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

## Completed Sub-Graphs (11 sub-graphs, 13 graph node replacements)

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

### 8. Predictor Step 1 (insert particles + viscosity + mass FD) — Pattern B
- **Task replaced**: PredStep1Task
- **Sequential pre-processing**: INSERT_ALL_PARTICLES (cross-mesh, global state)
- **Kernels**: COMPUTE_VISCOSITY_KERNEL, MASS_FINITE_DIFFERENCES_NEW_KERNEL
- **Files**: data/pred_step1_data.h, state/pred_step1_state.h, task/pred_step1_kernel_task.h

### 9. Corrector Condensation — Pattern A
- **Task replaced**: CorrCondensTask
- **Kernel**: CONDENSATION_EVAPORATION_KERNEL (new extraction from fire.f90)
- **Files**: data/corr_condens_data.h, state/corr_condens_state.h, task/corr_condens_kernel_task.h

### 10. Predictor Wall + Divergence — Pattern B
- **Task replaced**: PredWallDivTask
- **Sequential pre-processing**: WALL_BC (reads OMESH for ghost cells)
- **Kernels**: PARTICLE_MOMENTUM_TRANSFER_KERNEL (new extraction from part.f90), DIVERGENCE_PART_1_KERNEL
- **Files**: data/pred_wall_div_data.h, state/pred_wall_div_state.h, task/pred_wall_div_kernel_task.h

### 11. Corrector Particle Step — Pattern B
- **Task replaced**: CorrParticleTask
- **Sequential pre-processing**: PARTICLE_MASS_ENERGY_TRANSFER + MOVE_PARTICLES (cross-mesh transfer)
- **Kernel**: PARTICLE_MOMENTUM_TRANSFER_KERNEL
- **Files**: data/corr_particle_data.h, state/corr_particle_state.h, task/corr_particle_kernel_task.h

All verified byte-identical (DEVC) across 1-mesh and 4-mesh test configurations.

## New Kernel Extractions (this session)

| Kernel | Source Module | Kernel Module | Lines |
|--------|---------------|---------------|-------|
| PARTICLE_MOMENTUM_TRANSFER_KERNEL | part.f90 | part_kernels.f90 | 48 |
| CONDENSATION_EVAPORATION_KERNEL | fire.f90 | fire_kernels.f90 | 303 |
| CALCULATE_RHO_F_KERNEL | wall.f90 | wall_kernels.f90 | 60 |
| NEAR_SURFACE_GAS_VARIABLES_KERNEL | wall.f90 | wall_kernels.f90 | 142 |
| SCALAR_TO_POINT_K | wall.f90 | wall_kernels.f90 | 20 |
| GET_TRILINEAR_WEIGHTS_K | wall.f90 | wall_kernels.f90 | 55 |

## Remaining Sequential Tasks (not parallelizable)

| Task | Routines | Blocker |
|------|----------|---------|
| PredFinalTask | MATCH_VELOCITY, VELOCITY_BC, CC_END_STEP | MATCH_VELOCITY/VELOCITY_BC use OMESH |
| CorrWallBCTask | WALL_BC | 239-line orchestration, OMESH in ASSIGN_GHOST_VALUE, SURFACE_HEAT_TRANSFER |
| CorrRadiationTask | COMPUTE_RADIATION | Complex iterative solver, no kernel, low priority |
| CorrFinalTask | MATCH_VELOCITY, VELOCITY_BC, CC_END_STEP, outputs | OMESH access in velocity routines |
| All barrier tasks | MESH_EXCHANGE, PRESSURE_ITERATION, HVAC_CALC, etc. | Inherently global/sequential |

## Analysis of Remaining Tasks

### PredFinalTask / CorrFinalTask
MATCH_VELOCITY and VELOCITY_BC both read OMESH data. These routines coordinate velocity values across mesh boundaries — fundamentally sequential. No kernel extraction possible without redesigning the inter-mesh velocity matching.

### CorrWallBCTask
WALL_BC (wall.f90) is a 239-line orchestration routine. Key sub-routines:
- **ASSIGN_GHOST_VALUE**: Heavy OMESH access (ghost cell interpolation) → must stay sequential
- **SURFACE_HEAT_TRANSFER**: 95% cell-local, but INTERPOLATED_BC case uses OMESH → mixed
- **CALCULATE_ZZ_F**: 95% cell-local, but CONSUME_MASS uses OMESH → mixed
- **CALCULATE_RHO_F**: Pure cell-local → ✓ extracted as CALCULATE_RHO_F_KERNEL
- **NEAR_SURFACE_GAS_VARIABLES**: Pure cell-local → ✓ extracted as NEAR_SURFACE_GAS_VARIABLES_KERNEL (with helpers SCALAR_TO_POINT_K, GET_TRILINEAR_WEIGHTS_K)
- Already extracted kernels: CALCULATE_RHO_D_F, CALC_DEPOSITION, PYROLYSIS (in wall_kernels.f90)

Decomposition is possible but requires splitting individual sub-routines (e.g., SURFACE_HEAT_TRANSFER) into OMESH and non-OMESH parts. Medium complexity, moderate ROI.

### CorrRadiationTask
COMPUTE_RADIATION is an iterative solver with complex internal state management. Low priority for parallelization.

### Barrier Tasks
MESH_EXCHANGE, PRESSURE_ITERATION, HVAC_CALC are inherently global synchronization points that operate across all meshes simultaneously. They cannot be parallelized within the current architecture.

## Next Steps

### Phase 1: WALL_BC Decomposition — Partially Complete

**Completed kernel extractions:**
- ✓ CALCULATE_RHO_F → CALCULATE_RHO_F_KERNEL (60 lines, wall_kernels.f90)
- ✓ NEAR_SURFACE_GAS_VARIABLES → NEAR_SURFACE_GAS_VARIABLES_KERNEL (142 lines, wall_kernels.f90)
  - Helpers: SCALAR_TO_POINT_K (20 lines), GET_TRILINEAR_WEIGHTS_K (55 lines)
- ✓ Old routines removed from wall.f90, all call sites updated

**Deferred — low ROI:**
- SURFACE_HEAT_TRANSFER (378 lines): INTERPOLATED_BC case (lines 616-775) deeply interleaves OMESH data with local computation. Splitting would require conditional dispatch at call sites. Large effort, small parallelizable fraction.
- CALCULATE_ZZ_F (413 lines): CONSUME_MASS (24 lines) is the only OMESH-dependent part, but the remaining 389 lines use many module-level pointer aliases (WALL, BOUNDARY_PROP1/2, CELL_INDEX, etc.) requiring M% conversion. Large effort, marginal gain.
- SOLID_HEAT_TRANSFER (1177 lines): Uses BACK_MESH (cross-mesh) for back-to-back wall cells. Not a kernel candidate.

**Conclusion:** A full WALL_BC sub-graph is not feasible without redesigning the ASSIGN_GHOST_VALUE call chain. The extracted kernels (CALCULATE_RHO_F_KERNEL, NEAR_SURFACE_GAS_VARIABLES_KERNEL) are available for future use if the WALL_BC orchestration is refactored.

### Phase 2: CC_IBM Integration ✓ COMPLETED

All parallelized sub-graphs now include CC_IBM processing:

**DivSetup (predictor + corrector):**
- Orchestrators call `CC_VELOCITY_BC` sequentially (OMESH access)
- Kernel wrapper calls `CUTFACE_VELOCITIES` (pre/post) and `CC_VELOCITY_FLUX` from `ccib_velocity_kernels.f90`

**Velocity Predictor:**
- Collector calls `CC_PROJECT_VELOCITY(STORE=.FALSE.)` after kernel execution

**Velocity Corrector:**
- Orchestrator calls `CC_PROJECT_VELOCITY(STORE=.TRUE.)` before kernel dispatch
- Collector calls `CC_PROJECT_VELOCITY(STORE=.FALSE.)` after kernel execution

**Already embedded in kernels (no changes needed):**
- `DIVERGENCE_PART_1_KERNEL`: includes `CC_DIVERGENCE_PART_1`, `SET_EXIMDIFFLX_3D`, etc.
- `DIVERGENCE_PART_2_KERNEL`: includes `GET_CUTCELL_DDDT`, solid cell zeroing
- `COMPUTE_VISCOSITY_KERNEL`: includes `CUTFACE_VELOCITIES`, `CC_COMPUTE_KRES`, `CC_COMPUTE_VISCOSITY`
- `DENSITY_KERNEL`: includes `SET_EXIMADVFLX_3D`
- `PARTICLE_MOMENTUM_TRANSFER_KERNEL`: includes `CUTFACE_VELOCITIES`

All verified byte-identical on 1-mesh and 4-mesh tests (non-CC_IBM cases).

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

## Complete Kernel Module Inventory

### Directly Used in Sub-Graphs

| Kernel | File | Sub-Graph |
|--------|------|-----------|
| VELOCITY_PREDICTOR_KERNEL | velo_kernels.f90 | VelocityPredictor |
| CHECK_STABILITY_KERNEL | velo_kernels.f90 | VelocityPredictor |
| VELOCITY_CORRECTOR_KERNEL | velo_kernels.f90 | VelocityCorrector |
| CHECK_DIVERGENCE_KERNEL | divg_kernels.f90 | VelocityCorrector |
| DIVERGENCE_PART_1_KERNEL | divg_kernels.f90 | CorrDivPart1, PredWallDiv |
| DIVERGENCE_PART_2_KERNEL | divg_kernels.f90 | DivPart2 (pred+corr) |
| COMPUTE_VISCOSITY_KERNEL | velo_kernels.f90 | CorrStep1, PredStep1 |
| VELOCITY_FLUX_KERNEL | velo_kernels.f90 | DivSetup (pred+corr) |
| MASS_FINITE_DIFFERENCES_NEW_KERNEL | mass_kernels.f90 | CorrStep1, PredStep1 |
| DENSITY_KERNEL | mass_kernels.f90 | CorrStep1, DensityPred |
| CONDENSATION_EVAPORATION_KERNEL | fire_kernels.f90 | CorrCondens |
| PARTICLE_MOMENTUM_TRANSFER_KERNEL | part_kernels.f90 | PredWallDiv, CorrParticle |

### Indirectly Used (called by other kernels)

| Kernel | File | Called By |
|--------|------|-----------|
| BAROCLINIC_CORRECTION_KERNEL | velo_kernels.f90 | Internal to VELOCITY_*_KERNEL |
| Turb kernels (WALE_VISCOSITY, WALL_MODEL, etc.) | turb_kernels.f90 | Called from COMPUTE_VISCOSITY_KERNEL |
| Wall kernels (PYROLYSIS, CALCULATE_RHO_D_F, CALCULATE_RHO_F_KERNEL, NEAR_SURFACE_GAS_VARIABLES_KERNEL, etc.) | wall_kernels.f90 | Called from WALL_BC orchestration |
| CCIB kernels (21 routines) | ccib_*_kernels.f90 | Called from CC_IBM orchestration paths |
| Fire kernels (COMBUSTION_MODEL, etc.) | fire_kernels.f90 | Called from COMBUSTION_LOAD_BALANCED |
| Pressure kernels (FFT, RHS, residuals) | pres_kernels.f90 | Called from PRESSURE_ITERATION |
