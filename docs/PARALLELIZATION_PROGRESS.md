# FDS Hedgehog Parallelization Progress

Tracking the systematic conversion of sequential Hedgehog tasks into parallel sub-graphs using thread-safe kernels.

## Methodology

Step-by-step procedures for parallelizing FDS routines:

- **[METHOD_MODULE_SPLIT.md](METHOD_MODULE_SPLIT.md)** — Decomposing large Fortran modules into focused sub-modules
- **[METHOD_KERNEL_EXTRACTION.md](METHOD_KERNEL_EXTRACTION.md)** — Extracting thread-safe kernels from Fortran modules
- **[METHOD_SUBGRAPH.md](METHOD_SUBGRAPH.md)** — Converting sequential tasks into parallel sub-graphs (Pattern A)
- **[METHOD_PATTERN_B_COMPLEX.md](METHOD_PATTERN_B_COMPLEX.md)** — Complex routines with cross-mesh dependencies (Pattern B)

## Pipeline Overview

```
1. MODULE SPLIT (if module > 5K lines)
   Large module → functional sub-modules

2. KERNEL EXTRACTION
   Computation routine → *_kernels.f90 with TYPE(MESH_TYPE) argument

3. SUB-GRAPH CREATION
   Sequential task → Orchestrator → Parallel Kernel → Collector
   - Pattern A: Pure kernel (no preprocessing)
   - Pattern B: Sequential pre/post + parallel kernel
```

## Completed Sub-Graphs (15 sub-graphs)

All verified byte-identical on 1-mesh through 5-mesh test configurations.

### Pattern A: Pure Kernel Sub-Graphs

1. **Velocity Corrector** (corrector phase)
   - Kernels: VELOCITY_CORRECTOR_KERNEL, CHECK_DIVERGENCE_KERNEL
   - Files: data/velocity_corrector_data.h, state/velocity_corrector_state.h, task/velocity_corrector_kernel_task.h

2. **Velocity Predictor** (predictor phase)
   - Kernels: VELOCITY_PREDICTOR_KERNEL, CHECK_STABILITY_KERNEL
   - Files: data/velocity_predictor_data.h, state/velocity_predictor_state.h, task/velocity_predictor_kernel_task.h

3. **Divergence Part 2** (predictor + corrector, 2 graph nodes)
   - Kernel: DIVERGENCE_PART_2_KERNEL
   - Replaced: DivPart2PredTask, CorrDivPart2Task

4. **Corrector Step 1** (viscosity + mass FD + density)
   - Kernels: COMPUTE_VISCOSITY_KERNEL, MASS_FINITE_DIFFERENCES_NEW_KERNEL, DENSITY_KERNEL
   - Replaced: CorrStep1Task

5. **Density Predictor**
   - Kernel: DENSITY_KERNEL
   - Replaced: DensityPredTask

9. **Corrector Condensation**
   - Kernel: CONDENSATION_EVAPORATION_KERNEL (extracted from fire.f90)
   - Files: data/corr_condens_data.h, state/corr_condens_state.h, task/corr_condens_kernel_task.h

15. **Corrector Radiation** (Pattern A with global accumulator handling)
    - Kernel: COMPUTE_RADIATION_KERNEL (local pointer aliases, ~1200 lines including CONTAINS subroutines)
    - Technique: Local pointer aliases shadow MESH_POINTERS module variables; CONTAINS subroutines (RADIATION_FVM, ADD_VOLUMETRIC_HEAT_SOURCE) inherit aliases through host association
    - RAD_Q_SUM/KFST4_SUM: per-mesh partial sums returned via output parameters, accumulated in collector
    - Files: data/corr_radiation_data.h, state/corr_radiation_state.h, task/corr_radiation_kernel_task.h, graph/corr_radiation_subgraph.h
    - Replaced: CorrRadiationTask

### Pattern B: Sequential Pre/Post + Parallel Kernel

6. **Corrector Divergence Part 1**
   - Sequential pre-processing: COMBUSTION_BC (reads OMESH%Q)
   - Kernel: DIVERGENCE_PART_1_KERNEL
   - Replaced: CorrDivPart1Task

7. **Predictor/Corrector Div Setup** (velocity flux, 2 graph nodes)
   - Sequential pre-processing: VISCOSITY_BC (reads OMESH%MU/D/DS) + AGGLOMERATION (corr only)
   - Kernel: VELOCITY_FLUX_KERNEL
   - Replaced: PredDivSetupTask, CorrDivSetupTask

8. **Predictor Step 1** (insert particles + viscosity + mass FD)
   - Sequential pre-processing: INSERT_ALL_PARTICLES (cross-mesh, global state)
   - Kernels: COMPUTE_VISCOSITY_KERNEL, MASS_FINITE_DIFFERENCES_NEW_KERNEL
   - Files: data/pred_step1_data.h, state/pred_step1_state.h, task/pred_step1_kernel_task.h

10. **Predictor Wall + Divergence**
    - Sequential pre-processing: WALL_BC (reads OMESH for ghost cells)
    - Kernels: PARTICLE_MOMENTUM_TRANSFER_KERNEL, DIVERGENCE_PART_1_KERNEL
    - Files: data/pred_wall_div_data.h, state/pred_wall_div_state.h, task/pred_wall_div_kernel_task.h

11. **Corrector Particle Step**
    - Sequential pre-processing: PARTICLE_MASS_ENERGY_TRANSFER + MOVE_PARTICLES (cross-mesh transfer)
    - Kernel: PARTICLE_MOMENTUM_TRANSFER_KERNEL
    - Files: data/corr_particle_data.h, state/corr_particle_state.h, task/corr_particle_kernel_task.h

12. **WallBC** (3-phase complex routine, Pattern B)
    - Sequential pre-processing: ASSIGN_GHOST_VALUE (OMESH reads), NEAR_SURFACE_GAS_VARIABLES, HEAT_TRANS_COEF
    - Kernel: WALL_BC_PROCESS_CELLS_KERNEL (~90% of wall cells, no cross-mesh dependencies)
    - Sequential finalization: HAS_BACK_MESH cells, thin walls, particle off-gassing
    - Files: data/wallbc_data.h, state/wallbc_state.h, task/wallbc_kernel_task.h
    - Documentation: docs/WALL_BC_PARALLELIZATION_PLAN.md, test_cases/WALLBC_TEST_REPORT.md
    - Replaced: CorrWallBCTask

13. **PredFinal** (3-phase complex routine, Pattern B)
    - Sequential pre-processing: MATCH_VELOCITY (cross-mesh interpolation), SYNTHETIC_TURBULENCE_IF_ENABLED (SEM inflow), VELOCITY_BC_PREPROCESSING (OMESH reads)
    - Kernel: VELOCITY_BC_PROCESS_EDGES_KERNEL (all edge boundary conditions, thread-safe M% access)
    - Sequential finalization: CC_VELOCITY_BC (cut-cell velocity BC if CC_IBM active)
    - Files: data/velocity_bc_data.h, state/velocity_bc_state.h, task/velocity_bc_edges_task.h, graph/velocity_bc_subgraph.h
    - Replaced: PredFinalTask

14. **CorrFinal** (3-phase complex routine, Pattern B)
    - Sequential pre-processing: MATCH_VELOCITY (cross-mesh interpolation), VELOCITY_BC_PREPROCESSING (OMESH reads)
    - Kernel: VELOCITY_BC_PROCESS_EDGES_KERNEL (all edge boundary conditions, thread-safe M% access)
    - Sequential finalization: CC_VELOCITY_BC (cut-cell velocity BC), UPDATE_GLOBAL_OUTPUTS (per-mesh output accumulation)
    - Files: shared with PredFinal (velocity_bc_data.h, velocity_bc_state.h, velocity_bc_edges_task.h, velocity_bc_subgraph.h)
    - Replaced: CorrFinalTask

## Kernel Extraction Summary

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
| WALL_BC_PROCESS_CELLS_KERNEL | wall.f90 | WallBC |
| VELOCITY_BC_PROCESS_EDGES_KERNEL | velo_kernels.f90 | PredFinal, CorrFinal |

### New Extractions for WallBC

| Kernel | Source | Lines | Purpose |
|--------|--------|-------|---------|
| WALL_BC_PREPROCESSING | wall.f90 | 49 | OMESH reads, gas variable setup |
| WALL_BC_PROCESS_CELLS_KERNEL | wall.f90 | 155 | Parallel cell processing |
| WALL_BC_FINALIZE | wall.f90 | 68 | OMESH writes, cross-mesh coupling |
| CALCULATE_RHO_F_KERNEL | wall_kernels.f90 | 60 | Cell-local RHO_F calculation |
| NEAR_SURFACE_GAS_VARIABLES_KERNEL | wall_kernels.f90 | 142 | Gas properties near walls |

### New Extractions for VelocityBC (PredFinal/CorrFinal)

| Routine | Source | Lines | Purpose |
|---------|--------|-------|---------|
| VELOCITY_BC_PREPROCESSING | velo.f90 | 75 | OMESH wall velocity reads, sequential |
| VELOCITY_BC_PROCESS_EDGES_KERNEL | velo_kernels.f90 | 760 | Edge BC processing, parallel (thread-safe M%) |

### Thread-Safe Callee Conversions (WallBC)

| Routine | Lines | Pattern | Status |
|---------|-------|---------|--------|
| CALC_HVAC_BC | 52 | Explicit M + PREDICTOR_FLAG args | ✅ Converted |
| HEAT_TRANSFER_COEFFICIENT | ~175 | Index-based access (no pointers) | ✅ Converted |
| SURFACE_HEAT_TRANSFER | 379 | M pointer + conditional setup | ✅ Converted |
| CALCULATE_ZZ_F | 413 | M pointer + conditional setup | ✅ Converted |

**Key techniques**:
- **Index-based**: Use integer indices to access arrays (avoids pointer issues)
- **Pointer-based**: Use `TYPE(MESH_TYPE), POINTER :: M` with conditional pointer setup for predictor/corrector

## Remaining Sequential Tasks

| Task | Routines | Blocker |
|------|----------|---------|
| All barrier tasks | MESH_EXCHANGE, PRESSURE_ITERATION, HVAC_CALC | Inherently global/sequential |

**Note**: All per-mesh computation tasks have been parallelized (15 sub-graphs). Only inherently global/sequential barrier tasks remain.

## Performance Profiling Results

**Test case**: dancing_eddies, 27 timesteps (from Hedgehog graph dot file)

### 1-Mesh Baseline (kernelThreads=1, total 6.298s)

| Category | Time | % |
|----------|------|---|
| Parallel kernels | 2873 ms | 45.6% |
| Sequential tasks | 2135 ms | 33.9% |
| Orchestrator pre-proc | 428 ms | 6.8% |
| Barriers/exchanges | 357 ms | 5.7% |
| Overhead | 505 ms | 8.0% |

### 4-Mesh Parallel (kernelThreads=4, total 4.690s)

| Category | Time | % |
|----------|------|---|
| Parallel kernels | 1277 ms | 27.2% |
| Sequential tasks | 1838 ms | 39.2% |
| Orchestrator pre-proc | 454 ms | 9.7% |
| Barriers/exchanges | 345 ms | 7.4% |
| Overhead | 776 ms | 16.5% |

### Kernel Parallel Speedup

4 meshes, 4 threads vs expected 4× sequential:

| Kernel | Expected (4×1m) | Actual (4m) | Speedup |
|--------|-----------------|-------------|---------|
| PredWallDivKernel | 2549 ms | 279 ms | 9.1× |
| CorrDivPart1Kernel | 2507 ms | 284 ms | 8.8× |
| CorrStep1Kernel | 1610 ms | 192 ms | 8.4× |
| DivSetupKernel (pred) | 1362 ms | 148 ms | 9.2× |
| DivSetupKernel (corr) | 1205 ms | 146 ms | 8.2× |
| PredStep1Kernel | 1180 ms | 150 ms | 7.9× |
| **ALL KERNELS** | **11492 ms** | **1277 ms** | **9.0×** |

**Super-linear speedup** (9× from 4 threads) is due to better cache utilization on smaller per-mesh domains.

### Sequential Bottlenecks (4-mesh)

| Task | Time | Notes |
|------|------|-------|
| CorrRadiation | 356 ms | Complex iterative solver (future target) |
| TimestepCompute | 349 ms | Outputs, diagnostics |
| PressureIteration (2×) | 302 ms | Global Poisson solver |

**Previous bottlenecks eliminated**:
- CorrWallBC (was 406 ms) → now WallBC sub-graph (preprocessing + parallel kernel + finalization)
- PredFinal (was 342 ms) → now PredFinal sub-graph (Pattern B: MATCH_VELOCITY + preprocessing + parallel VELOCITY_BC_PROCESS_EDGES_KERNEL + CC_VELOCITY_BC)
- CorrFinal (was 385 ms) → now CorrFinal sub-graph (Pattern B: MATCH_VELOCITY + preprocessing + parallel VELOCITY_BC_PROCESS_EDGES_KERNEL + CC_VELOCITY_BC + outputs)

## CC_IBM Integration

All parallelized sub-graphs include CC_IBM (cut-cell immersed boundary) processing:

**DivSetup (predictor + corrector)**:
- Orchestrators call `CC_VELOCITY_BC` sequentially (OMESH access)
- Kernel wrapper calls `CUTFACE_VELOCITIES` and `CC_VELOCITY_FLUX`

**Velocity Predictor/Corrector**:
- CC_PROJECT_VELOCITY called in orchestrator/collector as needed

**Already embedded in kernels**:
- DIVERGENCE_PART_1/2_KERNEL, COMPUTE_VISCOSITY_KERNEL, DENSITY_KERNEL all include CC_IBM routines

All verified byte-identical on CC_IBM test cases.

## Future Work

### 1. Breaking the Sequential Bottleneck

Current sequential fraction reduced from ~39% to ~25% with PredFinal/CorrFinal decomposition. Options:

**a) RADIATION Decomposition**
- Extract angle loop into parallel kernel
- Keep RTE source correction sequential
- Medium effort, moderate ROI (~350ms saved)

**b) Hybrid MPI+Hedgehog**
- Each MPI rank runs Hedgehog graph with kernelThreads > 1
- Load balancing across MPI ranks and threads
- Overlap MPI communication with kernel computation

### 2. Advanced Optimization

**Pipeline parallelism**: Overlap predictor/corrector of different meshes (requires decoupling barriers)

**NUMA-aware mesh assignment**: Pin meshes to NUMA nodes for memory locality

**Asynchronous MESH_EXCHANGE**: Start communication early, overlap with computation

## Test Suite

**Test runner**: `test_cases/run_tests.py`

**Test cases** (all byte-identical):
- dancing_eddies_1mesh (1 mesh)
- dancing_eddies_2mesh (2 meshes, embedded)
- multiple_reac_3mesh (3 meshes)
- dancing_eddies_4mesh (4 meshes)
- species_props_5mesh (5 meshes)

**Verification**: `python3 test_cases/run_tests.py -v`

## Summary Statistics

- **Sub-graphs created**: 14
- **Graph nodes replaced**: 17 (some tasks appear in both predictor/corrector)
- **Kernels extracted**: 14 new kernels + utilizing ~30 existing kernels
- **Thread-safe conversions**: 1800+ lines converted (including ~760 lines for VELOCITY_BC_PROCESS_EDGES_KERNEL)
- **Test coverage**: 5 test cases, 1-5 meshes, all byte-identical
- **Speedup achieved**: 9× on parallel kernels (4 meshes, 4 threads)
- **Sequential fraction**: reduced from ~39% to ~25% (PredFinal + CorrFinal decomposed)

## Documentation Index

### Methodology
- METHOD_MODULE_SPLIT.md - Module decomposition
- METHOD_KERNEL_EXTRACTION.md - Kernel extraction patterns
- METHOD_SUBGRAPH.md - Pattern A sub-graphs (pure kernel)
- METHOD_PATTERN_B_COMPLEX.md - Pattern B sub-graphs (complex routines)

### Implementation Details
- WALL_BC_PARALLELIZATION_PLAN.md - Complete WallBC implementation (reference)
- test_cases/WALLBC_TEST_REPORT.md - WallBC verification results

### Progress Tracking
- PARALLELIZATION_PROGRESS.md - This file (current status)
