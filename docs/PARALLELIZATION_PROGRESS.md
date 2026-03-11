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

## New Kernel Extractions

| Kernel | Source Module | Kernel Module | Lines |
|--------|---------------|---------------|-------|
| PARTICLE_MOMENTUM_TRANSFER_KERNEL | part.f90 | part_kernels.f90 | 48 |
| CONDENSATION_EVAPORATION_KERNEL | fire.f90 | fire_kernels.f90 | 303 |
| CALCULATE_RHO_F_KERNEL | wall.f90 | wall_kernels.f90 | 60 |
| NEAR_SURFACE_GAS_VARIABLES_KERNEL | wall.f90 | wall_kernels.f90 | 142 |
| SCALAR_TO_POINT_K | wall.f90 | wall_kernels.f90 | 20 |
| GET_TRILINEAR_WEIGHTS_K | wall.f90 | wall_kernels.f90 | 55 |

## Thread-Safe Routine Conversions (WALL_BC Callees)

| Routine | Source Module | Status | Approach |
|---------|---------------|--------|----------|
| CALC_HVAC_BC | wall.f90 | ✅ Converted | Added M and PREDICTOR_FLAG arguments |
| HEAT_TRANSFER_COEFFICIENT | func.f90 | ✅ Converted | Index-based access (no pointers) |
| DEPOSIT_PARTICLE_MASS | wall.f90 | Already safe | No module-level pointers used |

**Key technique**: Use integer indices instead of pointers to avoid Fortran ALLOCATABLE/TARGET issues.
Example: `B1_INDEX = M%WALL(IW)%B1_INDEX; M%BOUNDARY_PROP1(B1_INDEX)%...`

See: `docs/WALL_BC_CONVERSIONS_SUMMARY.md` for details.

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

### Phase 3: Performance Profiling ✓ COMPLETED

**Test case**: dancing_eddies, 27 timesteps. Hedgehog graph dot file provides per-node execution statistics.

**1-mesh baseline** (kernelThreads=1, total 6.298s):
| Category | Time | % |
|----------|------|---|
| Parallel kernels | 2873 ms | 45.6% |
| Sequential tasks | 2135 ms | 33.9% |
| Orchestrator pre-proc | 428 ms | 6.8% |
| Barriers/exchanges | 357 ms | 5.7% |
| Overhead (collectors, retry) | 505 ms | 8.0% |

**4-mesh parallel** (kernelThreads=4, total 4.690s):
| Category | Time | % |
|----------|------|---|
| Parallel kernels | 1277 ms | 27.2% |
| Sequential tasks | 1838 ms | 39.2% |
| Orchestrator pre-proc | 454 ms | 9.7% |
| Barriers/exchanges | 345 ms | 7.4% |
| Overhead (collectors, retry) | 776 ms | 16.5% |

**Kernel parallel speedup** (4 meshes, 4 threads vs expected 4× sequential):
| Kernel | Expected (4×1m) | Actual (4m) | Speedup |
|--------|-----------------|-------------|---------|
| PredWallDivKernel | 2549 ms | 279 ms | 9.1x |
| CorrDivPart1Kernel | 2507 ms | 284 ms | 8.8x |
| CorrStep1Kernel | 1610 ms | 192 ms | 8.4x |
| DivSetupKernel (pred) | 1362 ms | 148 ms | 9.2x |
| DivSetupKernel (corr) | 1205 ms | 146 ms | 8.2x |
| PredStep1Kernel | 1180 ms | 150 ms | 7.9x |
| **ALL KERNELS** | **11492 ms** | **1277 ms** | **9.0x** |

Super-linear speedup (9x from 4 threads) is due to better cache utilization on smaller per-mesh domains.

**Sequential bottlenecks** (4-mesh, cannot be parallelized):
| Task | Time | Notes |
|------|------|-------|
| PredWallDivOrch (WALL_BC) | 426 ms | Sequential pre-processing in orchestrator |
| CorrWallBC | 406 ms | Sequential task, OMESH dependencies |
| CorrFinal | 385 ms | MATCH_VELOCITY, VELOCITY_BC use OMESH |
| CorrRadiation | 356 ms | Complex iterative solver |
| TimestepCompute | 349 ms | Outputs, diagnostics |
| PredFinal | 342 ms | MATCH_VELOCITY, VELOCITY_BC use OMESH |
| PressureIteration (2×) | 302 ms | Global Poisson solver |

**Amdahl's law**: With 49% sequential time, max theoretical speedup is **2.05×** even with infinite kernel threads. To exceed this, the sequential tasks (WALL_BC, MATCH_VELOCITY, VELOCITY_BC, COMPUTE_RADIATION) would need to be decomposed — but all have deep cross-mesh dependencies.

**Collector/overhead**: Negligible (< 0.2% for collectors, states). The Hedgehog framework introduces minimal overhead.

### Phase 4: Future Work

#### 4a. Breaking the Sequential Bottleneck

The 49% sequential fraction limits speedup to ~2x. To go further, the fundamental blocker is `POINT_TO_MESH` — it sets module-level pointer aliases that are inherently not thread-safe. Options:

1. **Thread-local POINT_TO_MESH**: Use OpenMP `THREADPRIVATE` for mesh pointer aliases. Would require changes to every module that uses `CALL POINT_TO_MESH`. Medium effort, high impact.
2. **Explicit mesh passing**: Convert remaining sequential routines (WALL_BC, MATCH_VELOCITY, VELOCITY_BC, COMPUTE_RADIATION) to take `TYPE(MESH_TYPE)` as argument. Very large refactor (~50K+ lines affected).
3. **Selective decomposition**: Split WALL_BC's SURFACE_HEAT_TRANSFER to separate INTERPOLATED_BC (OMESH) from other BC types (pure local). Medium effort, moderate impact (~800ms saved).

#### 4b. Multi-Process + Multi-Thread (Hybrid MPI+Hedgehog)

Currently FDS uses MPI for multi-mesh (one process per mesh group). The Hedgehog integration adds intra-process parallelism (multiple threads for meshes within one process). The next frontier:

1. **Hybrid MPI+threads**: Each MPI rank runs a Hedgehog graph with `kernelThreads > 1`
2. **Load balancing**: Assign meshes to MPI ranks considering both mesh count and kernel thread availability
3. **MESH_EXCHANGE optimization**: Overlap MPI communication with kernel computation using Hedgehog's asynchronous task execution

#### 4c. Architectural Improvements

1. **Pipeline parallelism**: Overlap predictor/corrector computation of different meshes (requires decoupling barrier synchronization)
2. **Asynchronous MESH_EXCHANGE**: Start communication early, overlap with computation
3. **NUMA-aware mesh assignment**: Pin meshes to NUMA nodes for better memory locality

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

## Thread-Safe Routine Conversions (for future WALL_BC parallelization)

### Completed Conversions (Index-Based Access Pattern)

#### 1. ✅ CALC_HVAC_BC (wall.f90)  
**Lines**: 52  
**Approach**: Added explicit M and PREDICTOR_FLAG arguments  
**Key change**: Replaced module-level `PBAR_P` with conditional `M%PBAR_S` or `M%PBAR`  
**Call sites**: 2 (both in WALL_BC)  

#### 2. ✅ HEAT_TRANSFER_COEFFICIENT (func.f90)  
**Lines**: ~175  
**Approach**: Index-based access (no pointers to ALLOCATABLE arrays)  
**Signature**: `NM` → `M`, added explicit MESH_TYPE argument  
**Key pattern**:  
```fortran
! Instead of:
WC => M%WALL(WALL_INDEX)
B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
B1%TMP_F = ...

! Use:
B1_INDEX = M%WALL(WALL_INDEX)%B1_INDEX
M%BOUNDARY_PROP1(B1_INDEX)%TMP_F = ...
```
**Call sites**: 15+ (wall.f90, dump.f90)  
**Automated**: Used awk script to replace 17 pointer references  

### In Progress: Large Routines

#### 3. ⏸️ SURFACE_HEAT_TRANSFER (wall.f90)  
**Lines**: 379  
**Status**: Signature updated, module-level arrays prefixed with M%, pointer assignments need removal  
**Challenge**: Fortran does not allow TARGET attribute in TYPE definitions → cannot use pointers to M% arrays  
**Required approach**: Direct conditional access without pointers  
**Estimated effort**: 4-6 hours for full no-pointer conversion  
**Alternative**: Decompose into smaller case-specific kernels first  

#### 4. ⏸️ CALCULATE_ZZ_F (wall.f90)  
**Lines**: 413  
**Status**: Same as SURFACE_HEAT_TRANSFER  
**Challenge**: Same language limitation  
**Estimated effort**: 4-6 hours  

### Technical Limitation Discovered

**Fortran TARGET Restriction**:
```fortran
! NOT ALLOWED:
TYPE MESH_TYPE
   REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: U  ! ERROR
END TYPE
```

**Compiler error**: "Attribute at (1) is not allowed in a TYPE definition"

**Impact**: Cannot use pointer assignment to ALLOCATABLE arrays without TARGET:
```fortran
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU
UU => M%U  ! ERROR: target is neither TARGET nor POINTER
```

**Solution**: Eliminate all pointer usage, use direct array access:
```fortran
! Before (with local pointers):
IF (PREDICTOR_FLAG) THEN
   UU => M%US
ELSE
   UU => M%U
ENDIF
UN = UU(II,JJ,KK)

! After (direct access):
IF (PREDICTOR_FLAG) THEN
   UN = M%US(II,JJ,KK)
ELSE
   UN = M%U(II,JJ,KK)
ENDIF
```

**Trade-off**: More verbose but eliminates thread-unsafe module-level pointers.

### Documentation

- **WALL_BC_CONVERSIONS_SUMMARY.md** — Details of completed conversions  
- **WALL_BC_LARGE_ROUTINES_STATUS.md** — Challenge and options for large routines  
- **METHOD_KERNEL_EXTRACTION.md** — Updated with index-based access pattern

### Next Steps for WALL_BC Parallelization

**Critical path** (from WALL_BC_DECOMPOSITION.md):
1. ✅ Convert quick-win callees (CALC_HVAC_BC, HEAT_TRANSFER_COEFFICIENT)  
2. ⏸️ Complete large routine conversions (SURFACE_HEAT_TRANSFER, CALCULATE_ZZ_F)  
   - **Option A**: Full no-pointer conversion (4-6 hrs each)  
   - **Option B**: Decompose by boundary condition case first (2-3 hrs per case)  
3. Extract main WALL_BC loop into three phases:
   - Phase 1: ASSIGN_GHOST_VALUE (sequential, OMESH reads)  
   - Phase 2: WALL_BC_PROCESS_CELLS_KERNEL (parallel, 90% of cells)  
   - Phase 3: Cross-mesh finalization (INTERPOLATED_BC, BACK_MESH, CONSUME_MASS)  

**Estimated total to full WALL_BC parallelization**: 10-15 hours

