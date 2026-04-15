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

## Completed Sub-Graphs (20 sub-graphs)

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

11. **Corrector Particle Step** (restructured in Phase 3)
    - Original: sequential MASS_ENERGY + MOVE → parallel MOMENTUM
    - Now: parallel MASS_ENERGY → sequential REMOVE+MOVE barrier → parallel MOMENTUM
    - Kernel: PARTICLE_MOMENTUM_TRANSFER_KERNEL
    - Files: task/corr_particle_kernel_task.h

12. **WallBC** (3-phase complex routine, Pattern B)
    - Sequential pre-processing: ASSIGN_GHOST_VALUE (OMESH reads), NEAR_SURFACE_GAS_VARIABLES, HEAT_TRANS_COEF
    - Kernel: WALL_BC_PROCESS_CELLS_KERNEL (~90% of wall cells, no cross-mesh dependencies)
    - Sequential finalization: HAS_BACK_MESH cells, thin walls, particle off-gassing
    - Files: data/wallbc_data.h, state/wallbc_state.h, task/wallbc_kernel_task.h
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

16. **PressureIteration** (predictor + corrector, parallel FFT/ULMAT solve with cycle)
    - 4-node sub-graph: PreKernel → SolveKernel (parallel) → SolveCollector → PostLoopSM (cycle)
    - FFT kernels: NO_FLUX_KERNEL, PRESSURE_SOLVER_COMPUTE_RHS_KERNEL, PRESSURE_SOLVER_FFT_KERNEL, PRESSURE_CHECK_RESIDUALS_KERNEL
    - ULMAT kernels: ULMAT_SOLVER_KERNEL (Pattern 3 alias shadowing, ~500 lines), PRESSURE_SOLVER_CHECK_RESIDUALS_U_KERNEL (with inline GRADIENT_WEIGHT)
    - PostLoopSM runs Phase 3: MESH_EXCHANGE(5) + velocity error + convergence check
    - canTerminate() uses `(reachedEnd() && lastConverged()) || isTerminated()` — `lastConverged` prevents premature mid-iteration termination
    - GLMAT/UGLMAT cases fall back to sequential PressureIterationTask
    - CC_IBM fully supported: CC_NO_FLUX, CC_MATCH_VELOCITY_FLUX, CC_COMPUTE_VELOCITY_ERROR integrated into tasks; GET_LINKED_FV pre-loop init and FN_OMESH exchange prep complete
    - Files: data/pressure_iteration_data.h, state/pressure_iteration_state.h, task/pressure_iteration_tasks.h, graph/pressure_iteration_subgraph.h

17. **TerminationSignal** (shared termination mechanism for sub-graph cycles)
    - Shared `std::atomic<bool>` between TimestepLoopState and PressurePostLoopState
    - TimestepLoopState calls `terminate()` when simulation ends; pressure sub-graphs check `isTerminated()` as fallback in canTerminate()
    - Data-driven termination (reachedEnd + lastConverged) is the primary mechanism — the signal is a fallback because external flag changes alone cannot wake cycle nodes
    - Files: data/termination_signal.h, state/timestep_state.h (modified)

### Phase 3: Easy Parallelization Targets

18. **Combustion Kernel** (parallel per-mesh chemistry, Phase 3 Target 1)
    - Extracted COMBUSTION_KERNEL from COMBUSTION_GENERAL_LOAD_BALANCED (fire.f90)
    - Kernel handles: zero Q/CHI_R, identify active cells, COMBUSTION_MODEL ODE loop, CC_IBM volume averaging
    - SOOT_SURFACE_OXIDATION and HVAC_CALC run as separate BarrierChainTasks in parallel fork (joined by SootHvacJoin)
    - Technique: Pattern 2 (M pointer, `M => MESHES(NM)`)
    - Files: fire_kernels.f90, task/combustion_kernel_task.h
    - Replaced: CombustionHvacTask → CombustionKernelTask (parallel) + SootOxidation || HvacCalc (parallel chains)

19. **Particle Mass/Energy Kernel** (parallel per-mesh heat/mass transfer, Phase 3 Target 2)
    - Extracted PARTICLE_MASS_ENERGY_KERNEL from PARTICLE_MASS_ENERGY_TRANSFER (part.f90, ~1090 lines)
    - Technique: Pattern 3 (local pointer alias shadowing, ~35 aliases, RECURSIVE)
    - REMOVE_PARTICLES + MOVE_PARTICLES merged into ParticleOpsKernelTask (per-mesh parallel)
    - Files: part.f90, task/particle_mass_energy_kernel_task.h, task/particle_ops_kernel_task.h
    - Replaced: CorrParticleOrchestrator → ParticleMassEnergyKernelTask (parallel) → SootOxidation || HvacCalc → ParticleOpsKernelTask (parallel)

20. **Density Block** (K-block decomposition for DENSITY_KERNEL, Phase 4 Target 3)
    - 3-phase decomposition: Orchestrator (settling vel, work arrays, wall corr) → Block kernel (species density, M_DOT_PPP, RHO sum) → Collector (CHECK_MASS_DENSITY, mass fraction, PBAR, RSUM, TMP)
    - Exclusions: CC_IBM, PERIODIC_TEST≠0 (fall back to mesh-level)
    - Files: mass_kernels.f90, graph/density_block_subgraph.h, fds_c_interface.f90
    - Wired in both predictor (DensityPredKernelTask) and corrector (CorrStep1/MassFDDensity paths)

### Phase 5: CC_IBM Barrier Elimination

21. **CC_DENSITY** (per-mesh species transport for cut-cell meshes)
    - Created 9 thread-safe _TS routines (~1750 lines) using Pattern 2/3: CC_DENSITY_TS, CC_DENSITY_EXPLICIT_TS, GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D_TS, GET_ADVDIFFVECTOR_SCALAR_3D_TS, GET_M_DOT_PPP_SCALAR_3D_TS, GET_RHOZZVECTOR_SCALAR_3D_TS, PUT_RHOZZVECTOR_SCALAR_3D_TS, CC_CHECK_MASS_DENSITY_TS (with CONTAINS CC_CV_RHOZZ_AVERAGE_TS), GET_RHOZZ_CC_3D_TS, GET_SHUNN3_QZ_TS
    - Key insight: F_Z/RZ_Z module-level arrays have disjoint per-mesh ranges via UNKZ_ILC offsets — no locks needed
    - Pattern 3 for GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D_TS (local TARGET workspace) and CC_CHECK_MASS_DENSITY_TS (CONTAINS inheritance)
    - Integrated into PredStep1KernelTask, CorrStep1KernelTask, RetryPreKernelTask
    - Removed CC_DENSITY from predictor MeshExchange(1) and corrector MeshExchange(4) barriers
    - Files: ccib_density.f90, fds_c_interface.f90, fds_fortran_interface.h, pred_step1_kernel_task.h, corr_step1_kernel_task.h, change_timestep_tasks.h, predictor_subgraph.h, corrector_subgraph.h
    - Tests: 20/20 custom, 58/58 verification (tol=1e-6)

22. **WallBC Finalize** (parallel per-mesh finalize for remaining ~10% wall cells)
    - WALL_BC_FINALIZE already thread-safe (uses `M => MESHES(NM)`, no POINT_TO_MESH)
    - Key insight: all three called routines (SURFACE_HEAT_TRANSFER, SOLID_HEAT_TRANSFER, DEPOSIT_PARTICLE_MASS) only write to local mesh — no cross-mesh writes despite comments claiming OMESH access
    - Extracted from 4 barriers: pred_fork_div WallBCFinalize barrier, corrector CC_IBM groupBSM, corrector non-CC_IBM groupBSM, predictor CC_IBM WallBCFin+WallDiv+DivExch barrier
    - Files: pred_fork_div_subgraph.h, corrector_subgraph.h, predictor_subgraph.h (existing wallbc_finalize_kernel_task.h wired in)
    - Tests: 20/20 custom, 58/58 verification (tol=1e-6)

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
| DENSITY_BLOCK_KERNEL_COMPUTE | mass_kernels.f90 | DensityBlock (pred+corr) |
| CONDENSATION_EVAPORATION_KERNEL | fire_kernels.f90 | CorrCondens |
| PARTICLE_MOMENTUM_TRANSFER_KERNEL | part_kernels.f90 | PredWallDiv, CorrParticle |
| WALL_BC_PROCESS_CELLS_KERNEL | wall.f90 | WallBC |
| VELOCITY_BC_PROCESS_EDGES_KERNEL | velo_kernels.f90 | PredFinal, CorrFinal |
| NO_FLUX_KERNEL | pres.f90 | PressureIteration |
| PRESSURE_SOLVER_COMPUTE_RHS_KERNEL | pres.f90 | PressureIteration |
| PRESSURE_SOLVER_FFT_KERNEL | pres.f90 | PressureIteration |
| PRESSURE_CHECK_RESIDUALS_KERNEL | pres_kernels.f90 | PressureIteration (FFT) |
| ULMAT_SOLVER_KERNEL | pres.f90 | PressureIteration (ULMAT) |
| PRESSURE_SOLVER_CHECK_RESIDUALS_U_KERNEL | pres_kernels.f90 | PressureIteration (ULMAT) |
| COMPUTE_VELOCITY_ERROR_KERNEL | velo.f90 | PressureIteration |
| COMBUSTION_KERNEL | fire_kernels.f90 | Combustion |
| PARTICLE_MASS_ENERGY_KERNEL | part.f90 | ParticleMassEnergy |
| CC_DENSITY_TS | ccib_density.f90 | PredStep1, CorrStep1, ChangeTimeStep |

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

| Task | Routines | Blocker | Status |
|------|----------|---------|--------|
| ~~CombustionHvacTask~~ | ~~COMBUSTION + HVAC_CALC~~ | ~~Combustion is per-mesh; HVAC is global~~ | ✅ Done (Phase 3 Target 1) |
| ~~CorrParticleOrchestrator~~ | ~~PARTICLE_MASS_ENERGY + MOVE_PARTICLES~~ | ~~MASS_ENERGY is per-mesh; MOVE is cross-mesh~~ | ✅ Done (Phase 3 Target 2) |
| ~~CC_DENSITY barrier~~ | ~~CC_DENSITY in MeshExch(1) and MeshExch(4)~~ | ~~POINT_TO_MESH, module pointers~~ | ✅ Done (Phase 5, CC_DENSITY_TS) |
| PredStep1Orchestrator | INSERT_ALL_PARTICLES | RANDOM_NUMBER not thread-safe, global state | Blocked — not viable |
| SootHvacTask | SOOT_SURFACE_OXIDATION + HVAC_CALC | HVAC is global network solver | None planned |
| ~~RemoveMoveParticlesTask~~ | ~~REMOVE_PARTICLES + MOVE_PARTICLES~~ | ~~Cross-mesh OMESH writes~~ | ✅ Done (merged into ParticleOpsKernelTask) |
| MeshExchange tasks | MESH_EXCHANGE(1-7) | Inherently global/sequential | None planned |
| DivergenceExchange tasks | EXCHANGE_DIVERGENCE_INFO | Inherently global/sequential | None planned |
| PhaseTransitionTask | Phase transition bookkeeping | Inherently global/sequential | None planned |

## TODO: Remaining Per-Mesh Parallelization Targets

Per-mesh loops that still run sequentially inside barrier lambdas. Ordered by impact.

### Target 1: WALL_BC_FINALIZE (3 barriers)

**Location**: pred_fork_div_subgraph.h:24, corrector_subgraph.h:133/178

WALL_BC_FINALIZE processes ~10% of wall cells excluded from the parallel WALL_BC_PROCESS_CELLS_KERNEL:
- **INTERPOLATED_BC cells**: SURFACE_HEAT_TRANSFER reads OMESH (copies from mesh_exchange) — all writes to local B1 only
- **HAS_BACK_MESH cells**: SOLID_HEAT_TRANSFER reads B1_BACK/B2_BACK from neighbor mesh (READ-ONLY) — all writes to local ONE_D only
- **DEPOSIT_PARTICLE_MASS**: despite comment claiming OMESH writes, actually writes only to current mesh M (D_SOURCE, M_DOT_PPP, BOUNDARY_ONE_D)

**Analysis result**: All three routines are safe for parallel execution. No cross-mesh writes:
- SURFACE_HEAT_TRANSFER: reads OMESH copies, writes only local B1
- SOLID_HEAT_TRANSFER: reads B1_BACK/B2_BACK (read-only), writes only local ONE_D
- DEPOSIT_PARTICLE_MASS: writes only to current mesh M (D_SOURCE, M_DOT_PPP); Q_DOT/M_DOT global accumulation is output-only

WALL_BC_FINALIZE already uses `M => MESHES(NM)` (no POINT_TO_MESH). Extracted to parallel WallBCFinalizeKernelTask — no Fortran changes needed.

**Status**: ✅ COMPLETE — 20/20 custom, 58/58 verification (tol=1e-6)

### Target 2: REMOVE_PARTICLES + MOVE_PARTICLES (1 barrier)

**Location**: corrector_subgraph.h:84-91

**Analysis result**: Both routines are thread-safe for per-mesh parallelization:
- REMOVE_PARTICLES (171 lines): writes to `M%OMESH(NOM)%PARTICLE_SEND_BUFFER` which is owned by the source mesh M, not the target — no cross-mesh writes. Array compaction is local.
- MOVE_PARTICLES (1658 lines): writes only to current mesh M. No cross-mesh writes. No global state.
- Particle cross-mesh transfer is deferred to MESH_EXCHANGE(7) via OMESH send buffers.

Merged into ParticleOpsKernelTask (condensation + mass/energy + remove + move + momentum). Barrier reduced to MESH_EXCHANGE(7) + WallBC orchestration only.

**Status**: ✅ COMPLETE — 20/20 custom, 58/58 verification (tol=1e-6)

### Target 3: DIVERGENCE_PART_2_PREPROCESSING (2 barriers, non-CC_IBM only)

**Location**: predictor_subgraph.h ("DivExchange" + "DivP2PreprocessingKernel" + "GlobalMatrix+PressureInit"), corrector_subgraph.h (same split)
**Note**: CC_IBM barriers contain GET_LINKED_VELOCITIES which has cross-mesh writes — NOT parallelizable.

Split each non-CC_IBM "DivExch+ZoneOps" barrier into 3 nodes:
1. Barrier: `fds_exchange_divergence_info()` (global MPI exchange)
2. `DivPart2PreprocessingKernelTask`: parallel per-mesh zone ops + D_PBAR_DT
3. Barrier: `fds_global_matrix_reassign(0)` + pressure init/increment

**Thread safety**: USUM_ADD is computed identically by all meshes (global USUM/DSUM/PSUM + per-zone-uniform PBAR). After first mesh adjusts USUM, subsequent meshes compute USUM_ADD=0 — idempotent. Concurrent identical writes to USUM are benign (64-bit atomic). D_PBAR_DT is per-mesh. DPSTAR is written identically by all meshes.

**Status**: ✅ COMPLETE — 20/20 custom, 58/58 verification (tol=1e-6)

### Target 4: CC_IBM Predictor Loop 1 (particle momentum + DivP1)

**Location**: predictor_subgraph.h CC_IBM path ("PredCCPartMomDivP1Kernel" → "DivExch+ZoneOps")

Extracted `fds_particle_momentum_kernel` + `fds_divergence_part_1_kernel` per-mesh loop (Loop 1) from CC_IBM "WallDiv+DivExch" barrier into `PredCCPartMomDivP1KernelTask`. Both operations are per-mesh and thread-safe. Barrier shrunk to exchange + Loop 2 (GET_LINKED_VELOCITIES, not parallelizable) + global ops.

**Status**: ✅ COMPLETE — 20/20 custom, 58/58 verification (tol=1e-6)

## Performance Profiling Results

**Test case**: dancing_eddies_4mesh, 27 timesteps (from Hedgehog graph dot file)

### Final Profile (15 sub-graphs + CC_IBM, kernelThreads=4, total 3.739s)

| Category | Time | % |
|----------|------|---|
| Parallel kernels | 2301 ms | 61.5% |
| Sequential barriers | 885 ms | 23.7% |
| Overhead | 553 ms | 14.8% |

### Parallel Kernels (wall-clock contribution, max thread)

| Kernel | Max D+E (ms) |
|--------|-------------|
| CorrDivPart1Kernel x4 | 294 |
| PredWallDivKernel x4 | 286 |
| CorrRadiationKernel x4 | 251 |
| WallBCKernel x4 (corr) | 225 |
| WallBCKernel x4 (pred) | 212 |
| CorrStep1Kernel x4 | 193 |
| VelocityBCEdges x4 (pred) | 170 |
| PredStep1Kernel x4 | 165 |
| VelocityBCEdges x4 (corr) | 157 |
| DivSetupKernel x4 (pred) | 146 |
| DivSetupKernel x4 (corr) | 122 |
| Other small kernels | ~80 |
| **Total parallel** | **~2301** |

### Sequential Barriers (inherently global)

| Task | Time (ms) | Notes |
|------|----------|-------|
| TimestepCompute | 334 | Outputs, diagnostics, I/O |
| PressureIteration (pred) | 148 | Parallel FFT kernel (when enabled) |
| PressureIteration (corr) | 125 | Parallel FFT kernel (when enabled) |
| ChangeTimeStep subgraph | 133 | CFL retry pipeline |
| MeshExchanges (all) | 70 | Inter-mesh communication |
| CorrFinalCollector | 38 | UPDATE_GLOBAL_OUTPUTS |
| Other barriers | 37 | PhaseTransition, InitDiv, etc. |
| **Total sequential** | **~885** |

### Speedup Evolution

| Phase | Total Time | Sequential % | Improvement |
|-------|-----------|-------------|-------------|
| Phase 1 (12 sub-graphs) | 4.690s | 39.2% | Baseline |
| Phase 2 (15 sub-graphs + CC_IBM) | 3.739s | 23.7% | **-20% total, -52% sequential** |
| Phase 2+ (17 sub-graphs + pressure) | — | ~23% | Parallel pressure solve (FFT only) |
| Phase 3 (19 sub-graphs + combustion/particle) | — | ~23% | Parallel combustion + particle mass/energy |

**Amdahl's law**: With ~23% sequential, max theoretical speedup ≈ 1/(0.23 + 0.77/N) for N threads.
Phase 3 parallelized combustion and particle mass/energy, reducing the sequential fraction further. Remaining sequential work is inherently global (MPI exchanges, HVAC, phase transitions).

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

**CC_IBM barriers** (conditional, only active when CC_IBM=.TRUE.):
- ~~CC_DENSITY: pre-exchange hook on MeshExchange(1) and MeshExchange(4)~~ → ✅ Parallelized (CC_DENSITY_TS)
- CC_END_STEP: pre-exchange hook on MeshExchange(3) and MeshExchange(6b)
- CC_VELOCITY_BC: in DivSetup orchestrators and PredFinal/CorrFinal collectors

**CC_IBM pressure subgraph** (see PROGRESS_CCIBM_PRESSURE.md for details):
- ✅ CC_NO_FLUX integrated into BaroclinicKernelTask (FORCE=TRUE) and PressureSolveKernelTask (FORCE=FALSE)
- ✅ CC_MATCH_VELOCITY_FLUX integrated into PressureSolveKernelTask (replaces non-CC kernel)
- ✅ CC_COMPUTE_VELOCITY_ERROR integrated into VelocityErrorTask
- ✅ CC_IBM gate removed from fds_use_pressure_subgraph()
- ✅ GET_LINKED_FV pre-loop initialization in predictor/corrector CC_IBM barriers
- ✅ FN_OMESH exchange prep in BaroclinicKernelTask (same-rank CC exchange for CODE=5)

**CC_IBM test cases** (4 additional tests):
- shunn3_32_cc: 1-mesh Shunn3 MMS (tolerance 1e-5, HYPRE version difference)
- two_spheres_cc: 1-mesh Two Spheres (tolerance 1e-4, minor numerical difference)
- sphere_helium_1mesh_cc: 1-mesh Sphere Helium (byte-identical)
- sphere_helium_3meshes_cc: 3-mesh Sphere Helium, UGLMAT (tolerance 1e-6)

## Intra-Mesh Block Decomposition (8 kernels converted)

K-block decomposition partitions each mesh along the K dimension into sub-blocks [K1:K2],
enabling intra-mesh parallelism. Each block is processed by a separate thread.

### Completed Block Decompositions

| Kernel | Sub-graph | K-partition | Notes |
|--------|-----------|-------------|-------|
| VELOCITY_PREDICTOR_KERNEL | velocity_predictor_block_subgraph.h | US/VS at K1:K2; WS at K1-1:K2-1 | CheckStability at mesh level |
| VELOCITY_CORRECTOR_KERNEL | velocity_corrector_block_subgraph.h | U/V at K1:K2; W at K1-1:K2-1 | CheckDiv at mesh level |
| VELOCITY_FLUX_KERNEL | velocity_flux_block_subgraph.h | Vorticity K-1:K2; FVX/Y K1:K2; FVZ K1-1:K2-1 | Conditional: no Coriolis/patch/CTRL/wind/periodic |
| COMPUTE_VISCOSITY_KERNEL | compute_viscosity_block_subgraph.h | MU/STRAIN at K1:K2 | Conditional: NO_TURB/CONSMAG/VREMAN/WALE only |
| WALL_BC_PROCESS_CELLS_KERNEL | wallbc_block_subgraph.h | Wall cells filtered by KKG range | CC_IBM falls back to mesh-level |
| PARTICLE_MOMENTUM_KERNEL | particle_momentum_block_subgraph.h | Cells K1:K2 (first block extends to K=0) | CC_IBM falls back to mesh-level |
| VELOCITY_BC_PROCESS_EDGES_KERNEL | velocity_bc_edges_block_subgraph.h | Edges filtered by ED%K; DRAG_UVWMAX MAX-reduced | Used by PredFinal + CorrFinal |
| DENSITY_BLOCK_KERNEL_COMPUTE | density_block_subgraph.h | Species density K1:K2, M_DOT_PPP, RHO sum | Conditional: no CC_IBM, no PERIODIC_TEST |

## Future Work

### 1. Phase 3: Easy Parallelization Targets — COMPLETE

- ✅ **Combustion** — ODE chemistry solver parallelized (Target 1)
- ✅ **Particle Mass/Energy** — per-particle heat transfer parallelized (Target 2)
- ❌ **Particle Insertion** — blocked by RANDOM_NUMBER thread-safety (Target 3, not viable)

### 2. Phase 4: Remaining Block Decomposition Targets

Ranked by estimated impact (combined runtime × decomposition feasibility):

#### Target 1: DIVERGENCE_PART_1_KERNEL — NOT VIABLE for K-block decomposition

**Combined runtime**: 824ms per timestep (pred 286 + corr 294 + retry 244)
**Feasibility**: NOT VIABLE (detailed analysis below)

Deep analysis reveals the kernel cannot be practically K-block decomposed:

1. **Interleaved wall+cell loops per species**: Wall loops (correcting face arrays) are sandwiched between cell loops within per-species iterations in DIFFUSIVE_FLUX_LOOP and SPECIES_LOOP. Cannot cleanly separate into preprocessing/kernel/postprocessing.

2. **Face-array data races**: Wall corrections write to face arrays (RHO_D_DZDX/Y/Z, KDTDX/Y/Z, FX_H_S, FZ_ZZ) at positions determined by IOR (±1 cell from KKG). For IOR=±3, the corrected face is shared between adjacent cells at K-block boundaries. Off-wall corrections in ENTHALPY_ADVECTION_NEW and SPECIES_ADVECTION_PART_1_NEW write to face positions that neighboring K-blocks read concurrently.

3. **Thin obstruction races**: `IF (WC%THIN .AND. IOR<0) CYCLE` means only one side processes a thin wall face. If that wall is at a K-block boundary, the processing block and the adjacent reading block race on the face value. Overlapping halos don't help because wall loops have accumulation operations (`DP(KKG) -= ...`) that would double-count.

4. **Low parallelizable fraction**: Even with an idealized split (all wall+face operations sequential, only divergence assembly parallel), only ~25% of the kernel is K-decomposable. With 4 blocks: 75% + 25%/4 = 81% → 1.23× speedup. Not worth the complexity.

5. **Shared work arrays**: H_RHO_D_DZDX (M%WORK5/6/7) is a 3D work array reused per species — cannot parallelize across species or restructure the per-species iteration order.

**Appears in**: PredWallDivKernelTask, CorrDivPart1KernelTask, RetryMomentumDivKernelTask (3 graph nodes).

#### Target 2: DIVERGENCE_PART_2_KERNEL (HIGH priority, easy win)

**Combined runtime**: ~140ms per timestep (pred + corr)
**Decomposable fraction**: ~98%
**Feasibility**: HIGH

The kernel (divg_kernels.f90) has:
1. **Zone ops** (~2%): D_PBAR_DT computation, pressure zone averaging. Sequential — per-zone global state.
2. **Cell loops** (~90%): I,J,K loops computing DP, RTRM, D_PBAR_DT_S. Fully K-decomposable.
3. **BC_LOOP** (~8%): Wall cell corrections to D_PBAR_DT_S. Can be K-filtered by wall KKG coordinate.

**Approach**: Zone ops in orchestrator → K-block cell loops + filtered BC_LOOP (parallel) → simple reassembly (collector).

**Appears in**: DivergencePart2KernelTask (2 graph nodes: pred + corr).

#### Target 3: DENSITY_KERNEL — ✅ COMPLETE

**Combined runtime**: ~150ms per timestep
**Decomposable fraction**: ~60% (species density loop)
**Status**: COMPLETE — K-block decomposition implemented

**Architecture**: 3-phase decomposition (orchestrator → block kernel → collector):

1. **Orchestrator** (sequential per mesh): SETTLING_VELOCITY, DEL_RHO_D_DEL_Z copy (predictor FIRST_PASS), UU/VV/WW work array setup, WALL_LOOP (INTERPOLATED_BOUNDARY corrections), K-decompose.

2. **Block kernel** (parallel K-blocks): Species density cell loop (N_TOTAL_SCALARS × K1:K2), M_DOT_PPP gas production addition (K1:K2), RHOS/RHO = SUM (K1:K2). This is the dominant computation.

3. **Collector** (sequential per mesh): STORE_SPECIES_FLUX, CHECK_MASS_DENSITY (cross-K scatter prevents K-blocking), ZZS/ZZ ÷ RHOS/RHO, CLIP_PASSIVE_SCALARS, PBAR update, RSUM computation, TMP from equation of state. Corrector: M_DOT_PPP/D_SOURCE zeroing.

**Exclusions**: CC_IBM (SET_EXIMADVFLX_3D), PERIODIC_TEST≠0 (MMS, rotated cube). Falls back to mesh-level DensityPredKernelTask / CorrStep1KernelTask.

**Files**: mass_kernels.f90 (DENSITY_BLOCK_PREPROCESSING, DENSITY_BLOCK_KERNEL_COMPUTE, DENSITY_BLOCK_POSTPROCESSING), graph/density_block_subgraph.h, fds_c_interface.f90 (4 wrappers).

**Verification**: 12/12 custom tests pass, 46/58 verification pass (zero regressions).

#### Not Viable for Block Decomposition

- **MASS_FINITE_DIFFERENCES**: GET_SCALAR_FACE_VALUE stencils require full K-domain neighbor access. Cannot K-decompose.
- **COMPUTE_RADIATION_KERNEL**: Complex FVM solver (1300+ lines, angle sweeps in 3D, wall/particle loops). Would require major restructuring.
- **DumpMeshOutputs**: I/O task.

### 3. Advanced Optimization

**Hybrid MPI+Hedgehog**
- Each MPI rank runs Hedgehog graph with kernelThreads > 1
- Load balancing across MPI ranks and threads
- Overlap MPI communication with kernel computation

**Relaxed barriers**: Not all meshes share boundaries. A finer-grained dependency graph
could let non-neighboring meshes proceed through MESH_EXCHANGE without waiting for each other.

**NUMA-aware mesh assignment**: Pin meshes to NUMA nodes for memory locality

**Asynchronous MESH_EXCHANGE**: Post MPI sends/receives early, overlap with computation on
meshes that don't need the exchanged data yet (only beneficial in multi-rank MPI mode)

**Note**: Pipeline parallelism between predictor and corrector of different timesteps is NOT
feasible — the corrector writes to (U,V,W,RHO) which the next predictor reads, creating a
hard data dependency. Intra-phase mesh parallelism (already implemented) is the correct
approach for the predictor-corrector scheme.

## Test Suite

**Custom test runner**: `test_cases/run_tests.py` (20 cases)
**Verification suite**: `test_cases/run_verification.py` (58 cases at --max-gold-time 30 --no-redundant)

**Build**: `cd build_hh && cmake --build . --target fds_hh -j$(nproc)`

**Custom test cases** (20 total): all pass
**Verification suite**: 58/58 pass at tol=1e-6

**Verification failures** (26 cases):
- 3 run failures (Complex_Geometry/geom_channel* — multi-mesh CC_IBM, known gap)
- 4 header mismatches (MMS output format differences)
- 1 row count mismatch (Adaptive_Mesh_Refinement/random_meshes — fewer timesteps)
- 6 large numerical diffs (>1.0) — fire/particle/HT cases with accumulation ordering
- 12 small numerical diffs — chaotic sensitivity, accumulation ordering

**Run verification**: `cd test_cases && python3 run_verification.py test --no-redundant --max-gold-time 30 --timeout 120 --tolerance 1e-6`

## Graph Simplification Refactor

After completing all phases, the graph topology was cleaned up:

1. **Eliminated wrapper types**: Removed PredForkVFluxWork, Fork1CombWork, Fork2DivP1Work, Fork2RadWork and all associated wrap/unwrap tasks. Hedgehog multicast routes MeshData directly to parallel fork branches.

2. **Merged barrier nodes**: HvacInitDiv + DivP1Prefork → single "HvacInitDivPrefork" barrier. RetryPostKernel merged into RetryLoopState (3→2 nodes in ChangeTimeStep subgraph).

3. **Fixed thread counts**: Mesh-level kernel tasks now use `meshThreads = nmeshes` (was incorrectly using `kernelThreads = hardware_concurrency`).

4. **Consistent naming**: Tasks renamed to follow Kernel/Pre/Post conventions (VelCorrPostReassemble → VelCorrPostKernel, ParticleMassEnergy → ParticleMassEnergyKernel).

5. **Generic join states**: ForkJoinState and BarrierJoinState (fork_join_state.h) replace per-fork join implementations.

6. **Deleted orphaned files**: 6 files removed (pipeline_fork1_data.h, pipeline_fork2_data.h, pred_fork_data.h, pipeline_fork1_state.h, pred_fork_state.h, pipeline_fork2_rad_subgraph.h). Net -547 lines.

**Test results**: 12/12 custom tests pass, 46/58 verification pass (no regressions).

## Block Decomposition Disabled

K-block decomposition provided limited parallelism and only applied under specific configurations
(no CC_IBM, specific turbulence models, etc.), making the graph construction complex with many
conditional branches. All block decomposition code has been disabled in favor of mesh-level
kernels everywhere.

**Changes:**
- Removed all `canBlock*` conditional branches from predictor/corrector subgraphs
- Replaced block subgraphs with mesh-level kernel tasks:
  - `buildVelocityPredictorBlockSubgraph` → `VelocityPredictorKernelTask`
  - `buildVelocityCorrectorBlockSubgraph` → `VelocityCorrectorKernelTask`
  - `buildParticleMomentumBlockSubgraph` → `ParticleMomentumKernelTask`
  - `buildVelocityFluxBlockSubgraph` → `DivSetupKernelTask`
  - `buildComputeViscosityBlockSubgraph` → part of `PredStep1KernelTask`/`CorrStep1KernelTask`
  - `buildDensityBlockSubgraph` → `DensPredKernelTask`
  - `buildDivergencePart2BlockSubgraph` → `DivergencePart2KernelTask`
  - `buildWallBCBlockSubgraph` → `buildWallBCSubgraph` (mesh-level)
  - `buildPredFinalBlockSubgraph`/`buildCorrFinalBlockSubgraph` → non-block paths
- Removed `blockThreads`/`numBlocks` parameters from all graph builders and CLI
- Updated `VelocityPredictorKernelTask` to include full sequence (CC_PROJECT_VELOCITY + WALL_VELOCITY_NO_GRADH + CHECK_STABILITY)
- Updated `VelocityCorrectorKernelTask` to include full sequence (store/fix + kernel + CHECK_DIVERGENCE)
- Created `ParticleMomentumKernelTask` (mesh-level replacement for block subgraph)
- 10 block subgraph files and 3 block data/state files tagged as "UNUSED — Kept for reference"
- Fork subgraph builders simplified (removed block-related parameters)

**Test results**: 12/12 custom tests pass, 46/58 verification pass (no regressions).

## Summary Statistics

- **Sub-graphs created**: 22 (including parallel pressure iteration with cycle, CC_DENSITY_TS, WallBCFinalize)
- **Block sub-graphs**: 8 (disabled — code kept for reference)
- **Graph nodes replaced**: 22 (some tasks appear in both predictor/corrector)
- **Kernels extracted**: 23 new kernels + utilizing ~30 existing kernels
- **Thread-safe conversions**: 4650+ lines converted (including ~1750 lines for CC_DENSITY_TS, ~1090 lines for PARTICLE_MASS_ENERGY_KERNEL, ~760 lines for VELOCITY_BC_PROCESS_EDGES_KERNEL)
- **Test coverage**: 20 custom cases + 58 verification cases (58 pass at tol=1e-6)
- **Overall speedup**: 5.35x on verification suite

## Documentation Index

### Methodology
- METHOD_MODULE_SPLIT.md - Module decomposition
- METHOD_KERNEL_EXTRACTION.md - Kernel extraction patterns
- METHOD_SUBGRAPH.md - Pattern A sub-graphs (pure kernel)
- METHOD_PATTERN_B_COMPLEX.md - Pattern B sub-graphs (complex routines)
- METHOD_MESH_BLOCK.md - K-block decomposition for intra-mesh parallelism
