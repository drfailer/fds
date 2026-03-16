# Mesh Block Decomposition Progress

Tracking the classification and conversion of Hedgehog kernel tasks for intra-mesh block parallelism.

## Methodology

See [METHOD_MESH_BLOCK.md](METHOD_MESH_BLOCK.md) for the step-by-step procedure.

## Infrastructure

- [ ] Create `MeshBlockData` token (`data/mesh_block_data.h`)
- [ ] Create `MeshBlockDecomposeState` (`state/mesh_block_decompose_state.h`)
- [ ] Create `MeshBlockReassembleState` (`state/mesh_block_reassemble_state.h`)
- [ ] Add `fds_get_ibar`/`fds_get_jbar`/`fds_get_kbar` C bindings
- [ ] Verify infrastructure with first kernel conversion

## Kernel Classification Legend

- **Mesh Block**: All cell loops can be restricted to `[I1:I2, J1:J2, K1:K2]` sub-ranges
- **Mesh**: Must operate on entire mesh (wall loops, particle loops, zone loops, global writes)
- **Mixed**: Contains both block-decomposable and mesh-level operations; may need splitting

## Predictor Sub-Graph

### PredStep1KernelTask

**Graph location:** Predictor entry (after PredStep1Orchestrator scatters)
**File:** `task/pred_step1_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `COMPUTE_VISCOSITY_KERNEL` | velo_kernels.f90:921 | **Mesh** | CONTAINS subroutines (WALE, Deardorff, etc.), wall loops (IW), turbulence model dispatch |
| 2 | `MASS_FINITE_DIFFERENCES_NEW_KERNEL` | mass_kernels.f90:23 | **Mesh** | Cell loops (I,J,K) + wall face correction loops (WALL_LOOP_2, WALL_LOOP_3) |

**Task classification:** Mesh — both kernels have wall loops that iterate over all wall cells.
**Status:** DONE (classified, no conversion needed)

---

### DensityPredKernelTask

**Graph location:** After PredStep1KernelTask
**File:** `task/density_pred_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DENSITY_KERNEL` | mass_kernels.f90:350 | **Mesh** | I,J,K loops + wall loops + zone loops (PBAR updates) + CHECK_MASS_DENSITY (mesh-level flags) |

**Task classification:** Mesh — zone-dependent pressure updates, wall boundary corrections, CHECK_MASS_DENSITY mesh flags.
**Status:** DONE (classified, no conversion needed)

---

### DivSetupKernelTask (predictor instance)

**Graph location:** After MeshExchange(1), optional CC_IBM orchestrator
**File:** `task/div_setup_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `VELOCITY_FLUX_KERNEL` | velo_kernels.f90:325 | **Mesh** | Large kernel (590 lines) with CONTAINS, wall loops (DO IW at line 698), species loops, cell loops |

Note: `fds_set_baroclinic_false` and `fds_viscosity_bc_kernel` are called before this kernel but are sequential pre-processing (run in orchestrator or before dispatch).

**Task classification:** Mesh — wall loops interleaved with cell loops, CONTAINS subroutines with host association.
**Status:** DONE (classified, no conversion needed)

---

### PredWallDivKernelTask

**Graph location:** After WallBC sub-graph
**File:** `task/pred_wall_div_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MOMENTUM_TRANSFER_KERNEL` | part_kernels.f90:24 | **Mesh Block** | Pure I,J,K loop (0:IBAR, 0:JBAR, 0:KBAR). Block kernel exists. |
| 2 | `DIVERGENCE_PART_1_KERNEL` | divg_kernels.f90:27 | **Mesh** | Cell loops + multiple wall loops + per-mesh local accumulators (D_SUM_LOC, P_SUM_LOC, U_SUM_LOC) |

**Task classification:** Mixed — PARTICLE_MOMENTUM is Mesh Block (lightweight), DIVERGENCE_PART_1 is Mesh (dominates runtime). No conversion worthwhile.
**Status:** DONE (classified, no conversion needed)

---

### DivergencePart2KernelTask (predictor instance)

**Graph location:** After DivergenceExchange barrier
**File:** `task/divergence_part2_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_2_KERNEL` | divg_kernels.f90:1396 | **Mesh** | Zone loops (D_PBAR_DT, pressure averaging) + I,J,K cell loops + wall loops (BC_LOOP) |

**Task classification:** Mesh — zone-dependent pressure calculations must run before cell loops; wall loops for boundary corrections.
**Status:** DONE (classified, no conversion needed)

---

### PressureSolveKernelTask (predictor instance)

**Graph location:** Inside PressureIteration sub-graph cycle
**File:** `task/pressure_iteration_tasks.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `NO_FLUX_KERNEL` | pres_kernels.f90:737 | **Mesh** | Wall loop (IW=1 to N_EXTERNAL+N_INTERNAL) |
| 2 | `PRESSURE_SOLVER_COMPUTE_RHS` | pres_kernels.f90:24 | **Mesh** | Wall pre-processing + I,J,K RHS assembly loop |
| 3 | `PRESSURE_SOLVER_FFT` | pres_kernels.f90:308 | **Mesh** | FFT solver (tridiagonal + spectral transform) — global mesh operation |
| 4 | `PRESSURE_SOLVER_CHECK_RESIDUALS` | pres_kernels.f90:477 | **Mesh** | I,J,K residual loop + baroclinic correction with wall checks |

Or ULMAT variant:
| 3a | `ULMAT_SOLVER_KERNEL` | pres.f90 | **Mesh** | PARDISO/HYPRE direct solve per zone — global mesh operation |
| 4a | `PRESSURE_SOLVER_CHECK_RESIDUALS_U_KERNEL` | pres_kernels.f90:586 | **Mesh** | I,J,K loops + inline GRADIENT_WEIGHT |

**Task classification:** Mesh — FFT/ULMAT solvers operate on the full mesh; wall loops in NO_FLUX and COMPUTE_RHS.
**Status:** DONE (classified, no conversion needed)

---

### VelocityPredictorKernelTask

**Graph location:** After PressureIteration, before ChangeTimeStep
**File:** `task/velocity_predictor_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `VELOCITY_PREDICTOR_KERNEL` | velo_kernels.f90:112 | **Mesh Block** | Three pure I,J,K loops for US, VS, WS (staggered grid). No cross-cell deps. |
| 2 | `CHECK_STABILITY_KERNEL` | velo_kernels.f90:1253 | **Mesh** | I,J,K CFL/VN loops + wall loop, all with min/max reductions to mesh scalars |

**Task classification:** Mixed — VELOCITY_PREDICTOR block-decomposed along K, CHECK_STABILITY runs at mesh level after reassembly.
**Status:** DONE

**Implementation:**
- Block kernel: `VELOCITY_PREDICTOR_BLOCK_KERNEL(M, DT, K1, K2)` in velo_kernels.f90
- Sub-graph: `graph/velocity_predictor_block_subgraph.h` — Decompose → BlockKernel → Reassemble → CheckStability (if !skipCFL)
- CC_IBM path: uses original mesh-level `VelocityPredictorKernelTask` (CFL deferred to CC collector)
- Verified: 12/12 custom pass, verification pending

---

### VelocityBCEdgesTask (predictor instance, inside PredFinal sub-graph)

**Graph location:** PredFinal sub-graph kernel
**File:** `task/velocity_bc_edges_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `MATCH_VELOCITY_KERNEL` | velo_kernels.f90:2238 | **Mesh** | External wall loop (IW=1 to N_EXTERNAL_WALL_CELLS) |
| 2 | `VELOCITY_BC_PROCESS_EDGES_KERNEL` | velo_kernels.f90:1430 | **Mesh** | Edge loop with wall BC application, wall loops |

Note: `fds_velocity_bc_preprocessing` runs before the kernel but is sequential pre-processing.

**Task classification:** Mesh — both kernels iterate over wall/edge indices, not (I,J,K) grid.
**Status:** DONE (classified, no conversion needed)

---

### RetryMomentumDivKernelTask (inside ChangeTimeStep sub-graph)

**Graph location:** ChangeTimeStep CFL retry loop
**File:** `task/change_timestep_tasks.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MOMENTUM_TRANSFER_KERNEL` | part_kernels.f90:24 | **Mesh Block** | Same as PredWallDiv — block kernel exists |
| 2 | `DIVERGENCE_PART_1_KERNEL` | divg_kernels.f90:27 | **Mesh** | Same as PredWallDiv — cell+wall loops + local accumulators |

**Task classification:** Mixed — same as PredWallDiv (#6). PARTICLE_MOMENTUM lightweight; DIVERGENCE_PART_1 dominates.
**Status:** DONE (classified, no conversion needed)

---

## Corrector Sub-Graph

### CorrStep1KernelTask

**Graph location:** Corrector entry
**File:** `task/corr_step1_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `COMPUTE_VISCOSITY_KERNEL` | velo_kernels.f90:921 | **Mesh** | Same as PredStep1 — CONTAINS, wall loops, turbulence dispatch |
| 2 | `MASS_FINITE_DIFFERENCES_NEW_KERNEL` | mass_kernels.f90:23 | **Mesh** | Same as PredStep1 — cell loops + wall correction |
| 3 | `DENSITY_KERNEL` | mass_kernels.f90:350 | **Mesh** | Same as DensityPred — cell loops + CHECK_MASS_DENSITY |

**Task classification:** Mesh — all three kernels have wall loops and/or zone-level operations.
**Status:** DONE (classified, no conversion needed)

---

### DivSetupKernelTask (corrector instance)

**Graph location:** After MeshExchange(4), optional CC_IBM orchestrator
**File:** `task/div_setup_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `VELOCITY_FLUX_KERNEL` | velo_kernels.f90:325 | **Mesh** | Same as predictor instance — wall loops, CONTAINS |

Note: corrector also calls `fds_agglomeration(dt, nm)` after the kernel — this is a sequential operation.

**Task classification:** Mesh — wall loops interleaved with cell loops, CONTAINS with host association.
**Status:** DONE (classified, no conversion needed)

---

### CombustionKernelTask

**Graph location:** After DivSetup, before SootHvac barrier
**File:** `task/combustion_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `COMBUSTION_GENERAL_KERNEL` | fire.f90:532 | **Mesh** | Builds active cell list, calls COMBUSTION_MODEL per cell, STOP_STATUS global flag, CONTAINS with host association, CC_IBM cut-cell loops |

**Task classification:** Mesh — builds cell index list (not I,J,K grid loop), global STOP_STATUS flag, CONTAINS subroutines with host association.
**Status:** DONE (classified, no conversion needed)

---

### CorrCondensKernelTask

**Graph location:** After SootHvac barrier
**File:** `task/corr_condens_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `CONDENSATION_EVAPORATION_KERNEL` | fire_kernels.f90:1138 | **Mesh** | Species loop containing I,J,K cell loop + wall loop (WALL_LOOP), sharing ZZ_INTERIM/RHO_INTERIM/TMP_INTERIM. Wall loop modifies gas cells at arbitrary (I,J,K) locations. Cannot decompose — wall and cell loops interleaved within species iteration. |

**Task classification:** Mesh — wall loops interleaved with cell loops within SPEC_LOOP, shared interim arrays prevent clean K-decomposition.
**Status:** DONE (classified, no conversion needed)

---

### ParticleMassEnergyKernelTask

**Graph location:** After CorrCondens, before RemoveMoveParticles barrier
**File:** `task/particle_mass_energy_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MASS_ENERGY_TRANSFER_KERNEL` | part_kernels.f90 (wrapper) | | Particle loop (DO IP=1,NLP) — iterates over Lagrangian particles |

**Task classification:** Mesh — particle loop (DO IP=1,NLP) iterates over Lagrangian particles, not cells.
**Status:** DONE (classified, no conversion needed)

---

### CorrParticleKernelTask

**Graph location:** After RemoveMoveParticles barrier
**File:** `task/corr_particle_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MOMENTUM_TRANSFER_KERNEL` | part_kernels.f90:27 | **Mesh Block** | Pure I,J,K loop (0:IBAR, 0:JBAR, 0:KBAR), per-cell force accumulation. CC_IBM CUTFACE_VELOCITIES is mesh-level. |

**Task classification:** Mesh Block — entire kernel is a single I,J,K loop with no cross-cell dependencies.
**Status:** DONE

**Implementation:**
- Block kernel: `PARTICLE_MOMENTUM_BLOCK_KERNEL(M, DT, K1, K2)` in part_kernels.f90
- K-partition: first block extends to K=0 (covers full 0:KBAR range); cells K=K1:K2 for other blocks
- Sub-graph: `graph/particle_momentum_block_subgraph.h` — Decompose → BlockKernel → Reassemble
- CC_IBM path: uses original mesh-level `CorrParticleKernelTask` (CUTFACE_VELOCITIES requires full mesh)
- Verified: pending

---

### WallBCKernelTask (inside WallBC sub-graph)

**Graph location:** WallBC sub-graph kernel
**File:** `task/wallbc_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `WALL_BC_PROCESS_CELLS_KERNEL` | wall_kernels.f90:585 | **Mesh** | Wall loop (IW=1 to N_EXTERNAL+N_INTERNAL), iterates over wall cell index, not (I,J,K) grid |

**Task classification:** Mesh — wall cell iteration (IW index), cannot be restricted to K sub-ranges.
**Status:** DONE (classified, no conversion needed)

---

### CorrRadiationKernelTask (inside CorrRadiation sub-graph)

**Graph location:** CorrRadiation sub-graph kernel
**File:** `task/corr_radiation_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `COMPUTE_RADIATION_KERNEL` | radi.f90:3443 | **Mixed** | FVM angle sweeps (I,J,K loops per angle, decomposable) + wall loops + particle loops + global reductions (RAD_Q_SUM). 1300+ lines with complex CONTAINS structure. |

**Task classification:** Mixed — angle sweep cell loops could theoretically be K-decomposed, but wall/particle loops and complex 3D sweep pattern make it impractical. Would require major restructuring.
**Status:** DONE (classified, no conversion — complexity too high for benefit)

---

### CorrDivPart1KernelTask

**Graph location:** After MeshExchange(2)+InitDiv
**File:** `task/corr_div_part1_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_1_KERNEL` | divg_kernels.f90:27 | **Mesh** | Cell loops + multiple wall loops + per-mesh local accumulators (D_SUM_LOC, P_SUM_LOC, U_SUM_LOC) |

Note: `fds_combustion_bc_kernel` runs before this kernel in the same task.

**Task classification:** Mesh — same as predictor instance (#6).
**Status:** DONE (classified, no conversion needed)

---

### DivergencePart2KernelTask (corrector instance)

**Graph location:** After DivergenceExchange barrier
**File:** `task/divergence_part2_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_2_KERNEL` | divg_kernels.f90:1396 | **Mesh** | Same as predictor instance — zone loops + wall loops |

**Task classification:** Mesh — same as predictor instance.
**Status:** DONE (classified, no conversion needed)

---

### PressureSolveKernelTask (corrector instance)

**Graph location:** Inside PressureIteration sub-graph cycle
**File:** `task/pressure_iteration_tasks.h`

Same kernels as predictor instance (see above).

**Task classification:** Mesh — same as predictor instance.
**Status:** DONE (classified, no conversion needed)

---

### VelocityCorrectorKernelTask

**Graph location:** After PressureIteration
**File:** `task/velocity_corrector_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `VELOCITY_CORRECTOR_KERNEL` | velo_kernels.f90:159 | **Mesh Block** | Three pure I,J,K loops for U, V, W (staggered grid). No cross-cell deps beyond neighbor reads. |
| 2 | `CHECK_DIVERGENCE_KERNEL` | divg_kernels.f90:1608 | **Mesh** | I,J,K loop with max/min reductions to mesh scalars (RESMAX, DIVMX, DIVMN + index tracking) |

**Task classification:** Mixed — VELOCITY_CORRECTOR block-decomposed along K, CHECK_DIVERGENCE runs at mesh level after reassembly.
**Status:** DONE

**Implementation:**
- Block kernel: `VELOCITY_CORRECTOR_BLOCK_KERNEL(M, DT, K1, K2)` in velo_kernels.f90
- K-partition: cells [K1,K2] within [1,KBAR]; U/V at K=K1:K2, W at K=K1-1:K2-1 (+KBAR for last block)
- Sub-graph: `graph/velocity_corrector_block_subgraph.h` — Decompose → BlockKernel → Reassemble → CheckDiv
- CC_IBM path: uses original mesh-level `VelocityCorrectorKernelTask` (CC_PROJECT_VELOCITY requires full mesh)
- Verified: 12/12 custom pass, 45/59 verification pass (no regressions)

---

### VelocityBCEdgesTask (corrector instance, inside CorrFinal sub-graph)

**Graph location:** CorrFinal sub-graph kernel
**File:** `task/velocity_bc_edges_task.h`

Same kernels as predictor instance (see above).

**Task classification:** Mesh — same as predictor instance.
**Status:** DONE (classified, no conversion needed)

---

## Timestep Pipeline

### DumpMeshOutputsTask

**Graph location:** After TimestepGlobal barrier
**File:** `task/timestep_tasks.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `fds_dump_mesh_outputs` | dump.f90 | | I/O operations, file writes — not a computational kernel |

**Task classification:** Mesh (I/O, not parallelizable via blocks)
**Status:** Not applicable — I/O task, no block decomposition needed

---

## Summary

### Unique Kernel Tasks (deduplicated across predictor/corrector)

| # | Task | Kernels | Classification |
|---|------|---------|----------------|
| 1 | VelocityPredictorKernelTask | VELOCITY_PREDICTOR_KERNEL, CHECK_STABILITY_KERNEL | **DONE** — block kernel + mesh CheckStability |
| 2 | VelocityCorrectorKernelTask | VELOCITY_CORRECTOR_KERNEL, CHECK_DIVERGENCE_KERNEL | **DONE** — block kernel + mesh CheckDiv |
| 3 | DensityPredKernelTask | DENSITY_KERNEL | **Mesh** — zone loops, wall loops, CHECK_MASS_DENSITY |
| 4 | PredStep1/CorrStep1KernelTask | COMPUTE_VISCOSITY, MASS_FINITE_DIFFS, DENSITY | **Mesh** — all have wall loops |
| 5 | DivSetupKernelTask | VELOCITY_FLUX_KERNEL | **Mesh** — wall loops, CONTAINS host association |
| 6 | PredWallDivKernelTask | PARTICLE_MOMENTUM + DIVERGENCE_PART_1 | **Mixed** — PART_MOM block-able but lightweight; DIV_PART_1 mesh |
| 7 | DivergencePart2KernelTask | DIVERGENCE_PART_2_KERNEL | **Mesh** — zone loops, wall loops |
| 8 | CombustionKernelTask | COMBUSTION_GENERAL_KERNEL | **Mesh** — cell list, STOP_STATUS, CONTAINS host association |
| 9 | CorrCondensKernelTask | CONDENSATION_EVAPORATION_KERNEL | **Mesh** — wall+cell loops interleaved in species loop |
| 10 | ParticleMassEnergyKernelTask | PARTICLE_MASS_ENERGY_TRANSFER_KERNEL | **Mesh** — particle loop (DO IP=1,NLP) |
| 11 | CorrParticleKernelTask | PARTICLE_MOMENTUM_TRANSFER_KERNEL | **DONE** — block kernel, CC_IBM mesh fallback |
| 12 | WallBCKernelTask | WALL_BC_PROCESS_CELLS_KERNEL | **Mesh** — wall cell iteration (IW index) |
| 13 | CorrRadiationKernelTask | COMPUTE_RADIATION_KERNEL | **Mixed** — angle sweeps + wall/particle loops |
| 14 | PressureSolveKernelTask | NO_FLUX + COMPUTE_RHS + FFT/ULMAT + CHECK_RESIDUALS | **Mesh** — FFT/ULMAT global solvers, wall loops |
| 15 | VelocityBCEdgesTask | MATCH_VELOCITY, VELOCITY_BC_PROCESS_EDGES | **Mesh** — wall/edge loops |
| 16 | RetryMomentumDivKernelTask | PARTICLE_MOMENTUM + DIVERGENCE_PART_1 | **Mixed** — same as #6 |

### Progress Counters

- **Total unique kernel tasks:** 16
- **Classified:** 16 / 16
- **Converted to mesh block:** 3 (VelocityPredictor, VelocityCorrector, CorrParticleMomentum)
- **Confirmed mesh-only:** 11
- **Mixed (no conversion):** 2 (PredWallDiv, RetryMomentumDiv — PART_MOM lightweight relative to DIV_PART_1)
