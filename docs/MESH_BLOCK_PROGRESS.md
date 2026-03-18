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
| 1 | `COMPUTE_VISCOSITY_KERNEL` | velo_kernels.f90:1268 | **Mixed** | Cell loops (MU_DNS, STRAIN_RATE, turb MU, KRES) block-decomposable; wall loops + corner mirroring sequential |
| 2 | `MASS_FINITE_DIFFERENCES_NEW_KERNEL` | mass_kernels.f90:23 | **Mesh** | Cell loops (I,J,K) + wall face correction loops (WALL_LOOP_2, WALL_LOOP_3) |

**Task classification:** Mixed — viscosity cell loops block-decomposed for non-DEARDORFF/DYNSMAG/CC_IBM; wall loops sequential.
**Status:** DONE — `COMPUTE_VISCOSITY_BLOCK_KERNEL` + `COMPUTE_VISCOSITY_POST_BLOCK` (velo_kernels.f90). Block sub-graph: `graph/compute_viscosity_block_subgraph.h`. Supports NO_TURB, CONSMAG, VREMAN, WALE. Falls back to mesh-level for DEARDORFF (default), DYNSMAG, CC_IBM.

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
| 1 | `VELOCITY_FLUX_KERNEL` | velo_kernels.f90:325 | **Mesh Block (2-phase)** | Main cell loops (vorticity, FVX, FVY, FVZ, DIRECT_FORCE) are block-decomposable. Wall loop in CORIOLIS_FORCE (line 698) writes ghost cells only — excluded via runtime check. |

Note: `fds_set_baroclinic_false` and `fds_viscosity_bc_kernel` are called before this kernel but are sequential pre-processing (run in orchestrator before block dispatch).

**Task classification:** Mesh Block (conditional) — block-decomposed when no Coriolis/patch/CTRL/wind/periodic features; mesh-level fallback otherwise.
**Status:** DONE

**Implementation:**
- Block kernel: `VELOCITY_FLUX_BLOCK_KERNEL(M, T, DT, NM, APPLY_TO_ESTIMATED, K1, K2)` in velo_kernels.f90
- K-partition: vorticity at K=max(0,K1-1):K2 (extended for FVX/FVY dependency); FVX/FVY at K1:K2; FVZ at K1-1:K2-1 (staggered, KBAR for last block)
- Sub-graph: `graph/velocity_flux_block_subgraph.h` — Orchestrator(pre-proc + decompose) → BlockKernel → Collector(agglomeration)
- Runtime check: `fds_velocity_flux_can_block_decompose()` excludes CC_IBM, CTRL_DIRECT_FORCE, Coriolis, patch velocity, open wind, periodic tests
- CC_IBM path: uses original mesh-level `DivSetupKernelTask`
- Verified: 12/12 custom pass, 45/59 verification pass (no regressions)

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
**File:** `task/divergence_part2_kernel_task.h`, `graph/divergence_part2_block_subgraph.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_2_KERNEL` | divg_kernels.f90:1396 | **Block (2-phase)** | Zone ops sequential in orchestrator; cell loops + BC_LOOP block-decomposed along K |

**Task classification:** Block (conditional) — zone ops (R_PBAR, USUM, D_PBAR_DT_P) run sequentially per-mesh in orchestrator; pressure zone DP, solid zeroing, BC_LOOP, DIV+DDDT decomposed into K-blocks. CC_IBM falls back to mesh-level.
**Status:** DONE

**Implementation:**
- Preprocessing: `DIVERGENCE_PART_2_PREPROCESSING(M, DT, NM)` in divg_kernels.f90 — R_PBAR computation, zone ops (USUM_ADD, USUM modification), D_PBAR_DT_P computation, CC_IBM GET_LINKED_VELOCITIES
- Block kernel: `DIVERGENCE_PART_2_BLOCK_KERNEL(M, DT, NM, K1, K2)` in divg_kernels.f90 — pressure zone DP (K1:K2), solid cell zeroing (K1:K2), BC_LOOP (K-filtered), DIV+DDDT computation (K1:K2)
- K-partition: BC_LOOP wall cells filtered by BC%KK coordinate; first block includes K=0 ghost cells via `K_MIN = MERGE(0, K1, K1==1)`
- Sub-graph: `graph/divergence_part2_block_subgraph.h` — Orchestrator(preprocessing + decompose) → BlockKernel(parallel) → Reassemble
- CC_IBM path: uses original mesh-level `DivergencePart2KernelTask` (GET_LINKED_VELOCITIES/GET_CUTCELL_DDDT require full mesh)
- Verified: 8/12 custom pass (same 4 pre-existing failures), 42/59 verification pass (no regressions)

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
**File:** `task/velocity_bc_edges_task.h`, `graph/velocity_bc_edges_block_subgraph.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `MATCH_VELOCITY_KERNEL` | velo_kernels.f90:2238 | **Mesh** | External wall loop — sequential preprocessing in orchestrator |
| 2 | `VELOCITY_BC_PROCESS_EDGES_KERNEL` | velo_kernels.f90:2361 | **Block** | Edge loop filtered by ED%K range; DRAG_UVWMAX via local accumulator + MAX reduction |

**Task classification:** Block — edge loop filtered by K coordinate. MATCH_VELOCITY + VELOCITY_BC_PREPROCESSING run in orchestrator; DRAG_UVWMAX reduced in collector.
**Status:** DONE

**Implementation:**
- Block kernel: `VELOCITY_BC_PROCESS_EDGES_KERNEL(M,NM,T,...,K1_IN,K2_IN,DRAG_UVWMAX_LOCAL)` — OPTIONAL parameters in velo_kernels.f90
- K-partition: edges filtered by ED%K; first block includes K=0 boundary edges
- DRAG_UVWMAX: per-block local accumulator, MAX-reduced in collector, written back via `fds_set_drag_uvwmax`
- Sub-graph: `graph/velocity_bc_edges_block_subgraph.h` — Orchestrator(SYNTHETIC_TURBULENCE+MATCH+PREPROCESS+decompose) → BlockKernel → Collector(reassemble+DRAG reduce+CC_VELOCITY_BC)
- CC_VELOCITY_BC runs in collector (mesh-level)
- Verified: 12/12 custom pass, 45/59 verification pass (no regressions), WUI vegetation tests pass

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
| 1 | `VELOCITY_FLUX_KERNEL` | velo_kernels.f90:325 | **Mesh Block (2-phase)** | Same as predictor — block-decomposed with runtime feature check |

Note: corrector also calls `fds_agglomeration(dt, nm)` after the block kernel in the collector — sequential post-processing.

**Task classification:** Mesh Block (conditional) — same as predictor instance.
**Status:** DONE (same implementation as predictor)

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
**File:** `task/wallbc_kernel_task.h` (mesh-level), `graph/wallbc_block_subgraph.h` (block-level)

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `WALL_BC_PROCESS_CELLS_KERNEL` | wall.f90:1274 | **Block** | Wall cell loop filtered by KKG range (Approach A: Wall Cell K-Indexing) |

**Task classification:** Block — wall cells filtered by gas cell K-coordinate (BC%KKG) into K sub-ranges. CFACE cells skipped (CC_IBM excluded). Particle loop also filtered by KKG.
**Status:** DONE — `WALL_BC_PROCESS_CELLS_BLOCK_KERNEL(M,NM,...,K1,K2)` in wall.f90. Block sub-graph: `graph/wallbc_block_subgraph.h`. Orchestrator runs preprocessing per-mesh, K-decomposes, collector runs WALL_BC_FINALIZE. Falls back to mesh-level for CC_IBM.

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
**File:** `task/divergence_part2_kernel_task.h`, `graph/divergence_part2_block_subgraph.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_2_KERNEL` | divg_kernels.f90:1396 | **Block (2-phase)** | Same as predictor instance — block-decomposed with orchestrator preprocessing |

**Task classification:** Block (conditional) — same as predictor instance. Shares block kernel implementation.
**Status:** DONE

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
**File:** `task/velocity_bc_edges_task.h`, `graph/velocity_bc_edges_block_subgraph.h`

Same kernels as predictor instance (see above).

**Task classification:** Block — same as predictor instance. Shares block kernel implementation.
**Status:** DONE

**Implementation:** Same as predictor instance. CorrFinal block sub-graph additionally runs `UPDATE_GLOBAL_OUTPUTS` in collector (per-mesh output accumulation). See predictor VelocityBCEdgesTask entry for full details.

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
| 5 | DivSetupKernelTask | VELOCITY_FLUX_KERNEL | **DONE** — block kernel (conditional), mesh fallback for Coriolis/patch/CTRL/wind/periodic |
| 6 | PredWallDivKernelTask | PARTICLE_MOMENTUM + DIVERGENCE_PART_1 | **Mixed** — PART_MOM block-able but lightweight; DIV_PART_1 mesh |
| 7 | DivergencePart2KernelTask | DIVERGENCE_PART_2_KERNEL | **DONE** — block kernel (zone ops sequential), CC_IBM mesh fallback |
| 8 | CombustionKernelTask | COMBUSTION_GENERAL_KERNEL | **Mesh** — cell list, STOP_STATUS, CONTAINS host association |
| 9 | CorrCondensKernelTask | CONDENSATION_EVAPORATION_KERNEL | **Mesh** — wall+cell loops interleaved in species loop |
| 10 | ParticleMassEnergyKernelTask | PARTICLE_MASS_ENERGY_TRANSFER_KERNEL | **Mesh** — particle loop (DO IP=1,NLP) |
| 11 | CorrParticleKernelTask | PARTICLE_MOMENTUM_TRANSFER_KERNEL | **DONE** — block kernel, CC_IBM mesh fallback |
| 12 | WallBCKernelTask | WALL_BC_PROCESS_CELLS_KERNEL | **DONE** — block kernel (KKG filter), CC_IBM mesh fallback |
| 13 | CorrRadiationKernelTask | COMPUTE_RADIATION_KERNEL | **Mixed** — angle sweeps + wall/particle loops |
| 14 | PressureSolveKernelTask | NO_FLUX + COMPUTE_RHS + FFT/ULMAT + CHECK_RESIDUALS | **Mesh** — FFT/ULMAT global solvers, wall loops |
| 15 | VelocityBCEdgesTask | MATCH_VELOCITY, VELOCITY_BC_PROCESS_EDGES | **DONE** — block kernel (ED%K filter), DRAG_UVWMAX MAX reduction |
| 16 | RetryMomentumDivKernelTask | PARTICLE_MOMENTUM + DIVERGENCE_PART_1 | **Mixed** — same as #6 |

### Progress Counters

- **Total unique kernel tasks:** 16
- **Classified:** 16 / 16
- **Converted to mesh block:** 8 (VelocityPredictor, VelocityCorrector, CorrParticleMomentum, DivSetup/VelocityFlux, ComputeViscosity, WallBC, VelocityBCEdges, DivergencePart2)
- **Confirmed mesh-only:** 6
- **Mixed (no conversion):** 2 (PredWallDiv, RetryMomentumDiv — PART_MOM lightweight relative to DIV_PART_1)

### Phase 4: Remaining Block Decomposition Targets

Ranked by estimated impact (runtime × feasibility):

| Priority | Kernel | Combined Time (ms) | Decomposable % | Feasibility | Blocker |
|----------|--------|-------------------|----------------|-------------|---------|
| ~~1~~ | ~~DIVERGENCE_PART_1_KERNEL~~ | ~~824~~ | ~~15-25%~~ | ~~LOW~~ | GET_SCALAR_FACE_VALUE stencils require full K-domain; wall+cell loops interleaved within species loops; wall corrections create race conditions at block boundaries |
| ~~2~~ | **DIVERGENCE_PART_2_KERNEL** | 140 (pred + corr) | 98% | **DONE** | Zone ops sequential in orchestrator; cell loops + BC_LOOP K-decomposed |
| 3 | DENSITY_KERNEL | ~150 | 70% | MODERATE | CHECK_MASS_DENSITY needs K±1 halo; zone PBAR updates sequential |
| 4 | MASS_FINITE_DIFFERENCES | ~100 | LOW | LOW | GET_SCALAR_FACE_VALUE stencils require full K-domain neighbor access |
| 5 | COMPUTE_RADIATION_KERNEL | ~250 | ~30% | LOW | Complex FVM solver (1300+ lines, CONTAINS, angle sweeps, wall+particle loops) |
| 6 | DumpMeshOutputs | — | 0% | NONE | I/O task, not parallelizable |
