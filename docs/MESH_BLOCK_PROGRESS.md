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
| 1 | `COMPUTE_VISCOSITY_KERNEL` | velo_kernels.f90:807 | | CONTAINS subroutines (WALE, Deardorff, etc.), wall loops, turbulence model dispatch |
| 2 | `MASS_FINITE_DIFFERENCES_NEW_KERNEL` | mass_kernels.f90:28 | | Cell loops (I,J,K) for scalar transport + wall face correction loop |

**Task classification:** Pending
**Status:** Not started

---

### DensityPredKernelTask

**Graph location:** After PredStep1KernelTask
**File:** `task/density_pred_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DENSITY_KERNEL` | mass_kernels.f90:358 | | Multiple I,J,K loops (density update, CHECK_MASS_DENSITY redistribution) |

**Task classification:** Pending
**Status:** Not started

---

### DivSetupKernelTask (predictor instance)

**Graph location:** After MeshExchange(1), optional CC_IBM orchestrator
**File:** `task/div_setup_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `VELOCITY_FLUX_KERNEL` | velo_kernels.f90:211 | | Large kernel with CONTAINS subroutines, wall loops, species loops, cell loops |

Note: `fds_set_baroclinic_false` and `fds_viscosity_bc_kernel` are called before this kernel but are sequential pre-processing (run in orchestrator or before dispatch).

**Task classification:** Pending
**Status:** Not started

---

### PredWallDivKernelTask

**Graph location:** After WallBC sub-graph
**File:** `task/pred_wall_div_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MOMENTUM_TRANSFER_KERNEL` | part_kernels.f90:27 | | Single I,J,K loop (0:IBAR, 0:JBAR, 0:KBAR), pure stencil |
| 2 | `DIVERGENCE_PART_1_KERNEL` | divg_kernels.f90:43 | | Cell loops + wall loops + per-mesh local accumulators (D_SUM_LOC, P_SUM_LOC, U_SUM_LOC) |

**Task classification:** Pending
**Status:** Not started

---

### DivergencePart2KernelTask (predictor instance)

**Graph location:** After DivergenceExchange barrier
**File:** `task/divergence_part2_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_2_KERNEL` | divg_kernels.f90:1404 | | Pressure update I,J,K loop + zone loop for pressure averaging |

**Task classification:** Pending
**Status:** Not started

---

### PressureSolveKernelTask (predictor instance)

**Graph location:** Inside PressureIteration sub-graph cycle
**File:** `task/pressure_iteration_tasks.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `NO_FLUX_KERNEL` | pres_kernels.f90:737 | | Wall loop (IW=1 to N_EXTERNAL+N_INTERNAL) |
| 2 | `PRESSURE_SOLVER_COMPUTE_RHS` | pres_kernels.f90:24 | | Wall pre-processing + I,J,K RHS assembly loop |
| 3 | `PRESSURE_SOLVER_FFT` | pres_kernels.f90:308 | | FFT solver (tridiagonal + spectral transform) — global mesh operation |
| 4 | `PRESSURE_SOLVER_CHECK_RESIDUALS` | pres_kernels.f90:477 | | I,J,K residual loop + baroclinic correction with wall checks |

Or ULMAT variant:
| 3a | `ULMAT_SOLVER_KERNEL` | pres.f90 | | PARDISO/HYPRE direct solve per zone — global mesh operation |
| 4a | `PRESSURE_SOLVER_CHECK_RESIDUALS_U_KERNEL` | pres_kernels.f90:586 | | I,J,K loops + inline GRADIENT_WEIGHT |

**Task classification:** Pending
**Status:** Not started

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
| 1 | `MATCH_VELOCITY_KERNEL` | velo_kernels.f90:2238 | | External wall loop (IW=1 to N_EXTERNAL_WALL_CELLS) |
| 2 | `VELOCITY_BC_PROCESS_EDGES_KERNEL` | velo_kernels.f90:1430 | | Edge loop with wall BC application |

Note: `fds_velocity_bc_preprocessing` runs before the kernel but is sequential pre-processing.

**Task classification:** Pending
**Status:** Not started

---

### RetryMomentumDivKernelTask (inside ChangeTimeStep sub-graph)

**Graph location:** ChangeTimeStep CFL retry loop
**File:** `task/change_timestep_tasks.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MOMENTUM_TRANSFER_KERNEL` | part_kernels.f90:27 | | Same as PredWallDiv — pure I,J,K loop |
| 2 | `DIVERGENCE_PART_1_KERNEL` | divg_kernels.f90:43 | | Same as PredWallDiv — cell+wall loops + local accumulators |

**Task classification:** Pending
**Status:** Not started

---

## Corrector Sub-Graph

### CorrStep1KernelTask

**Graph location:** Corrector entry
**File:** `task/corr_step1_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `COMPUTE_VISCOSITY_KERNEL` | velo_kernels.f90:807 | | Same as PredStep1 — CONTAINS, wall loops, turbulence dispatch |
| 2 | `MASS_FINITE_DIFFERENCES_NEW_KERNEL` | mass_kernels.f90:28 | | Same as PredStep1 — cell loops + wall correction |
| 3 | `DENSITY_KERNEL` | mass_kernels.f90:358 | | Same as DensityPred — cell loops + CHECK_MASS_DENSITY |

**Task classification:** Pending
**Status:** Not started

---

### DivSetupKernelTask (corrector instance)

**Graph location:** After MeshExchange(4), optional CC_IBM orchestrator
**File:** `task/div_setup_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `VELOCITY_FLUX_KERNEL` | velo_kernels.f90:211 | | Same as predictor instance |

Note: corrector also calls `fds_agglomeration(dt, nm)` after the kernel — this is a sequential operation.

**Task classification:** Pending
**Status:** Not started

---

### CombustionKernelTask

**Graph location:** After DivSetup, before SootHvac barrier
**File:** `task/combustion_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `COMBUSTION_MODEL_KERNEL` | fire_kernels.f90:1143 | | Cell loops for combustion + species reaction loops per cell |

**Task classification:** Pending
**Status:** Not started

---

### CorrCondensKernelTask

**Graph location:** After SootHvac barrier
**File:** `task/corr_condens_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `CONDENSATION_EVAPORATION_KERNEL` | fire_kernels.f90:1539 | | Cell loop with local species chemistry per cell |

**Task classification:** Pending
**Status:** Not started

---

### ParticleMassEnergyKernelTask

**Graph location:** After CorrCondens, before RemoveMoveParticles barrier
**File:** `task/particle_mass_energy_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MASS_ENERGY_TRANSFER_KERNEL` | part_kernels.f90 (wrapper) | | Particle loop (DO IP=1,NLP) — iterates over Lagrangian particles |

**Task classification:** Pending
**Status:** Not started

---

### CorrParticleKernelTask

**Graph location:** After RemoveMoveParticles barrier
**File:** `task/corr_particle_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `PARTICLE_MOMENTUM_TRANSFER_KERNEL` | part_kernels.f90:27 | | Pure I,J,K loop (0:IBAR, 0:JBAR, 0:KBAR) |

**Task classification:** Pending
**Status:** Not started

---

### WallBCKernelTask (inside WallBC sub-graph)

**Graph location:** WallBC sub-graph kernel
**File:** `task/wallbc_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `WALL_BC_PROCESS_CELLS_KERNEL` | wall_kernels.f90:585 | | Wall loop (IW=1 to N_EXTERNAL+N_INTERNAL), per-wall-cell operations |

**Task classification:** Pending
**Status:** Not started

---

### CorrRadiationKernelTask (inside CorrRadiation sub-graph)

**Graph location:** CorrRadiation sub-graph kernel
**File:** `task/corr_radiation_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `COMPUTE_RADIATION_KERNEL` | radi.f90 | | FVM radiation solver with angle loops + cell loops, complex CONTAINS structure |

**Task classification:** Pending
**Status:** Not started

---

### CorrDivPart1KernelTask

**Graph location:** After MeshExchange(2)+InitDiv
**File:** `task/corr_div_part1_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_1_KERNEL` | divg_kernels.f90:43 | | Cell loops + wall loops + per-mesh local accumulators |

Note: `fds_combustion_bc_kernel` runs before this kernel in the same task.

**Task classification:** Pending
**Status:** Not started

---

### DivergencePart2KernelTask (corrector instance)

**Graph location:** After DivergenceExchange barrier
**File:** `task/divergence_part2_kernel_task.h`

| # | Fortran Kernel | Source | Classification | Notes |
|---|----------------|--------|----------------|-------|
| 1 | `DIVERGENCE_PART_2_KERNEL` | divg_kernels.f90:1404 | | Same as predictor instance |

**Task classification:** Pending
**Status:** Not started

---

### PressureSolveKernelTask (corrector instance)

**Graph location:** Inside PressureIteration sub-graph cycle
**File:** `task/pressure_iteration_tasks.h`

Same kernels as predictor instance (see above).

**Task classification:** Pending
**Status:** Not started

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

**Task classification:** Pending
**Status:** Not started

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

| # | Task | Kernels | Expected Classification |
|---|------|---------|------------------------|
| 1 | VelocityPredictorKernelTask | VELOCITY_PREDICTOR_KERNEL, CHECK_STABILITY_KERNEL | **DONE** — block kernel + mesh CheckStability |
| 2 | VelocityCorrectorKernelTask | VELOCITY_CORRECTOR_KERNEL, CHECK_DIVERGENCE_KERNEL | **DONE** — block kernel + mesh CheckDiv |
| 3 | DensityPredKernelTask | DENSITY_KERNEL | Likely mixed (I,J,K + CHECK_MASS_DENSITY) |
| 4 | PredStep1/CorrStep1KernelTask | COMPUTE_VISCOSITY_KERNEL, MASS_FINITE_DIFFERENCES_NEW_KERNEL, DENSITY_KERNEL | Likely mixed (wall loops in viscosity) |
| 5 | DivSetupKernelTask | VELOCITY_FLUX_KERNEL | Likely mixed (large, CONTAINS, wall loops) |
| 6 | PredWallDivKernelTask | PARTICLE_MOMENTUM_TRANSFER_KERNEL, DIVERGENCE_PART_1_KERNEL | Likely mixed (block + wall + accumulators) |
| 7 | DivergencePart2KernelTask | DIVERGENCE_PART_2_KERNEL | Likely mixed (I,J,K + zone loop) |
| 8 | CombustionKernelTask | COMBUSTION_MODEL_KERNEL | Pending analysis |
| 9 | CorrCondensKernelTask | CONDENSATION_EVAPORATION_KERNEL | Likely mesh block (per-cell chemistry) |
| 10 | ParticleMassEnergyKernelTask | PARTICLE_MASS_ENERGY_TRANSFER_KERNEL | Likely mesh (particle loop) |
| 11 | CorrParticleKernelTask | PARTICLE_MOMENTUM_TRANSFER_KERNEL | Likely mesh block (pure I,J,K) |
| 12 | WallBCKernelTask | WALL_BC_PROCESS_CELLS_KERNEL | Likely mesh (wall loop) |
| 13 | CorrRadiationKernelTask | COMPUTE_RADIATION_KERNEL | Likely mesh (FVM angle sweep, complex) |
| 14 | PressureSolveKernelTask | NO_FLUX + COMPUTE_RHS + FFT/ULMAT + CHECK_RESIDUALS | Likely mesh (FFT/ULMAT are global solvers) |
| 15 | VelocityBCEdgesTask | MATCH_VELOCITY_KERNEL, VELOCITY_BC_PROCESS_EDGES_KERNEL | Likely mesh (wall/edge loops) |
| 16 | RetryMomentumDivKernelTask | PARTICLE_MOMENTUM_TRANSFER_KERNEL, DIVERGENCE_PART_1_KERNEL | Same as PredWallDiv (#6) |

### Progress Counters

- **Total unique kernel tasks:** 16
- **Classified:** 2 / 16
- **Converted to mesh block:** 2
- **Confirmed mesh-only:** 0
