# Phase 3: Easy Parallelization Targets

Three sequential nodes in the Hedgehog graph contain per-mesh loops that could be
parallelized with moderate effort. This document tracks the plan and progress.

## Overview

| Target | Current Location | Per-Mesh Cost | Mesh Count Benefit | Priority |
|--------|-----------------|---------------|-------------------|----------|
| Combustion (serial chemistry) | CombustionHvacTask | Moderate-High (ODE solve) | High (fire cases) | 1 |
| Particle Mass/Energy | CorrParticleOrchestrator | Moderate (per-particle HT) | Medium (particle cases) | 2 |
| Particle Insertion | PredStep1Orchestrator | Low-Moderate | Low (particle cases) | 3 |

All three follow the same pattern: a per-mesh loop in an orchestrator or barrier
task that currently runs sequentially but operates on independent mesh data.

---

## Target 1: Combustion Kernel Extraction

**Current code**: `CombustionHvacTask` in `barrier_tasks.h` calls `fds_combustion(t,dt)`.

**Fortran routine**: `COMBUSTION_LOAD_BALANCED` in `fire.f90:38-89`, which calls
`COMBUSTION_GENERAL_LOAD_BALANCED` in `fire.f90:92-340`.

**Structure to parallelize** (serial chemistry path, `fire.f90:243-304`):
```
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   POINT_TO_MESH(NM)
   DO NC=1,NCHEM_ACTIVE_CELLS  ! per-cell ODE solve
      CALL COMBUSTION_MODEL(...)
   ENDDO
ENDDO
```

**Why it matters**: COMBUSTION_MODEL is an expensive ODE solver. For fire
simulations with many chemically active cells, this loop can take significant time.

**HVAC constraint**: `HVAC_CALC` must remain sequential (global network solver,
zero per-mesh work). The current merged CombustionHvacTask must be split.

### Steps

- [x] **Step 0: Verify the plan against current code**
  - Read `fire.f90` COMBUSTION_LOAD_BALANCED and COMBUSTION_GENERAL_LOAD_BALANCED fully
  - Read `fire.f90` COMBUSTION_MODEL and its callees
  - Confirm the serial chemistry path (`fire.f90:243-304`) is truly per-mesh independent
  - Check what POINT_TO_MESH sets and which module-level pointers COMBUSTION_MODEL reads
  - Identify all global arrays written by the chemistry loop (Q, CHI_R, ZZ, etc.)
  - Verify the MPI load-balanced path is NOT taken in single-process mode
  - Check the cut-cell volume averaging loop (`fire.f90:307-331`) for independence
  - Check SOOT_SURFACE_OXIDATION for thread safety
  - Document findings: list all module-level variables accessed, note any race conditions

- [x] **Step 1: Create thread-safe combustion kernel**
  - Extract `COMBUSTION_KERNEL(M, NM, T, DT)` in a new `fire_kernels.f90` module
  - Replace POINT_TO_MESH usage with `M => MESHES(NM)` (Pattern 2)
  - Thread-safe conversion: local pointers for Q, CHI_R, ZZ, RHOS, etc.
  - Handle the three phases:
    1. Pre-chemistry: zero Q/CHI_R, CC_IBM cut-cell prep (per-mesh)
    2. Chemistry: identify active cells + COMBUSTION_MODEL loop (per-mesh)
    3. Post-chemistry: CC_IBM volume averaging + SOOT_SURFACE_OXIDATION (per-mesh)
  - Add C binding wrapper: `fds_combustion_kernel(nm, t, dt)`

- [x] **Step 2: Split CombustionHvacTask into Combustion sub-graph + HVAC**
  - Create `CombustionOrchestrator` (collects N meshes, runs global pre-processing if any)
  - Create `CombustionKernelTask` (parallel, calls `fds_combustion_kernel`)
  - Create `CombustionCollector` (collects N results, no post-processing needed)
  - Keep `HvacTask` as separate sequential barrier after the combustion sub-graph
  - Wire in corrector_subgraph.h: `DivSetup -> CombustionSubgraph -> HvacCollector -> HvacTask -> CorrCondens`

- [x] **Step 3: Verify correctness**
  - 12/12 custom tests pass, 45/59 verification pass (no regressions)

### Verification checklist (all verified)
- [x] Global arrays: Q_DOT accumulation — diagnostic only, already ignored in comparison
- [x] SOOT_SURFACE_OXIDATION: separated into SootHvacTask barrier (sequential)
- [x] MPI load-balanced path: bypassed in single-process mode
- [x] CHI_R array: per-mesh, thread-safe through kernel
- [x] CC_IBM cut-cell volume averaging: per-mesh, handled in kernel
- [x] Combustion pre-processing: kernel zeroes Q/CHI_R before chemistry for each mesh

---

## Target 2: Particle Mass/Energy Kernel Extraction

**Current code**: `CorrParticleOrchestrator` in `corr_particle_state.h` calls
`fds_particle_mass_energy(t, dt, nm)` + `fds_move_particles(t, dt, nm)` sequentially
for each mesh, then dispatches parallel `PARTICLE_MOMENTUM_KERNEL`.

**Key insight**: PARTICLE_MASS_ENERGY_TRANSFER is per-mesh independent (no cross-mesh
writes). MOVE_PARTICLES has cross-mesh particle transfer (writes to OMESH send buffers)
and must stay sequential.

**Proposed split**:
```
Current:  [sequential: MASS_ENERGY + MOVE for all meshes] -> [parallel: MOMENTUM]
Proposed: [parallel: MASS_ENERGY] -> [sequential: MOVE for all meshes] -> [parallel: MOMENTUM]
```

### Steps

- [x] **Step 0: Verify the plan against current code**
  - PARTICLE_MASS_ENERGY_TRANSFER (~1000 lines) only operates on MESHES(NM) data
  - REMOVE_PARTICLES has cross-mesh writes (OMESH send buffers) — kept sequential
  - Q_DOT/M_DOT global accumulators — diagnostic only, already ignored in comparison
  - All BOUNDARY_PROP1/PROP2 arrays are mesh-local (indexed by B1_INDEX/B2_INDEX)

- [x] **Step 1: Create thread-safe mass/energy kernel**
  - Extracted `PARTICLE_MASS_ENERGY_KERNEL(NM, T, DT)` in `part.f90` (same module)
  - Used Pattern 3 (local alias shadowing) — ~35 local aliases shadow MESH_POINTERS
  - Added C binding: `fds_particle_mass_energy_kernel(nm, t, dt)`
  - Added `fds_remove_particles(t, nm)` binding (calls POINT_TO_MESH + REMOVE_PARTICLES)

- [x] **Step 2: Restructure CorrParticleOrchestrator**
  - Replaced CorrParticleOrchestrator with direct wiring:
    1. ParticleMassEnergyKernelTask (parallel) -> Collector -> RemoveMoveParticlesTask (barrier)
    2. RemoveMoveParticlesTask -> CorrParticleKernelTask (parallel MOMENTUM)
  - Eliminated CorrParticleOrchestrator entirely

- [x] **Step 3: Verify correctness**
  - 12/12 custom tests pass (including bucket_test_1_short, activate_sprinklers)
  - 45/59 verification pass (no regressions, same as before)
  - All Sprinklers_and_Sprays/ cases pass (terminal_velocity, flat_fire)

### Verification checklist (all verified)
- [x] Thin wall heat transfer: handled by REMOVE_PARTICLES in sequential barrier
- [x] Particle-gas coupling arrays: M_DOT_PPP, D_SOURCE per-mesh — thread-safe
- [x] LAGRANGIAN_PARTICLE: per-mesh arrays, no cross-mesh manipulation in kernel
- [x] Wall cell updates: LP_CPUA, LP_MPUA in BOUNDARY_PROP2 — per-mesh, thread-safe
- [x] Species source terms: SPECIES_MIXTURE global, read-only — safe

---

## Target 3: Particle Insertion Kernel Extraction

**Current code**: `PredStep1Orchestrator` in `pred_step1_state.h` calls
`fds_insert_particles(t, nm)` sequentially for each mesh, then dispatches parallel
`COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES` kernels.

**Key insight**: INSERT_ALL_PARTICLES operates only on mesh NM. No cross-mesh data
writes — only sets a global `EXCHANGE_INSERTED_PARTICLES` flag afterward.

**Proposed change**: Merge INSERT_ALL_PARTICLES into the PredStep1KernelTask so it
runs in parallel with COMPUTE_VISCOSITY and MASS_FINITE_DIFFERENCES.

### Status: BLOCKED — Not viable for parallelization

**Analysis completed**: INSERT_ALL_PARTICLES has fundamental thread-safety blockers:

1. **RANDOM_NUMBER** — Fortran's intrinsic PRNG is NOT thread-safe. Used extensively
   for spray angles, diameters, particle placement (INSERT_SPRAY_PARTICLES,
   INSERT_VENT_PARTICLES, INSERT_VOLUMETRIC_PARTICLES). Would need per-mesh PRNG
   implementation.
2. **Global state modifications** — EXCHANGE_INSERTED_PARTICLES flag,
   N_ACTUATED_SPRINKLERS counter, INITIALIZATION%ALREADY_INSERTED.
3. **DEVICE_VARIABLES** — Sprinkler activation reads/writes global device state.
4. **File I/O** — Reads binary bulk density files.
5. **Cross-mesh OMESH access** — Reads OMESH send buffer counts.

Parallelizing this would require either:
- A thread-safe PRNG per mesh (significant refactoring)
- Or accepting non-reproducible particle placement

**Recommendation**: Skip this target. The cost/benefit ratio is poor — particle
insertion is a small fraction of total runtime, and the refactoring is high-risk.

---

## Test Coverage Requirements

The quick test suite (`test_cases/run_tests.py`) must include cases that exercise:

| Feature | Test Case | Meshes | Status |
|---------|-----------|--------|--------|
| Combustion (REAC) | multiple_reac_3mesh | 3 | Existing |
| Combustion (multi-mesh fire) | circular_burner_short | 8 | **New** |
| Particles (sprinkler) | bucket_test_1_short | 4 | **New** |
| Particles (activation/control) | activate_sprinklers | 1 | **New** |
| ULMAT pressure solver | dancing_eddies_ulmat | 4 | **New** |
| Non-reacting multi-mesh | dancing_eddies_4mesh_short | 4 | Existing |
| Non-reacting single-mesh | dancing_eddies_1mesh_short | 1 | Existing |

---

## Execution Order

Recommended order of implementation:

1. **Combustion** (highest payoff for fire cases, well-isolated from HVAC)
2. **Particle Mass/Energy** (moderate payoff, requires sub-graph restructuring)
3. **Particle Insertion** (lowest payoff, simplest change but RANDOM_NUMBER risk)

Each target should be implemented, tested, and committed independently before
starting the next one.
