# Parallelization Progress

## Current State

**Branch**: `hedgehog-integration`
**Tests**: 20/20 custom, 58/58 verification (tol=1e-6)

### Profile (dancing_eddies_4mesh, 27 steps, total 3.739s)

| Category | Time | % |
|----------|------|---|
| Parallel kernels | 2301 ms | 61.5% |
| Sequential barriers | 885 ms | 23.7% |
| Overhead (dispatch/wake) | 553 ms | 14.8% |

Sequential breakdown:

| Barrier | ms | Contents |
|---------|----|----------|
| TimestepState | 334 | INSERT_PARTICLES (per-mesh loop), STOP_CHECK, ADJUST_DT |
| PressureIteration (pred) | 148 | Parallel when pressure subgraph enabled |
| ChangeTimeStep | 133 | DIV_P2 + VEL_PRED per-mesh loops (sequential) |
| PressureIteration (corr) | 125 | Parallel when pressure subgraph enabled |
| MeshExchanges | 70 | MPI communication |
| CorrFinalCollector | 38 | HRR/mass reduce |
| Other | 37 | PhaseTransition, InitDiv |

### MPI vs Hedgehog Gap

MPI (1 process/mesh) is faster because:
- **~5-10 sync points/step** vs HH's **~15-20 barriers/step**
- Between MPI exchanges, per-mesh loops run independently with zero synchronization
- HH breaks these into kernel->barrier->kernel chains

Worst case (pre-Phase 6): Corrector MeshExch(4)->MeshExch(7) had **0 MPI barriers** but **4 HH barriers** (Fork1Collector, SootOxidation, HVAC, SootHvacJoin). Phase 6 eliminated 3 of these.

## Phase 6: Barrier Elimination (COMPLETE)

### Target A: SOOT_OXIDATION -> per-mesh kernel ✓

**Status**: DONE

Merged soot oxidation into `Fork1CombKernelTask` (per-mesh, thread-safe).
Eliminated 3 barriers: Fork1Collector, SootChainTask, SootHvacJoin.
Corrector flow: `DivSetup || (Combustion+Soot) -> ForkJoin -> Fork(HVAC || ParticleOps) -> Join`

### Target B: INSERT_ALL_PARTICLES -> PredStep1Kernel ✓

**Status**: DONE

Moved INSERT_ALL_PARTICLES into `PredStep1KernelTask` (parallel per-mesh).
EXCHANGE_INSERTED_PARTICLES handled in MeshExch1 barrier.
Reduced TimestepState sequential time.

### Target C: HVAC || ParticleOps fork ✓

**Status**: DONE

HVAC_CALC (global network solve) and ParticleOps (per-mesh kernel) run in parallel.
New corrector flow: `ForkJoin -> Fork(HVAC || ParticleOps) -> HvacPartJoin -> MeshExch7 -> WallBC`

### Target D: RetryLoopState restructure -> parallel DivP2 + VelPred ✓

**Status**: DONE

Replaced monolithic `RetryLoopState` (sequential per-mesh loops) with pipeline:
```
RetryMomDivKernel → RetryDivExchSM → RetryDivP2Kernel → RetryPressureSM → RetryVelPredKernel → RetryCheckSM
```

Key changes:
- `RetryCheckState` replaces `RetryLoopState` — only does stop_check + retry decision
- DivP2: barrier (exchange + preprocessing) → parallel block kernel
- VelPred: parallel kernel (kernel_only → cc_project → wall_velocity → check_stability)
- Reuses existing `DivergencePart2KernelTask` and `VelocityPredictorKernelTask`

Bugs found during implementation:
1. MATCH_VELOCITY incorrectly added — VELOCITY_PREDICTOR orchestrator does NOT call it
2. CHECK_STABILITY ordering: must come AFTER CC_PROJECT + WALL_VELOCITY_NO_GRADH (ULMAT sensitivity)

### Completed Impact

| Target | Change | Effect |
|--------|--------|--------|
| A: Soot -> kernel | 3 barriers eliminated | Combustion+soot fully parallel |
| B: InsertPart -> kernel | TimestepState reduced | INSERT parallel per-mesh |
| C: HVAC \|\| ParticleOps | Fork added | ParticleOps overlaps HVAC |
| D: Retry parallel | Pipeline replaces monolith | DivP2+VelPred parallel per-mesh |

## Completed Phases (1-5)

### Phase 1-2: Kernel extraction + CC_IBM thread safety
- 9 CC_DENSITY_TS routines, WALL_BC_FINALIZE extraction
- POINT_TO_MESH removed from all parallel tasks
- T_USED/T_CC_USED timing races eliminated
- 15 sub-graphs, CC_IBM fully thread-safe

### Phase 3: Combustion + particle mass/energy parallelization
- Fork1: DivSetup || Combustion (parallel kernels)
- ParticleOps: CONDENSATION + PARTICLE_MASS_ENERGY + MOVE + MOMENTUM (parallel kernel)

### Phase 4: Pressure subgraph
- Pipeline: Baroclinic -> Exchange -> Solve -> Convergence (cycle)
- CC_IBM integrated (CC_NO_FLUX, CC_MATCH_VELOCITY_FLUX, CC_COMPUTE_VELOCITY_ERROR)
- Dependency-managed parallel exchange (CODE 5)

### Phase 5: Dump phase + barrier merging
- Fork-join: DumpGlobal || DumpMeshOutputs
- TimestepState: merged join+loop
- BarrierChainTask, dual-output orchestrators, fork optimizations
- DivP1 split (prefork + early/late branches)
- DivP2 preprocessing extraction

## Irreducible Sequential Work

- MESH_EXCHANGE (MPI sync)
- HVAC_CALC (global network solve)
- PRESSURE_ITERATION convergence check (MPI reduction)
- GLOBAL_MATRIX_REASSIGN, EXCHANGE_DIVERGENCE_INFO
- STOP_CHECK, ADJUST_DT (global timestep)
