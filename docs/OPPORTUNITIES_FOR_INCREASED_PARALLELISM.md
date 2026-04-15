# Opportunities for Increased Parallelism in the Hedgehog Graph

This document outlines remaining opportunities to increase parallelism in the FDS Hedgehog dataflow graph by optimizing sequential loops within barrier states. The goal is to convert O(N) sequential mesh loops to parallel operations while maintaining numerical equivalence and minimizing overhead.

**Last updated**: 2026-04-15

## Analysis Approach

We examined the predictor and corrector subgraphs for barriers containing sequential loops over meshes. These loops typically:
1. Perform global operations (e.g., `fds_mesh_exchange`, `fds_hvac_calc`) that must remain sequential
2. Execute per-mesh function calls that are independent across meshes and safe for parallelization

The recommended strategy is to:
- Reuse existing kernel task infrastructure where possible
- Integrate parallel work into existing fork structures
- Focus on barriers with significant per-mesh computational work
- Avoid creating excessive new tasks for trivial operations

## Previously Identified — Now Complete

These opportunities from the original analysis have been implemented:

1. ~~WallBCFinalize barriers~~ → ✅ Extracted to WallBCFinalizeKernelTask (parallel per-mesh)
2. ~~REMOVE_PARTICLES + MOVE_PARTICLES barrier~~ → ✅ Merged into ParticleOpsKernelTask
3. ~~CC_IBM "WallBCFin+WallDiv+DivExch" barrier~~ → ✅ WallBCFinalize extracted, particle momentum + divergence in restructured CC_IBM path
4. ~~"Join1+RemoveMove+WallBCOrch" barrier~~ → ✅ Restructured: ParticleOps parallel, Soot||HVAC forked
5. ~~CC_IBM "WallBCFin+MeshExch6a" barrier~~ → ✅ WallBCFinalize extracted

## Remaining Opportunities

### Predictor Subgraph

#### 1. Barrier "DivExch+ZoneOps" (predictor_subgraph.h:109-121)
- **Current:** Per-mesh loop calling `fds_divergence_part_2_preprocessing(md->nm, md->dt)` inside barrier
- **Global operations:** `fds_exchange_divergence_info()`, `fds_global_matrix_reassign(0)`, pressure init/increment
- **Optimization:** Extract per-mesh zone ops to a parallel kernel task, keep global ops sequential
- **Notes:** Non-CC_IBM path only. No cross-mesh writes — thread-safe.
- **Impact:** Parallelizes zone operations for each mesh

#### 2. CC_IBM Barrier "WallDiv+DivExch" (predictor_subgraph.h:159-176)
- **Current:** Two sequential per-mesh loops:
  * Loop 1: `fds_particle_momentum_kernel()` + `fds_divergence_part_1_kernel()` per mesh
  * Loop 2: `fds_divergence_part_2_preprocessing()` per mesh (includes GET_LINKED_VELOCITIES)
- **Global operations:** `fds_exchange_divergence_info()`, `fds_global_matrix_reassign(0)`, pressure init/increment
- **Optimization:**
  - Loop 1: Could use existing ParticleMomentumKernelTask + DivP1 kernel in parallel
  - Loop 2: Contains GET_LINKED_VELOCITIES (cross-mesh writes) — **NOT parallelizable**
- **Impact:** Partial — only Loop 1 is parallelizable

### Corrector Subgraph

#### 3. Barrier "DivExch+ZoneOps" (corrector_subgraph.h:217-229)
- **Current:** Per-mesh loop calling `fds_divergence_part_2_preprocessing(md->nm, md->dt)` inside barrier
- **Global operations:** Same as predictor (#1)
- **Optimization:** Same as predictor (#1) — extract to parallel kernel
- **Notes:** Non-CC_IBM path only. Thread-safe.
- **Impact:** Same as #1

#### 4. CC_IBM Barrier "CorrDivExchange" (corrector_subgraph.h:162-175)
- **Current:** Per-mesh loop calling `fds_divergence_part_2_preprocessing()` (includes GET_LINKED_VELOCITIES)
- **Status:** **NOT parallelizable** — GET_LINKED_VELOCITIES has cross-mesh writes
- **Impact:** None

### CC_IBM Pressure Subgraph — COMPLETE (see PROGRESS_CCIBM_PRESSURE.md)

All 6 items done: GET_LINKED_FV pre-loop init, CC_NO_FLUX (baroclinic + solve),
CC_COMPUTE_VELOCITY_ERROR, FN_OMESH exchange prep, gate removed.

## Implementation Guidelines

### Task Creation Principles
1. **Granularity:** Each kernel task should handle 1-3 related Fortran functions to balance parallelism with overhead
2. **Thread Budget:** Allocate threads from the existing `ThreadBudget` system rather than creating new thread pools
3. **Global Operations:** Keep MPI/global operations sequential within barrier lambdas
4. **Task Integration:**
   - For barriers with 1*N token expectations: Replace lambda loop with direct task launch
   - For barriers with 2*N token expectations: Integrate into existing fork branches
5. **Naming Convention:** Use `[Operation]KernelTask` suffix for consistency

## Priority Recommendations

**High Impact, Low Effort:**
1. Extract `divergence_part_2_preprocessing` per-mesh loop from "DivExch+ZoneOps" barriers (#1, #3)

**Not viable:**
- CC_IBM "CorrDivExchange" GET_LINKED_VELOCITIES loop (#4) — cross-mesh writes
