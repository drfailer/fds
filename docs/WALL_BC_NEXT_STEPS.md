# WALL_BC Parallelization - Next Steps

## Summary

I've completed the foundational work for WALL_BC decomposition:

✅ **Added cross-mesh tracking flags to WALL_TYPE**
- `HAS_INTERPOLATED_BC`: Marks cells requiring OMESH access for interpolated boundaries
- `HAS_BACK_MESH`: Marks cells with back-side mesh coupling for thin walls
- Flags initialized in `FIND_WALL_BACK_INDICES` (init.f90)

## Key Finding

Analysis revealed that WALL_BC's sub-routines (`SURFACE_HEAT_TRANSFER`, `SOLID_HEAT_TRANSFER`, `HEAT_TRANSFER_COEFFICIENT`, etc.) all use `POINT_TO_MESH(NM)` internally, making them NOT thread-safe in their current form.

**Implication**: Creating a parallel kernel requires converting these sub-routines to accept `TYPE(MESH_TYPE)` arguments first — a large refactoring effort (~2000+ lines affected).

## Recommended Path Forward

### Option 1: Incremental Sub-Routine Conversion (High Effort, High Reward)
1. Convert `HEAT_TRANSFER_COEFFICIENT` to accept `TYPE(MESH_TYPE)` (~200 lines)
2. Convert `SURFACE_HEAT_TRANSFER` to accept `TYPE(MESH_TYPE)` (~379 lines)
3. Convert `SOLID_HEAT_TRANSFER` to accept `TYPE(MESH_TYPE)` (~1500 lines)
4. Extract cell-local portions of `CALCULATE_ZZ_F` to wall_kernels.f90
5. Create `WALL_BC_PROCESS_CELLS_KERNEL` using converted routines
6. Integrate with Hedgehog sub-graph

**Estimated effort**: 4-6 hours for conversion + testing
**Estimated speedup**: 2.3× on WALL_BC tasks (reduces 49% sequential bottleneck to ~30%)

### Option 2: Focus on Other Parallelizable Tasks (Lower Effort, Moderate Reward)
The profiling results show other sequential bottlenecks:
- **CorrFinal** (385 ms): MATCH_VELOCITY, VELOCITY_BC
- **PredFinal** (342 ms): MATCH_VELOCITY, VELOCITY_BC
- **CorrRadiation** (356 ms): COMPUTE_RADIATION

These routines may have simpler cross-mesh dependencies that could be easier to parallelize.

### Option 3: Hybrid MPI+Hedgehog (Architectural)
Instead of trying to parallelize WALL_BC within a single process, leverage multiple MPI ranks:
- Each MPI rank handles a subset of meshes
- Each rank runs Hedgehog with `kernelThreads=1` (no intra-rank parallelism for WALL_BC)
- Parallelism comes from multiple ranks running simultaneously

This avoids the thread-safety issue entirely while still improving performance on multi-core systems.

## My Recommendation

Given the complexity of converting WALL_BC's sub-routines and the diminishing returns (WALL_BC is only ~10-13% of total runtime based on profiling), I recommend **Option 3** for near-term gains, with **Option 1** as a longer-term investment.

The flags I've added (`HAS_INTERPOLATED_BC`, `HAS_BACK_MESH`) are still useful — they'll enable Option 1 when/if we pursue it.

## What Would You Like To Do Next?

1. **Proceed with Option 1** - I'll start converting sub-routines to be thread-safe
2. **Explore Option 3** - Set up hybrid MPI+Hedgehog architecture
3. **Focus on other tasks** - Identify and parallelize simpler sequential bottlenecks
4. **Continue current approach** - Create a demonstration three-phase structure (not yet parallelized)

Please let me know how you'd like to proceed.
