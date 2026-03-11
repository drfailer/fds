# WALL_BC Decomposition - Implementation Notes

## Current Status (2026-03-10)

### Completed
✅ Task 1: Added HAS_INTERPOLATED_BC and HAS_BACK_MESH flags to WALL_TYPE (type.f90:455-456)
✅ Task 1: Initialized flags in FIND_WALL_BACK_INDICES (init.f90:3778-3793)

### Challenge Discovered

Many routines called by WALL_BC use `POINT_TO_MESH(NM)` internally and access module-level variables:
- `SURFACE_HEAT_TRANSFER(NM, ...)` - uses POINT_TO_MESH
- `SOLID_HEAT_TRANSFER(NM, ...)` - uses POINT_TO_MESH
- `CALCULATE_ZZ_F(T, DT, ...)` - uses module-level pointers (RHOP, ZZP, etc.)
- `HEAT_TRANSFER_COEFFICIENT(NM, ...)` - uses POINT_TO_MESH

These routines are NOT thread-safe in their current form.

### Revised Approach

Instead of creating a fully thread-safe kernel immediately, we'll take a phased approach:

**Phase A (Current): Structure Setup**
1. Add flags to identify cross-mesh cells ✅
2. Create basic three-phase structure in WALL_BC
   - Phase 1: ASSIGN_GHOST_VALUE (sequential, OMESH reads)
   - Phase 2: Bulk processing (will call existing routines, skip flagged cells)
   - Phase 3: Cross-mesh finalization (INTERPOLATED_BC, BACK_MESH, CONSUME_MASS)

**Phase B (Next): Incremental Thread-Safety**
1. Convert SURFACE_HEAT_TRANSFER to take TYPE(MESH_TYPE)
2. Convert SOLID_HEAT_TRANSFER to take TYPE(MESH_TYPE)
3. Convert HEAT_TRANSFER_COEFFICIENT to take TYPE(MESH_TYPE)
4. Extract cell-local portions of CALCULATE_ZZ_F

**Phase C (Future): Full Parallelization**
1. Create wall_bc_kernels.f90 with WALL_BC_PROCESS_CELLS_KERNEL
2. Integrate with Hedgehog sub-graph
3. Benchmark and verify byte-identical results

## Why Not Split Routines Now?

- SURFACE_HEAT_TRANSFER: 379 lines, complex INTERPOLATED_BC case (158 lines)
- SOLID_HEAT_TRANSFER: ~1500 lines, deep BACK_MESH integration
- CALCULATE_ZZ_F: ~420 lines, CONSUME_MASS section interleaved with other logic

Splitting these requires careful analysis to avoid breaking existing behavior. Better to:
1. Get the structure right first (flags + three-phase orchestration)
2. Test that the reorganization is byte-identical
3. Then incrementally convert sub-routines to be thread-safe

## Current Implementation Strategy

For now, I'll reorganize WALL_BC into three phases that call the existing (non-thread-safe) routines. This establishes the structure and lets us verify correctness before tackling thread safety.

Once the structure is proven, we can tackle Phase B by converting routines one at a time to accept TYPE(MESH_TYPE) arguments.
