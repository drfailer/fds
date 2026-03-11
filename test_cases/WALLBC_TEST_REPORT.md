# WallBC Sub-Graph Integration Test Report

**Date**: 2026-03-11
**Branch**: hedgehog-integration
**Component**: WallBC Pattern B Sub-Graph

## Summary

✅ **All tests passed** - WallBC sub-graph successfully integrated and verified

The WallBC routine has been successfully parallelized using Hedgehog's Pattern B architecture (sequential preprocessing + parallel kernel + sequential finalization). All test cases produce byte-identical results compared to baseline, confirming correctness of the three-phase decomposition.

## Test Results

### Automated Test Suite (run_tests.py)

All 5 test cases passed with byte-identical output:

| Test Case | Meshes | Runtime | Status |
|-----------|--------|---------|--------|
| dancing_eddies_1mesh | 1 | 5.94s | ✅ PASS |
| dancing_eddies_2mesh | 2 | 13.09s | ✅ PASS |
| dancing_eddies_4mesh | 4 | 5.58s | ✅ PASS |
| multiple_reac_3mesh | 3 | 6.20s | ✅ PASS |
| species_props_5mesh | 5 | 1.05s | ✅ PASS |

**Total**: 5 tests, 5 passed, 0 failed

### Manual Verification Tests

Additional manual tests performed during development:

1. **1-mesh test (dancing_eddies_1mesh_short)**:
   - Runtime: ~0.1s simulation time (27 time steps)
   - Comparison: `_devc.csv` and `_hrr.csv` byte-identical to orig_1mesh baseline
   - Status: ✅ PASS

2. **4-mesh test (dancing_eddies_4mesh_short)**:
   - Runtime: ~0.1s simulation time (27 time steps)
   - Comparison: `_devc.csv` and `_hrr.csv` byte-identical to orig_4mesh baseline
   - Status: ✅ PASS

## WallBC Three-Phase Architecture

### Phase 1: Sequential Preprocessing (Orchestrator)

**Implementation**: `WallBCOrchestrator::execute()` in `wallbc_state.h`

**Operations**:
- Collects all N MeshData tokens
- Computes global parameters: DT_BC (from BC_CLOCK) and CALL_HT_1D (from WALL_COUNTER)
- Updates BC_CLOCK if calling 1-D heat transfer
- Runs sequential preprocessing for each mesh:
  - `ASSIGN_GHOST_VALUE`: Handles INTERPOLATED_BOUNDARY cells (OMESH reads)
  - `NEAR_SURFACE_GAS_VARIABLES_KERNEL`: Setup for all wall cells
  - `HEAT_TRANS_COEF`: Computation for thermally-thick surfaces

**Why sequential**: Cross-mesh dependencies via OMESH reads

### Phase 2: Parallel Kernel Execution

**Implementation**: `WallBCKernelTask::execute()` calling `WALL_BC_PROCESS_CELLS_KERNEL`

**Operations** (155 lines, wall.f90:1370-1524):
- Processes ~90% of wall cells (skips HAS_INTERPOLATED_BC and HAS_BACK_MESH)
- Calls thread-safe routines:
  - `SURFACE_HEAT_TRANSFER`
  - `CALCULATE_ZZ_F`
  - `CALC_HVAC_BC`
  - `HEAT_TRANSFER_COEFFICIENT`
- Handles CFACE cells and particles
- Independent per-mesh execution

**Why parallelizable**: No cross-mesh dependencies, cell-local processing

### Phase 3: Sequential Finalization (Collector)

**Implementation**: `WallBCCollector::execute()` in `wallbc_state.h`

**Operations**:
- Gathers all N WallBCWork tokens
- Sorts results by mesh index for deterministic ordering
- Runs sequential finalization for each mesh:
  - `HAS_BACK_MESH` cells (thin walls spanning meshes)
  - Thin wall lateral heat transfer
  - `DEPOSIT_PARTICLE_MASS` (particle off-gassing, updates neighboring mesh via OMESH)

**Why sequential**: Cross-mesh writes and BACK_MESH coupling

## Implementation Details

### Files Modified

1. **Source/wall.f90**:
   - Added `WALL_BC_PREPROCESSING` (49 lines, 1597-1645)
   - Added `WALL_BC_PROCESS_CELLS_KERNEL` (155 lines, 1370-1524)
   - Added `WALL_BC_FINALIZE` (68 lines, 1527-1594)
   - Restructured main `WALL_BC` to three-phase architecture (lines 141-163)

2. **Source/hedgehog/fds_c_interface.f90**:
   - Added C wrapper: `C_FDS_WALL_BC_PREPROCESSING` (BIND(C))
   - Added C wrapper: `C_FDS_WALL_BC_PROCESS_CELLS_KERNEL` (BIND(C), RECURSIVE)
   - Added C wrapper: `C_FDS_WALL_BC_FINALIZE` (BIND(C))
   - Added helper: `C_FDS_COMPUTE_WALL_BC_DT_BC` (computes T - BC_CLOCK)
   - Added helper: `C_FDS_CHECK_CALL_HT_1D` (checks WALL_COUNTER == WALL_INCREMENT)
   - Added helper: `C_FDS_UPDATE_BC_CLOCK` (updates BC_CLOCK = T)

3. **Source/hedgehog/fds_fortran_interface.h**:
   - Added 6 C declarations for Fortran functions

4. **Source/hedgehog/data/wallbc_data.h**:
   - Created `WallBCWork` struct with nm, t, dt, dt_bc, call_ht_1d, originalMeshData

5. **Source/hedgehog/state/wallbc_state.h**:
   - Created `WallBCOrchestrator` (Pattern B orchestrator)
   - Created `WallBCCollector` (Pattern B collector)

6. **Source/hedgehog/task/wallbc_kernel_task.h**:
   - Created `WallBCKernelTask` (parallel kernel execution)

7. **Source/hedgehog/graph/fds_graph.h**:
   - Added WallBC sub-graph component declarations (lines 194-202)
   - Wired sub-graph into corrector flow (lines 386-390)
   - Replaced sequential `corrWallBC` task

### Thread Safety

All callee routines verified thread-safe:
- ✅ `CALC_HVAC_BC` (52 lines, thread-safe conversion completed)
- ✅ `HEAT_TRANSFER_COEFFICIENT` (~175 lines, thread-safe conversion completed)
- ✅ `SURFACE_HEAT_TRANSFER` (379 lines, thread-safe conversion completed)
- ✅ `CALCULATE_ZZ_F` (413 lines, thread-safe conversion completed)

WALL_BC_PROCESS_CELLS_KERNEL marked RECURSIVE for thread safety.

## Performance Characteristics

**Estimated parallelization benefit**: ~90% of WALL_BC processing can run concurrently across meshes

**Key parameters computed once per time step**:
- `DT_BC` = T - BC_CLOCK (boundary condition time step)
- `CALL_HT_1D` = (WALL_COUNTER == WALL_INCREMENT) (1-D heat transfer flag)

**Cross-mesh coordination**:
- Phase 1: OMESH reads (INTERPOLATED_BOUNDARY cells) - ~5-10% of cells
- Phase 3: OMESH writes (particle off-gassing, BACK_MESH coupling) - ~5% of cells

## Integration Status

**Sub-graph count**: 12 total sub-graphs in hedgehog-integration branch

Previous sub-graphs (all verified):
1. Velocity Corrector
2. Velocity Predictor
3. DivPart2 (predictor + corrector)
4. CorrStep1 (viscosity + mass FD + density)
5. DensityPred
6. CorrDivPart1
7. DivSetup (predictor + corrector)
8. PredStep1
9. CorrCondens
10. PredWallDiv (WALL_BC sequential + momentum+div1 parallel)
11. CorrParticle

**New**: 12. WallBC (Pattern B: preprocessing + parallel kernel + finalization)

## Conclusions

1. **Correctness**: All test cases produce byte-identical results, confirming the three-phase decomposition correctly handles all cross-mesh dependencies.

2. **Thread Safety**: All prerequisite thread-safe conversions completed and verified. WALL_BC_PROCESS_CELLS_KERNEL safely executes in parallel.

3. **Architecture**: Pattern B (sequential pre/post-processing with parallel kernel) successfully applied to WALL_BC. The three-phase architecture cleanly separates:
   - Cross-mesh reads (preprocessing)
   - Independent cell processing (parallel kernel)
   - Cross-mesh writes (finalization)

4. **Maintainability**: Clear separation of concerns with well-defined interfaces. Each phase has a specific responsibility and can be tested independently.

5. **Scalability**: The ~90% parallelizable fraction enables significant speedup potential for multi-mesh simulations.

## Next Steps

Recommended future parallelization targets:
- RADIATION (corrRadiation) - large sequential task in corrector
- Additional predictor/corrector bottlenecks identified via profiling

## References

- Implementation plan: `docs/WALL_BC_PARALLELIZATION_PLAN.md`
- Thread-safe conversions: `docs/WALL_BC_CONVERSIONS_SUMMARY.md`
- Three-phase decomposition: `docs/WALL_BC_DECOMPOSITION.md`
- Sub-graph methodology: `docs/METHOD_SUBGRAPH.md`
