# VELOCITY_BC Hedgehog Integration Plan

## Status

✅ **Fortran refactoring complete** (Phases 1-4)
⏭️ **Next**: Hedgehog C++ graph integration

## Completed Work (Phases 1-4)

### Available Subroutines

1. **VELOCITY_BC_PREPROCESSING(M,NM,T,APPLY_TO_ESTIMATED_VARIABLES)**
   - ~77 lines
   - Sequential execution required (OMESH cross-mesh reads)
   - Sets up wall boundary velocities from neighboring meshes
   - Initializes M%DRAG_UVWMAX

2. **VELOCITY_BC_PROCESS_EDGES_KERNEL(M,NM,T,APPLY_TO_ESTIMATED_VARIABLES)**
   - ~762 lines
   - **Parallelizable per-mesh** (thread-safe)
   - Processes all cell edges (EDGE_LOOP)
   - Safe to run in parallel after preprocessing completes
   - Contains ~90% of VELOCITY_BC computation

3. **VELOCITY_BC_KERNEL(M,NM,T,APPLY_TO_ESTIMATED_VARIABLES)**
   - ~53 lines (simplified orchestrator)
   - Calls preprocessing + process_edges
   - Maintains backward compatibility

4. **VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES)**
   - Wrapper for backward compatibility
   - Calls VELOCITY_BC_KERNEL(MESHES(NM),NM,T,APPLY_TO_ESTIMATED_VARIABLES)

## Current Fortran Main Loop Structure

### Predictor Final Section (main.f90:780-787)

```fortran
! MESH_EXCHANGE(3) - velocity/pressure exchange barrier

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL MATCH_VELOCITY(NM)
ENDDO

VELOCITY_BC_LOOP: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (SYNTHETIC_EDDY_METHOD) CALL SYNTHETIC_TURBULENCE(DT,T,NM)
   CALL VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.TRUE.)
ENDDO VELOCITY_BC_LOOP
```

**Sequential time**: ~342 ms (4-mesh test)
**Parallelizable fraction**: ~75-80%

### Corrector Final Section (main.f90:966-975)

```fortran
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL MATCH_VELOCITY(NM)
ENDDO

VELOCITY_BC_LOOP_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)
   CALL UPDATE_GLOBAL_OUTPUTS(T,DT,NM)
ENDDO VELOCITY_BC_LOOP_2

CALL EXCHANGE_GLOBAL_OUTPUTS
```

**Sequential time**: ~385 ms (4-mesh test)
**Parallelizable fraction**: ~70-75%

## Hedgehog Integration Strategy

### Option A: Simple Task-Based Parallelization

**Easiest implementation** - minimal C++ code changes:

1. **Sequential preprocessing phase** (Orchestrator):
   ```cpp
   for (int nm = 0; nm < nmeshes; nm++) {
       match_velocity_(&nm);  // Cross-mesh sync
       velocity_bc_preprocessing_(&meshes[nm], &nm, &t, &apply_to_estimated);
   }
   ```

2. **Parallel kernel phase** (Task):
   ```cpp
   // Dispatched in parallel (one task per mesh)
   if (synthetic_eddy_method) synthetic_turbulence_(&dt, &t, &nm);
   velocity_bc_process_edges_kernel_(&meshes[nm], &nm, &t, &apply_to_estimated);
   ```

3. **Sequential finalization** (Collector):
   ```cpp
   // Optional: UPDATE_GLOBAL_OUTPUTS for corrector
   if (corrector) {
       for (int nm = 0; nm < nmeshes; nm++) {
           update_global_outputs_(&t, &dt, &nm);
       }
       exchange_global_outputs_();
   }
   ```

**Estimated speedup**:
- Predictor: 342ms → ~140ms (2.4× faster)
- Corrector: 385ms → ~170ms (2.3× faster)
- Combined savings: ~420ms on 4-mesh test

### Option B: Dedicated Sub-Graph (Pattern B - like WallBC)

**More complex** - follows established pattern:

Create `Source/hedgehog/graph/velocitybc_subgraph.h`:

```cpp
inline auto buildVelocityBCSubgraph(int nmeshes, size_t kernelThreads,
                                     bool predictor, double t, double dt) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityBC");

    // Orchestrator: MATCH_VELOCITY + VELOCITY_BC_PREPROCESSING
    auto orchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityBCWork>>(
        std::make_shared<VelocityBCOrchestrator>(nmeshes, predictor, t), "VelBCOrch");

    // Kernel Task: VELOCITY_BC_PROCESS_EDGES_KERNEL (+ SYNTHETIC_TURBULENCE if predictor)
    auto kernelTask = std::make_shared<VelocityBCKernelTask>(kernelThreads, predictor, t, dt);

    // Collector: gather results (+ UPDATE_GLOBAL_OUTPUTS if corrector)
    auto collectorSM = std::make_shared<hh::StateManager<1, VelocityBCWork, MeshData>>(
        std::make_shared<VelocityBCCollector>(nmeshes, predictor, t, dt), "VelBCCollector");

    // Wire components
    subgraph->inputs(orchSM);
    subgraph->edges(orchSM, kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}
```

Create supporting files:
- `Source/hedgehog/data/velocitybc_data.h` - VelocityBCWork data structure
- `Source/hedgehog/state/velocitybc_state.h` - Orchestrator + Collector
- `Source/hedgehog/task/velocitybc_kernel_task.h` - Kernel task
- Update `Source/hedgehog/fds_fortran_interface.h` - Add extern "C" declarations

**Benefit**: Named sub-graph appears in dot files for easier profiling/debugging

## Required C++ Fortran Interface Updates

Add to `Source/hedgehog/fds_fortran_interface.h`:

```cpp
extern "C" {
    void velocity_bc_preprocessing_(void* M, int* nm, double* t, bool* apply_to_estimated);
    void velocity_bc_process_edges_kernel_(void* M, int* nm, double* t, bool* apply_to_estimated);
    void match_velocity_(int* nm);
    void synthetic_turbulence_(double* dt, double* t, int* nm);
    void update_global_outputs_(double* t, double* dt, int* nm);
}
```

## Integration Steps

### For Option A (Recommended - quickest win):

1. ✅ Fortran refactoring complete
2. Add VelocityBC task to `Source/hedgehog/task/velocitybc_task.h`
3. Update `fds_graph.h` to replace predictor/corrector VELOCITY_BC loops
4. Update `fds_fortran_interface.h` with new function declarations
5. Compile and test (byte-identical required)
6. Profile 4-mesh test to measure speedup
7. Commit with performance results

**Estimated effort**: 2-3 hours (simple C++ integration)

### For Option B (More elaborate):

1. ✅ Fortran refactoring complete
2. Create `velocitybc_data.h` (work structure)
3. Create `velocitybc_state.h` (orchestrator + collector)
4. Create `velocitybc_kernel_task.h` (kernel task)
5. Create `velocitybc_subgraph.h` (dedicated sub-graph wrapper)
6. Update `fds_graph.h` to use new sub-graph
7. Update `fds_fortran_interface.h`
8. Compile and test (byte-identical required)
9. Profile and commit

**Estimated effort**: 4-6 hours (full Pattern B implementation)

## Expected Performance Impact

**Current bottleneck** (4-mesh, kernelThreads=4):
- Sequential tasks: 1838 ms (39.2%)
- Predictor final: ~342 ms
- Corrector final: ~385 ms
- Total sequential in final sections: ~727 ms

**After integration** (Option A):
- Predictor final: ~140 ms (75% parallelized)
- Corrector final: ~170 ms (75% parallelized)
- Sequential reduction: ~420 ms saved
- New sequential fraction: ~30-32%
- **Target overall speedup**: 1.5-1.7× (up from current 1.3×)

**After full Phase 2** (including PredFinal + CorrFinal + other optimizations):
- Sequential fraction: ~25%
- **Target overall speedup**: 2.0-2.5×

## Testing Protocol

After integration:

1. **Compile**: `cmake --build build_hh --target fds_hh -j$(nproc)`
2. **Run tests**: `cd test_cases && python3 run_tests.py -v`
3. **Verify**: All 5 test cases byte-identical
4. **Profile**: Run 4-mesh test with timing output
5. **Generate dot file**: Check graph structure
6. **Commit**: With performance comparison

## Notes

- MATCH_VELOCITY is inherently sequential (cross-mesh synchronization) - cannot parallelize
- VELOCITY_BC_PREPROCESSING must run sequentially before parallel kernel (OMESH dependencies)
- VELOCITY_BC_PROCESS_EDGES_KERNEL is the main parallelization target (~762 lines, ~90% of work)
- SYNTHETIC_TURBULENCE could be parallelized if converted to thread-safe pattern (future work)
- UPDATE_GLOBAL_OUTPUTS needs analysis for potential parallelization (future work)

## References

- WallBC Pattern B implementation: `Source/hedgehog/graph/wallbc_subgraph.h`
- Existing sub-graphs: `velocitypredictor_subgraph.h`, `velocitycorrector_subgraph.h`
- Main graph: `Source/hedgehog/graph/fds_graph.h`
- Fortran decomposition: `docs/VELOCITY_BC_REFACTORING_PLAN.md`
- Overall progress: `docs/PHASE2_PARALLELIZATION_PROGRESS.md`
