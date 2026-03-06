# ChangeTimeStep Sub-Graph Refactoring

**Date**: March 6, 2026
**Purpose**: Split monolithic ChangeTimeStepTask into focused tasks with dedicated state management

---

## Summary

The `ChangeTimeStepTask` previously contained the entire CFL retry loop logic (~70 lines) in a single execute() method. This has been refactored into a dedicated sub-graph with 11 small, focused tasks and a state-managed retry loop.

**Result**: ✅ **SUCCESS** - Test passes (27 timesteps completed)

---

## Architecture

### Before: Monolithic Task

```cpp
class ChangeTimeStepTask {
    void execute(BarrierData) {
        // Check if retry needed
        while (needRetry) {
            // 60+ lines of retry logic inline
            // - Restore UVW
            // - Density calculation
            // - Velocity flux
            // - HVAC, divergence, pressure, etc.
            // - Check again
        }
    }
};
```

**Problems:**
- Too many responsibilities in one task
- Retry loop managed imperatively (while loop)
- Hard to understand data flow
- No visibility into retry iterations

### After: Dedicated Sub-Graph

```
BarrierData
    ↓
CheckRetryTask (determines if retry needed)
    ├──done=true──→ RetryExitTask → MeshData (bypass)
    └──done=false─→ RetryDensityTask (start retry sequence)
                        ↓
                    RetryCCDensityTask
                        ↓
                    RetryVelocityFluxTask
                        ↓
                    RetryHvacTask
                        ↓
                    RetryInitDivTask
                        ↓
                    RetryDivergencePart1Task
                        ↓
                    RetryDivExchangeTask
                        ↓
                    RetryDivergencePart2Task
                        ↓
                    RetryPressureTask
                        ↓
                    RetryVelocityPredictorTask
                        ↓
                    RetryLoopState (check if another retry needed)
                    ├──done=false─→ RetryDensityTask (cycle)
                    └──done=true──→ RetryExitTask → MeshData (exit)
```

**Benefits:**
- Each task has a single, clear responsibility
- Retry loop managed declaratively (dataflow)
- Visible retry iterations (each cycle through the graph)
- Easy to extend or modify individual steps

---

## Components

### Data Structures

**`RetrySequenceData`** (data/change_timestep_data.h)
- Carries mesh data through retry sequence
- Fields: meshes, t, dt, iteration, done
- `done` flag controls bypass (true) vs. retry (false)

### Tasks

**`CheckRetryTask`** (task/change_timestep_tasks.h)
- Entry point: receives BarrierData
- Calls `fds_check_change_time_step()`
- Emits RetrySequenceData with done=true (no retry) or done=false (retry)

**Retry Sequence Tasks** (all sequential, numThreads=1)
1. **RetryDensityTask**: Restore UVW + density calculation
2. **RetryCCDensityTask**: CC_DENSITY + MESH_EXCHANGE(1)
3. **RetryVelocityFluxTask**: Velocity flux computation
4. **RetryHvacTask**: HVAC calculation
5. **RetryInitDivTask**: Initialize divergence integrals
6. **RetryDivergencePart1Task**: Wall BC + particle + divergence part 1
7. **RetryDivExchangeTask**: Exchange divergence info
8. **RetryDivergencePart2Task**: Divergence part 2
9. **RetryPressureTask**: Pressure iteration
10. **RetryVelocityPredictorTask**: Velocity predictor

**`RetryExitTask`** (task/change_timestep_tasks.h)
- Exit point: converts RetrySequenceData → MeshData
- Only processes when done=true (filters out done=false)

### State

**`RetryLoopState`** (state/change_timestep_state.h)
- Manages retry loop decision
- After retry sequence completes:
  - Checks stop status
  - Calls `fds_check_change_time_step()` again
  - Emits done=false (another retry) or done=true (exit)

### Sub-Graph Builder

**`buildChangeTimeStepSubgraph()`** (graph/change_timestep_subgraph.h)
- Constructs: `Graph<1, BarrierData, MeshData>`
- Wires all tasks and states
- Returns shared_ptr to complete sub-graph

---

## Key Design Patterns

### Bypass Path

When no retry is needed, the token bypasses the retry sequence entirely:

```cpp
// CheckRetryTask emits to BOTH paths
subgraph->edges(checkRetry, retryDensity);  // Main path
subgraph->edges(checkRetry, retryExit);     // Bypass path

// RetryDensityTask drops done=true tokens
if (data->done) { return; }  // Let retryExit handle it

// RetryExitTask processes done=true tokens
if (data->done) { emit MeshData; }
```

### Loop Management

The retry loop cycles through the graph declaratively:

```cpp
// RetryLoopState decides: cycle or exit
if (needRetry) {
    data->done = false;  // Cycle back to retryDensity
} else {
    data->done = true;   // Forward to retryExit
}

// Both edges exist from retryLoopSM
subgraph->edges(retryLoopSM, retryDensity);  // Cycle
subgraph->edges(retryLoopSM, retryExit);     // Exit
```

The `done` flag acts as a routing signal:
- `done=false`: retryDensity processes, retryExit drops
- `done=true`: retryDensity drops, retryExit processes

---

## Integration with Main Graph

The sub-graph is integrated into the main FDS graph as a single node:

```cpp
// In buildFDSGraph() (graph/fds_graph.h)

// Create sub-graph
auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph();

// Wire into main graph (same as before)
graph->edges(velPredictor, changeTimeStepCollectorSM);
graph->edges(changeTimeStepCollectorSM, changeTimeStepSubgraph);
graph->edges(changeTimeStepSubgraph, collector3SM);
```

Externally, it behaves exactly like the old `ChangeTimeStepTask`:
- Input: BarrierData (all meshes collected)
- Output: MeshData (individual mesh tokens)
- Semantics: Check CFL, retry if needed, emit when compliant

---

## Files Created

### New Files
- `Source/hedgehog/data/change_timestep_data.h` - RetrySequenceData definition
- `Source/hedgehog/task/change_timestep_tasks.h` - All retry tasks
- `Source/hedgehog/state/change_timestep_state.h` - RetryLoopState
- `Source/hedgehog/graph/change_timestep_subgraph.h` - Sub-graph builder

### Modified Files
- `Source/hedgehog/graph/fds_graph.h` - Use sub-graph instead of task
- `Source/hedgehog/task/barrier_tasks.h` - (old ChangeTimeStepTask can be removed)

---

## Testing

**Test case**: `dancing_eddies_1mesh_short.fds`
**Command**: `mpiexec -n 1 fds_hh dancing_eddies_1mesh_short.fds`
**Result**: ✅ Completes 27 timesteps successfully

**Output**:
```
Time Step:       1, Simulation Time: 0.0039916 s
Time Step:       2, Simulation Time: 0.0079832 s
...
Time Step:      27, Simulation Time: 0.1000000 s
[FDS-HH] Graph terminated.
STOP: FDS completed successfully
```

---

## Next Steps

### Completed
- ✅ Extract retry logic into sub-graph
- ✅ Create focused tasks for each operation
- ✅ Implement state-managed loop
- ✅ Test with 1-mesh case

### Future Work
1. **Remove old ChangeTimeStepTask** from barrier_tasks.h (now obsolete)
2. **Test with 4-mesh case** (parallel velocity kernel)
3. **Add retry counter** to track number of CFL retries per timestep
4. **Consider parallelizing per-mesh operations** in retry sequence (if beneficial)
5. **Apply same pattern to other complex tasks** (if any)

---

## Lessons Learned

1. **Dual-edge routing requires filtering** - When a state emits to multiple edges, downstream tasks must filter based on flags (done field).

2. **Bypass paths need direct connections** - To skip a sequence, emit to both the sequence entry (which drops) and the exit (which processes).

3. **Sub-graphs compose naturally** - A sub-graph can replace a single task in the main graph without changing the external interface.

4. **Declarative loops are clearer** - Managing loops via dataflow (cycle edges) is more transparent than imperative while loops.

5. **Small tasks are easier to debug** - When something goes wrong, it's clear which step failed.

---

## Bottom Line

The ChangeTimeStep refactoring successfully demonstrates how to extract complex logic from a monolithic task into a clean, composable sub-graph with state-managed loops. The new architecture is more maintainable, easier to understand, and provides better visibility into execution flow.
