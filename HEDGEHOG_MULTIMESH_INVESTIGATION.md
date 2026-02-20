# Hedgehog Multi-Mesh Crash Investigation

## Date: February 20, 2026

## Problem Summary

The Hedgehog dataflow integration for FDS crashes when running with multiple MPI processes (multi-mesh simulation), but works correctly with a single process.

## Test Results

### ✓ WORKING Cases
1. **Original FDS - 1 mesh**: 9.2s - Completes successfully
2. **Original FDS - 4 meshes** (4 processes): 4.1s - Completes successfully
3. **Hedgehog FDS - 1 mesh**: 9.0s - Completes successfully
4. **Hedgehog FDS - 1 process, 4 meshes**: Completes successfully (single process owns all meshes)

### ✗ FAILING Case
5. **Hedgehog FDS - 4 meshes** (4 processes): Crashes with segfault around time step 10

## Root Cause Analysis

### Initial Bug (FIXED)
**Problem**: Each MPI process was creating a graph for ALL meshes instead of only its local meshes.

**Evidence**:
- In multi-process mode, FDS distributes meshes across processes via `LOWER_MESH_INDEX` and `UPPER_MESH_INDEX`
- The original code called `fds_get_nmeshes()` which returns the TOTAL number of meshes
- Each process was pushing tokens for meshes 1-4, even though each process only owns 1 mesh
- Barriers were configured to wait for `nmeshes=4` tokens per process

**Fix Applied**:
1. Added C interface functions: `fds_get_lower_mesh_index()` and `fds_get_upper_mesh_index()`
2. Calculate `local_nmeshes = upper_mesh_index - lower_mesh_index + 1`
3. Build graph with `local_nmeshes` instead of total `nmeshes`
4. Only push MeshData tokens for meshes owned by current process: `for (nm = lower_mesh_index; nm <= upper_mesh_index; ++nm)`

**Result**: The simulation now starts correctly and all 4 processes recognize their local meshes:
```
Process 0: Local meshes=1 (range: 1-1)
Process 1: Local meshes=1 (range: 2-2)
Process 2: Local meshes=1 (range: 3-3)
Process 3: Local meshes=1 (range: 4-4)
```

### Remaining Bug (ACTIVE)
**Problem**: Simulation crashes with segmentation fault around time step 10.

**Evidence**:
- Simulation runs for ~10 time steps successfully
- Process rank 1 (mesh 2) consistently crashes first
- Crash occurs in FDS code (offset +0x14b220 and +0x6c4a7d in executable)
- Crash signature:
  ```
  Signal: Segmentation fault (11)
  Signal code: Address not mapped (1)
  Failing at address: 0x5bf877dfb080
  ```
- Output files are partially written (only 8 time steps instead of 27 expected)
- Device output file is empty (data not flushed)

**Hypothesis**: MPI synchronization issue

The Hedgehog graph on each process is independent. The barrier states in the graph only synchronize tokens within a single process's graph, NOT across MPI processes.

FDS's global operations (mesh_exchange, pressure_iteration, etc.) use MPI collective operations that require all processes to call them simultaneously. If processes are advancing through the graph at different rates, they can get out of sync, causing:
- MPI deadlocks (one process waiting for message from another)
- Race conditions (accessing shared data structures)
- Memory corruption from out-of-order operations

## Critical Code Locations

### Files Modified
1. `Source/hedgehog/fds_c_interface.f90` - Added mesh index getters
2. `Source/hedgehog/fds_fortran_interface.h` - Added C declarations
3. `Source/hedgehog/main_hh.cpp` - Fixed mesh initialization loop

### Key Graph Components
- `Source/hedgehog/graph/fds_graph.h` - Graph construction
- `Source/hedgehog/state/mesh_barrier_state.h` - Barrier implementations
- `Source/hedgehog/task/predictor_tasks.h` - Predictor task implementations

## Next Steps for Debugging

### Option 1: Add Diagnostic Output
Add print statements in:
- Each task's execute() method to log mesh number and time step
- Barrier states to log when barriers complete
- MPI operations to see synchronization points

### Option 2: Test with 2 Processes
Run with `mpiexec -n 2` to simplify the problem and see if crash pattern changes.

### Option 3: Review MPI Barrier Semantics
Check if fds_mesh_exchange() and other global routines are truly collective operations that require ALL processes to participate.

### Option 4: Add MPI Barriers Between Graph Steps
The Hedgehog graph may need explicit MPI barriers at key synchronization points to ensure all processes stay in lockstep.

### Option 5: Stack Trace Analysis
Use addr2line or gdb to identify exactly which function is crashing:
```bash
addr2line -e build_hh/Source/hedgehog/fds_hh 0x14b220
```

## Architectural Consideration

The fundamental challenge is that Hedgehog provides THREAD-level parallelism (via TBB), while FDS uses PROCESS-level parallelism (via MPI).

Each MPI process runs its own independent Hedgehog graph. The graphs don't communicate with each other - only the FDS Fortran code communicates via MPI. This creates a potential mismatch where:

1. Hedgehog tries to maximize throughput by pipelining tasks
2. FDS expects all processes to execute in strict lockstep for collective operations
3. If one process's graph gets ahead, it may call MPI operations before other processes are ready

**Potential Solutions**:
- Add MPI barriers in critical state managers
- Ensure graph structure prevents any process from getting more than 1 time step ahead
- Review all collective MPI operations and ensure proper synchronization

## Test Infrastructure

Test cases are in `/test_cases/`:
- `dancing_eddies_1mesh_short.fds` - Single mesh (T_END=0.1s)
- `dancing_eddies_4mesh_short.fds` - 4 meshes (T_END=0.1s)
- `run_tests.sh` - Automated test script
- `README.md` - Test documentation

Build commands:
```bash
# Regular build
cd build_hh && cmake -DUSE_OPENMP=OFF .. && cmake --build . --target fds_hh -j$(nproc)

# Debug build
cd build_hh && cmake -DCMAKE_BUILD_TYPE=Debug .. && cmake --build . --target fds_hh -j$(nproc)
```

Test commands:
```bash
# Run all tests
cd test_cases && ./run_tests.sh

# Run 4-mesh hedgehog only
cd test_cases/hh_4mesh && mpiexec -n 4 ../../build_hh/Source/hedgehog/fds_hh dancing_eddies_4mesh_short.fds
```
