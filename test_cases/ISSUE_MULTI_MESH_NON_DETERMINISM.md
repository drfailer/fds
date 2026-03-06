# Issue: Non-Deterministic Results in Multi-Mesh Tests

**Date**: March 6, 2026
**Status**: 🔴 OPEN - Blocking multi-mesh test reliability
**Priority**: HIGH - Affects 2 of 5 test cases
**Related**: Multi-mesh parallel execution within single process (Hedgehog parallelism)

---

## Problem Statement

Multi-mesh tests (dancing_eddies_2mesh, dancing_eddies_4mesh) produce **non-deterministic results** across consecutive runs. The simulations complete successfully but generate different numerical values in output files, causing test comparisons to fail.

**Important Context**: We are NOT using MPI for multi-mesh parallelism. The objective is to process multiple meshes **in parallel within a single process** using Hedgehog's dataflow parallelism framework.

---

## Symptoms

### Failing Tests

1. **dancing_eddies_2mesh** (2 meshes, embedded configuration)
   - Simulation completes successfully in ~18s
   - HRR values differ across runs (up to 3.36e-02 relative error)
   - 104 out of 636 comparisons fail
   - Max difference: 1.62e+01

2. **dancing_eddies_4mesh** (4 meshes, parallel configuration)
   - Simulation completes successfully in ~9s
   - HRR values differ across runs (up to 6.03e+12 relative error for small values)
   - 45 out of 348 comparisons fail
   - Max difference: 4.86e+13

### Example Non-Determinism

```
Row 9, Col 'kW':
  Run 1: -1.4523978E-003
  Run 2: -1.4036232E-003
  Diff: 3.36e-02

Row 5, Col 'kW':
  Run 1: -1.6054519E-016
  Run 2:  9.6829089E-004
  Diff: 6.03e+12 (huge relative error on small values)
```

---

## Working Tests (Deterministic)

- **dancing_eddies_1mesh** (1 mesh): ✅ PASS - Byte-for-byte reproducible
- **multiple_reac_3mesh** (3 meshes): ✅ PASS - Byte-for-byte reproducible

**Key Observation**: The 3-mesh test works deterministically despite being multi-mesh, while 2-mesh and 4-mesh do not. This suggests the issue is test-specific or mesh-configuration-specific, not a fundamental multi-mesh problem.

---

## Current Architecture

### Single-Process Multi-Mesh Execution

The current implementation runs all meshes in a **single MPI process** with Hedgehog managing parallelism:

```cpp
// From main_hh.cpp:
// Execute with -n 1 (single MPI process)
mpiexec -n 1 fds_hh input.fds

// Hedgehog graph processes meshes in parallel:
auto graph = buildFDSGraph(nmeshes, t, dt, tEnd, velCorrKernelThreads);
```

### Parallel Components

Currently, only the **velocity corrector kernel** is parallelized:

```cpp
// From fds_graph.h (line 49):
size_t velCorrKernelThreads = local_nmeshes;  // Parallel velocity kernel
auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(velCorrKernelThreads);
```

**All other tasks are sequential** (1 thread):
- PredStep1Task(1)
- DensityPredTask(1)
- CorrStep1Task(1)
- etc.

### Data Flow

Each mesh has its own `MeshData` token that flows through the graph:

```cpp
// Push one token per mesh
for (int nm = lower_mesh_index; nm <= upper_mesh_index; ++nm) {
    auto md = std::make_shared<MeshData>(nm, t, dt, 0);
    graph->pushData(md);
}
```

---

## Potential Root Causes

### 1. Race Conditions in Fortran Global State

**Hypothesis**: Sequential Hedgehog tasks still allow multiple mesh tokens to be in-flight simultaneously. If Fortran routines access shared global state without proper synchronization, this could cause races.

**Evidence**:
- Most FDS tasks use `POINT_TO_MESH(NM)` which sets global module pointers
- Sequential tasks (1 thread) might not prevent concurrent execution of different meshes
- Fortran global variables (e.g., in OUTPUT_DATA module) may be accessed by multiple mesh computations

**Example from POINT_TO_MESH**:
```fortran
! mesh.f90 - Sets global pointers
SUBROUTINE POINT_TO_MESH(NM)
  M => MESHES(NM)
  RHO => M%RHO
  U => M%U
  ! ... many global pointers set
END SUBROUTINE
```

If Task A for mesh 1 calls `POINT_TO_MESH(1)` while Task B for mesh 2 calls `POINT_TO_MESH(2)`, the global pointers will conflict.

### 2. Mesh Exchange / Barrier Synchronization

**Hypothesis**: Multi-mesh tests have mesh-to-mesh communication (embedded meshes, neighbor data exchange) that may not be properly synchronized in the single-process Hedgehog execution.

**Evidence**:
- dancing_eddies_2mesh uses embedded mesh configuration (fine mesh inside coarse)
- dancing_eddies_4mesh has 4 adjacent meshes that exchange boundary data
- multiple_reac_3mesh works (may have different communication pattern)

**Barrier points** in the graph:
- MESH_EXCHANGE collectors should ensure all meshes reach synchronization point
- But if Fortran mesh exchange logic expects MPI barriers, single-process execution might behave differently

### 3. Floating-Point Operation Ordering

**Hypothesis**: Parallel execution changes the order of floating-point operations, leading to different rounding.

**Evidence**:
- Differences are small but not within machine epsilon
- Some differences are large relative errors on very small values (near-zero HRR)

**Likelihood**: MEDIUM - This could explain small variations but not the large relative errors observed.

### 4. Velocity Corrector Kernel Parallelism

**Hypothesis**: The only parallel task (VelocityCorrectorKernelTask) has thread-safety issues.

**Evidence**:
- This is the only task with > 1 thread
- Uses Hedgehog parallel execution to process multiple meshes concurrently
- Extracted kernel code may still have shared state dependencies

**Code location**: `Source/hedgehog/task/velocity_corrector_kernel_task.h`

### 5. Test-Specific Physics Differences

**Hypothesis**: Dancing eddies tests (2-mesh, 4-mesh) have specific physics that expose concurrency issues, while multiple_reac doesn't.

**Observations**:
- multiple_reac_3mesh: Chemical reactions, very short T_END=0.0005s
- dancing_eddies_2/4mesh: Fluid dynamics, longer T_END=0.1s
- Different physics modules may have different thread-safety properties

---

## Reproduction Steps

### Reproduce Non-Determinism

```bash
cd test_cases

# Run 4-mesh test 3 times
for i in 1 2 3; do
    echo "=== Run $i ==="
    rm -rf run/dancing_eddies_4mesh/
    ./run_tests.py --test dancing_eddies_4mesh
done

# Each run will show different comparison results
```

### Generate New Gold and Compare

```bash
# Generate fresh gold file
./run_tests.py --generate-gold --use-hedgehog-gold --test dancing_eddies_4mesh

# Run test immediately after (should pass if deterministic)
./run_tests.py --test dancing_eddies_4mesh
# Result: FAILS - outputs differ from gold just generated
```

### Direct CSV Comparison

```bash
# Run twice and compare
mkdir -p test_run1 test_run2

cd test_run1
timeout 15 mpiexec -n 1 ../../build_hh/Source/hedgehog/fds_hh \
    ../inputs/dancing_eddies_4mesh_short.fds

cd ../test_run2
timeout 15 mpiexec -n 1 ../../build_hh/Source/hedgehog/fds_hh \
    ../inputs/dancing_eddies_4mesh_short.fds

# Compare
../compare_csv.py \
    test_run1/dancing_eddies_4mesh_short_hrr.csv \
    test_run2/dancing_eddies_4mesh_short_hrr.csv
# Result: Files differ (non-deterministic)
```

---

## Impact

### Affected Components

- ❌ **Multi-mesh testing**: Cannot reliably test 2-mesh or 4-mesh configurations
- ❌ **Regression testing**: Cannot detect actual bugs vs. normal variation
- ❌ **CI/CD integration**: Tests will randomly fail
- ✅ **Single-mesh execution**: Works perfectly (dancing_eddies_1mesh)
- ⚠️ **3-mesh execution**: Works (multiple_reac_3mesh) - unclear why this is different

### Test Results Summary

```
Test                    Status   Time    Deterministic
=====================  ========  ======  =============
dancing_eddies_1mesh   ✅ PASS   9s     YES
dancing_eddies_2mesh   ❌ FAIL   18s    NO (varies run-to-run)
multiple_reac_3mesh    ✅ PASS   5s     YES
dancing_eddies_4mesh   ❌ FAIL   9s     NO (varies run-to-run)
species_props_5mesh    ❌ FAIL   2s     N/A (crashes)
```

**Pass rate**: 2/5 (40%) - Only single-mesh and 3-mesh work

---

## Investigation Plan

### Phase 1: Isolate Root Cause

1. **Test 3-mesh configuration in detail**
   - Why does multiple_reac_3mesh work deterministically?
   - Compare mesh configuration vs. dancing_eddies tests
   - Check if it's physics-dependent (reactions vs. fluid flow)

2. **Disable VelocityCorrectorKernel parallelism**
   ```cpp
   // In fds_graph.h, change:
   size_t velCorrKernelThreads = 1;  // Force sequential instead of local_nmeshes
   ```
   - If this fixes non-determinism → problem is in velocity corrector kernel
   - If still non-deterministic → problem is elsewhere (likely Fortran global state)

3. **Add synchronization logging**
   - Log when each mesh enters/exits critical sections
   - Check if mesh tokens are truly executing sequentially or overlapping
   - Add timestamps to identify concurrent execution

4. **Review Fortran global state access**
   - Audit modules for shared variables: OUTPUT_DATA, GLOBAL_CONSTANTS, etc.
   - Check if Q_DOT, M_DOT accumulation is thread-safe
   - Look for SAVE variables in computational routines

### Phase 2: Fix Strategy Based on Findings

**If cause is Velocity Corrector parallelism:**
- Audit VelocityCorrectorKernelTask for shared state
- Ensure each mesh's data is truly independent
- Add thread-local copies if needed

**If cause is Fortran global state:**
- Option A: Force strict sequential execution (lose parallelism)
- Option B: Add mutex/locks around global state access
- Option C: Complete POINT_TO_MESH removal (long-term solution)

**If cause is mesh exchange synchronization:**
- Review CollectorState implementation
- Ensure barriers truly wait for all meshes
- Check MPI vs. single-process execution differences

### Phase 3: Validation

1. Run each test 10 times, verify byte-for-byte identical results
2. Compare gold files generated on different days
3. Test with different numbers of threads
4. Validate against original FDS (if available)

---

## Workarounds (Temporary)

### For Development

Use only deterministic tests:
```bash
./run_tests.py --test dancing_eddies_1mesh --test multiple_reac_3mesh
```

### For Testing Multi-Mesh

Use looser tolerance (not recommended for CI):
```bash
./run_tests.py --tolerance 1e-6  # Instead of 1e-10
```

### Manual Gold Regeneration

Generate gold immediately before each test run:
```bash
./run_tests.py --generate-gold --use-hedgehog-gold --test dancing_eddies_4mesh
./run_tests.py --test dancing_eddies_4mesh
```

---

## References

### Related Code

- `Source/hedgehog/graph/fds_graph.h:49` - Velocity corrector parallelism configuration
- `Source/hedgehog/task/velocity_corrector_kernel_task.h` - Only parallel task
- `Source/hedgehog/state/collector_state.h` - Mesh synchronization
- `Source/mesh.f90` - POINT_TO_MESH global state
- `test_cases/run_tests.py` - Test runner
- `test_cases/compare_csv.py` - Comparison tool

### Related Documentation

- `test_cases/TERMINATION_FIX.md` - Graph termination issue (resolved)
- `docs/architecture/FDS_ARCHITECTURE.md` - Overall architecture
- `MEMORY.md` - POINT_TO_MESH removal progress

### Test Files

- `test_cases/inputs/dancing_eddies_2mesh.fds` - 2-mesh embedded (non-deterministic)
- `test_cases/inputs/dancing_eddies_4mesh_short.fds` - 4-mesh parallel (non-deterministic)
- `test_cases/inputs/multiple_reac_3mesh.fds` - 3-mesh reactions (deterministic ✓)

---

## Success Criteria

Issue is resolved when:

1. ✅ All 4 multi-mesh tests pass consistently (2-mesh, 3-mesh, 4-mesh, 5-mesh)
2. ✅ Consecutive runs produce byte-for-byte identical output files
3. ✅ Gold files remain valid across multiple test runs
4. ✅ Tests pass with default tolerance (1e-10)
5. ✅ Non-determinism eliminated regardless of thread count

---

## Notes

- The fact that **multiple_reac_3mesh works** suggests the problem is solvable
- Single-mesh tests prove the core Hedgehog integration works correctly
- The issue is specifically about **multi-mesh parallelism within a single process**
- This is NOT about MPI multi-process execution (which we're not using)
- Focus should be on understanding why 3-mesh works when 2-mesh and 4-mesh don't

**Last Updated**: March 6, 2026
**Next Action**: Investigate why multiple_reac_3mesh is deterministic while dancing_eddies_2/4mesh are not
