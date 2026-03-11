# Parallel Execution Findings and Thread-Safety Analysis

**Date**: March 6, 2026
**Status**: ⚠️ **THREAD-SAFETY ISSUES IDENTIFIED**

## Executive Summary

Attempted to enable parallel multi-mesh execution (numThreads=N) in the velocity corrector sub-graph. **Segmentation faults occurred** due to Fortran thread-safety issues when calling kernels from multiple C++ threads simultaneously.

**Current Status**:
- ✅ **Sequential execution** (numThreads=1): **WORKS PERFECTLY**
- ❌ **Parallel execution** (numThreads>1): **SEGFAULTS**

## Test Results

### ✅ Test 1: Single Mesh, Parallel Enabled (numThreads=1)

Since single mesh means only 1 token, numThreads=1 effectively runs sequentially even with parallel setting.

**Result**: ✅ **BYTE-IDENTICAL** with baseline

### ❌ Test 2: Four Meshes, Parallel Enabled (numThreads=4)

**Command**:
```bash
# main_hh.cpp: size_t numThreads = local_nmeshes; (=4)
mpiexec -n 1 fds_hh dancing_eddies_4mesh_short.fds
```

**Result**: ❌ **SEGMENTATION FAULT**

**Stack Trace**:
```
[FDS-HH] Initialization complete.
[FDS-HH] Total meshes=4 Local meshes=4 (range: 1-4) t=0 dt=0.00399162 tEnd=0.1
[FDS-HH] Waiting for graph termination...
[pn115038:1060117] *** Process received signal ***
[pn115038:1060117] Signal: Segmentation fault (11)
[pn115038:1060117] Signal code: Address not mapped (1)
```

Segfault occurs immediately after graph execution starts (during first time step).

### ❌ Test 3: Four Meshes with OMP_NUM_THREADS=1

**Hypothesis**: Maybe nested parallelism (Hedgehog + OpenMP) causes the crash.

**Test**:
```bash
OMP_NUM_THREADS=1 mpiexec -n 1 fds_hh dancing_eddies_4mesh_short.fds
```

**Result**: ❌ **STILL SEGFAULTS**

Disabling OpenMP threads did not resolve the issue. The problem is not nested parallelism.

### ❌ Test 4: Four Meshes with numThreads=2

**Hypothesis**: Maybe 4 threads is too many.

**Result**: ❌ **STILL SEGFAULTS**

Issue occurs with any numThreads > 1.

## Root Cause Analysis

### Issue 1: OpenMP Directives in Kernels

**Location**: `Source/velo_kernels.f90:189-221`

```fortran
SUBROUTINE VELOCITY_CORRECTOR_KERNEL(M,DT)
    ...
    !$OMP PARALLEL PRIVATE(I,J,K)
    !$OMP DO SCHEDULE(STATIC)
    DO K=1,M%KBAR
        ...
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
END SUBROUTINE
```

**Problem**: When Hedgehog spawns N threads and each calls this kernel, we get nested parallelism:
- Hedgehog thread pool: N threads
- Each thread tries to spawn OpenMP threads
- Total: N × OpenMP_threads

**Impact**: This contributes to the crash but is not the sole cause (OMP_NUM_THREADS=1 still crashes).

### Issue 2: Global Module Variable Access (PRIMARY CAUSE)

**Location**: `Source/divg_kernels.f90:1669-1679`

```fortran
SUBROUTINE CHECK_DIVERGENCE_KERNEL(M)
    ...
    IF (PREDICTOR) THEN          ! <-- GLOBAL MODULE VARIABLE
        UU=>M%US
        ...
    ELSEIF (CORRECTOR) THEN      ! <-- GLOBAL MODULE VARIABLE
        UU=>M%U
        ...
    ENDIF

    SELECT CASE(CYLINDRICAL)     ! <-- GLOBAL MODULE VARIABLE
        CASE(.FALSE.)
            DIV = ...
        CASE(.TRUE.)
            DIV = ...
    END SELECT
    ...
```

**Globals accessed**:
- `PREDICTOR` - module variable from GLOBAL_CONSTANTS
- `CORRECTOR` - module variable from GLOBAL_CONSTANTS
- `CYLINDRICAL` - module variable from GLOBAL_CONSTANTS
- `STORE_CARTESIAN_DIVERGENCE` - module variable
- `FREEZE_VELOCITY` - module variable (in VELOCITY_CORRECTOR_KERNEL)

**Problem**:
1. Fortran module variables are **global state**
2. When multiple C++ threads call these kernels simultaneously, they all read from the same global variables
3. Fortran's module system and runtime may not be thread-safe when accessed from multiple foreign (C++) threads
4. Even though these are read-only accesses, Fortran's memory model may not guarantee thread-safety

### Issue 3: Fortran Runtime Thread-Safety

**Problem**: FDS was not compiled with thread-safe Fortran runtime flags.

**Typical thread-safe compilation** would require:
```bash
# GNU Fortran
gfortran -frecursive -fno-automatic ...

# Intel Fortran
ifort -recursive -reentrancy threaded ...
```

**Current FDS build**: Does not include explicit thread-safety flags for multi-threaded C++ calls.

### Issue 4: MESHES Array Access

**Location**: `Source/hedgehog/fds_c_interface.f90`

```fortran
SUBROUTINE C_FDS_VELOCITY_CORRECTOR_KERNEL(NM, T, DT) BIND(C, ...)
    USE MESH_VARIABLES, ONLY: MESHES    ! <-- MODULE VARIABLE (ARRAY)
    ...
    CALL VELOCITY_CORRECTOR_KERNEL(MESHES(NM), DT)  ! <-- ARRAY ACCESS FROM MULTIPLE THREADS
END SUBROUTINE
```

**Problem**:
- `MESHES` is a module-level allocatable array
- Multiple threads access `MESHES(NM)` with different NM values
- While accessing different array elements should be safe, Fortran's array descriptor handling may not be thread-safe

## Why Sequential Works

**Sequential execution** (numThreads=1):
- Only 1 Hedgehog thread exists
- No concurrent access to global module variables
- No concurrent access to MESHES array
- OpenMP threads within kernel are fine (standard OpenMP usage)

**Result**: Perfect byte-identical results.

## Attempted Solutions (All Failed)

1. **Disable OpenMP**: `OMP_NUM_THREADS=1` → Still segfaults
2. **Reduce threads**: `numThreads=2` → Still segfaults
3. **Single mesh**: Works, but numThreads=1 anyway (only 1 token)

## Architectural Implications

### Current Architecture Limitations

The velocity corrector sub-graph was designed for **data parallelism** (process N meshes in parallel). However:

**Fortran Constraints**:
- Kernels were extracted to remove `POINT_TO_MESH`
- But kernels still access global module variables
- Fortran module system is not designed for multi-threaded C++ calls
- FDS codebase assumes single-threaded execution or OpenMP-only parallelism

**Hedgehog-Fortran Impedance Mismatch**:
- Hedgehog: Modern C++ thread pool, expects thread-safe tasks
- FDS Fortran: Traditional HPC code, assumes MPI or OpenMP parallelism only
- Mixing paradigms causes undefined behavior

### What Works

**Current working configuration**:
```cpp
size_t numThreads = 1;  // Sequential per-mesh processing
```

**Parallelism achieved**:
- ✅ **MPI parallelism**: Multiple processes, each with numThreads=1
- ✅ **OpenMP parallelism**: Within each kernel, OpenMP threads parallelize loops
- ❌ **Hedgehog data parallelism**: Cannot process multiple meshes concurrently

**Performance**:
- Single mesh: No benefit from Hedgehog parallelism (only 1 token anyway)
- Multi-mesh, single MPI process: Sequential mesh processing (slow)
- Multi-mesh, multiple MPI processes: MPI parallelism + OpenMP within kernels (current FDS model)

## Recommended Solutions

### Option 1: Accept Sequential Hedgehog (SHORT-TERM)

**Status**: ✅ **CURRENT IMPLEMENTATION**

**Approach**:
- Keep `numThreads = 1` in all tasks
- Use Hedgehog for **pipeline parallelism** only (not data parallelism)
- Rely on MPI + OpenMP for actual parallelism

**Pros**:
- Works correctly now
- Byte-identical results
- Clean architecture (orchestrator/kernel/collector pattern)

**Cons**:
- No performance benefit from Hedgehog on single MPI process
- Multi-mesh cases still sequential per mesh

**Use Case**:
- Validates architecture
- Establishes pattern for future modules
- Works well with MPI (each process has its own graph)

### Option 2: Thread-Safe Fortran Compilation (MEDIUM-TERM)

**Approach**:
1. Recompile FDS with thread-safe Fortran flags:
   ```cmake
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecursive -fno-automatic")
   ```
2. Ensure Fortran runtime is thread-safe
3. Test with `numThreads > 1`

**Pros**:
- May enable true parallel execution
- Minimal code changes

**Cons**:
- Not guaranteed to work (Fortran standard doesn't require thread-safety)
- May impact performance (recursive allocation overhead)
- Requires full rebuild and testing

**Risk**: MODERATE (may not solve all issues)

### Option 3: Remove Global Variable Access (LONG-TERM)

**Approach**:
1. Pass PREDICTOR, CORRECTOR, CYLINDRICAL as kernel arguments
2. Store these in MESH_TYPE structure
3. Eliminate all global module variable reads in kernels

**Example**:
```fortran
SUBROUTINE CHECK_DIVERGENCE_KERNEL(M, IS_PREDICTOR, IS_CYLINDRICAL)
    TYPE(MESH_TYPE), INTENT(INOUT) :: M
    LOGICAL, INTENT(IN) :: IS_PREDICTOR, IS_CYLINDRICAL
    ...
    IF (IS_PREDICTOR) THEN
        UU => M%US
    ELSE
        UU => M%U
    ENDIF
    ...
END SUBROUTINE
```

**Pros**:
- True thread-safety
- Clean functional programming style
- Enables Hedgehog data parallelism

**Cons**:
- Significant refactoring required
- Changes all kernel signatures
- Must update all call sites

**Effort**: HIGH (weeks of work)

### Option 4: Hybrid Approach (ALTERNATIVE)

**Approach**:
- Use Hedgehog for **pipeline parallelism** (different stages concurrently)
- Use **OpenMP** for **data parallelism** (within each kernel)
- Keep `numThreads = 1` for mesh-processing tasks

**Architecture**:
```
Pipeline Stage 1 (Thread 1): Mesh 1 → VELOCITY_CORRECTOR_KERNEL (spawns OpenMP threads)
    ↓
Pipeline Stage 2 (Thread 2): Mesh 1 → DIVERGENCE_KERNEL (spawns OpenMP threads)

Meanwhile:
Pipeline Stage 1 (Thread 1): Mesh 2 → VELOCITY_CORRECTOR_KERNEL (spawns OpenMP threads)
```

**Pros**:
- Works with current code
- Achieves parallelism through pipelining
- No thread-safety issues

**Cons**:
- Different parallelism model than originally intended
- Limited by pipeline depth

**Use Case**: Small number of large meshes

## Performance Analysis

### Current Performance (Sequential)

**Single Mesh**:
- Hedgehog overhead: Minimal
- Performance: Identical to original FDS (OpenMP within kernels)

**Four Meshes, Single MPI Process**:
- Mesh processing: Sequential (Mesh1 → Mesh2 → Mesh3 → Mesh4)
- Within each mesh: OpenMP parallelism
- Total: Same as original FDS sequential mesh processing

**Four Meshes, Four MPI Processes** (current FDS model):
- Each process: 1 mesh, OpenMP parallelism
- MPI exchanges: Between processes
- Total: Full parallelism via MPI

### Potential Performance (If Parallel Worked)

**Four Meshes, Single MPI Process, numThreads=4**:
- Mesh processing: Parallel (all 4 meshes simultaneously)
- Within each mesh: OpenMP parallelism
- Total: 4× speedup on mesh-level operations

**Theoretical Speedup**:
- If velocity operations are 30% of time step: ~1.4× overall
- If velocity operations are 50% of time step: ~2.0× overall

**Reality**: ❌ Crashes, so 0× speedup

## Conclusions

### Key Findings

1. ✅ **Sub-graph architecture is sound** - Sequential execution proves the pattern works
2. ❌ **Fortran thread-safety blocks parallel execution** - Global module variables prevent multi-threading
3. ✅ **Correctness validated** - Byte-identical results with baseline
4. ⚠️ **Performance limited** - Sequential mesh processing (same as original)

### Current Value Proposition

**What the velocity corrector sub-graph provides**:
- ✅ Clean architectural pattern (orchestrator/kernel/collector)
- ✅ Separation of concerns (data flow vs computation)
- ✅ Foundation for future parallelization
- ✅ Works correctly with MPI parallelism
- ❌ **No performance benefit** on single MPI process with multiple meshes

### Future Work

**Short-term** (Recommended):
1. Document thread-safety limitations
2. Complete other sub-graph conversions (DIVERGENCE, etc.)
3. Use Hedgehog for pipeline parallelism only
4. Rely on MPI + OpenMP for data parallelism

**Medium-term** (If desired):
1. Attempt thread-safe Fortran compilation
2. Test with numThreads > 1
3. Profile to see if worth the effort

**Long-term** (Major effort):
1. Refactor kernels to remove global variable access
2. Pass all configuration as arguments
3. Enable true thread-safe data parallelism
4. Achieve within-node multi-mesh parallelism

## Recommendation

**For Production**: **Keep numThreads=1** (sequential)

**Rationale**:
- Works correctly
- No risk of crashes
- Maintains MPI + OpenMP parallelism model
- Provides architectural benefits without performance regression
- Allows systematic conversion of other modules
- Future-proofed for when thread-safety can be achieved

**For Research**: **Investigate Option 2** (thread-safe compilation)
- Low effort, moderate chance of success
- Could unlock parallelism without code changes
- Worth trying if performance is critical

**For Long-term**: **Plan Option 3** (remove globals)
- High effort, guaranteed success
- Aligns with modern parallel programming practices
- Enables full Hedgehog capabilities
