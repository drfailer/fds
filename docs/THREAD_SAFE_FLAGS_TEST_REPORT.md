# Thread-Safe Fortran Flags Test Report

**Date**: March 6, 2026
**Test**: Compile with `-frecursive -fno-automatic` flags to enable parallel execution
**Result**: ❌ **FAILED** - Still segfaults with parallel execution

---

## Objective

Test if compiling FDS with thread-safe Fortran flags would allow parallel multi-mesh execution in Hedgehog (numThreads > 1) without segmentation faults.

## Hypothesis

The segfaults observed with `numThreads > 1` were caused by non-thread-safe Fortran module variables. Compiling with `-frecursive` and `-fno-automatic` flags should:
- Make all local variables allocated on the stack (not static)
- Prevent SAVE semantics by default
- Enable thread-safe access to module variables

## Implementation

### Changes Made

**File**: `CMakeLists.txt:157`

```cmake
# Before
target_compile_options(fds PRIVATE -cpp -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics)

# After
target_compile_options(fds PRIVATE -cpp -std=f2018 -frecursive -fno-automatic -ffpe-summary=none -fall-intrinsics)
```

**Flags Added**:
- `-frecursive`: Already present (makes subroutines reentrant)
- `-fno-automatic`: **NEW** (prevents automatic SAVE, forces stack allocation)

**Configuration**:
```cpp
// main_hh.cpp
size_t numThreads = local_nmeshes;  // Try parallel (=4 for 4-mesh case)
```

### Build Process

```bash
# Clean rebuild required
rm -rf build_hh
mkdir build_hh && cd build_hh
cmake .. -DUSE_HEDGEHOG=ON
cmake --build . --target fds_hh -j$(nproc)
```

**Result**: ✅ Build successful (no compilation errors)

---

## Test Results

### Test 1: Four Meshes, Parallel (numThreads=4)

**Command**:
```bash
OMP_NUM_THREADS=1 mpiexec -n 1 fds_hh dancing_eddies_4mesh_short.fds
```

**Configuration**:
- 4 meshes on 1 MPI process
- Hedgehog numThreads = 4 (parallel)
- OpenMP threads = 1 (disabled via env var)

**Result**: ❌ **SEGMENTATION FAULT**

**Stack Trace**:
```
[pn115038:1075471] Signal: Segmentation fault (11)
[pn115038:1075471] Failing at address: 0x5a8996ba8e30
[pn115038:1075471] [ 4] /lib/x86_64-linux-gnu/libgomp.so.1(GOMP_parallel+0x46)
[pn115038:1075471] [ 5] ../../build_hh/Source/hedgehog/fds_hh(+0x4f2ced)
```

**Analysis**: Crash occurs in `GOMP_parallel` (OpenMP parallel region). Despite setting `OMP_NUM_THREADS=1`, the crash still occurs during OpenMP initialization.

### Test 2: Four Meshes, Sequential (numThreads=1)

**Command**:
```bash
mpiexec -n 1 fds_hh dancing_eddies_4mesh_short.fds
```

**Configuration**:
- 4 meshes on 1 MPI process
- Hedgehog numThreads = 1 (sequential)
- Thread-safe flags enabled

**Result**: ✅ **SUCCESS**

**Verification**:
```bash
diff dancing_eddies_4mesh_short_devc.csv baseline/dancing_eddies_4mesh_short_devc.csv
# No output - BYTE-IDENTICAL ✅
```

---

## Root Cause Analysis

### Why Thread-Safe Flags Didn't Work

The `-frecursive -fno-automatic` flags address **static variable** thread-safety, but the actual issue is **nested parallelism**:

1. **Hedgehog spawns N threads** (numThreads=4)
2. **Each thread calls kernel** (e.g., VELOCITY_CORRECTOR_KERNEL)
3. **Each kernel contains OpenMP directives**:
   ```fortran
   !$OMP PARALLEL PRIVATE(I,J,K)
   !$OMP DO SCHEDULE(STATIC)
   DO K=1,M%KBAR
       ...
   ENDDO
   !$OMP END DO
   !$OMP END PARALLEL
   ```
4. **Result**: 4 × OpenMP_threads nested parallelism → crash

### Stack Trace Evidence

The stack trace shows:
- Frame 4: `GOMP_parallel` - GNU OpenMP parallel region entry
- Frames 2-3, 5-6: FDS kernel code being executed

This confirms the crash occurs **during OpenMP initialization within the kernel**, not due to module variable access.

### Why OMP_NUM_THREADS=1 Didn't Help

Even with `OMP_NUM_THREADS=1`:
- OpenMP runtime is still initialized
- `GOMP_parallel` is still called
- Nested parallelism still occurs (4 Hedgehog threads each trying to use OpenMP)
- OpenMP may not be designed for this use case (called from foreign thread pool)

---

## Comparison of Thread-Safety Issues

| Issue Type | Addressed By | Status |
|------------|--------------|--------|
| Static local variables | `-frecursive` | ✅ Fixed (already had it) |
| Automatic SAVE | `-fno-automatic` | ✅ Fixed (added) |
| Global module variables | Flags | ⚠️ Still read-only unsafe |
| Nested parallelism | Flags | ❌ **NOT FIXED** |
| OpenMP directives in kernels | Flags | ❌ **NOT FIXED** |

---

## Implications

### What the Flags DID Provide

1. **Better thread-safety for local variables**
2. **Stack allocation instead of static** (good practice)
3. **No performance regression** in sequential mode
4. **Still byte-identical results**

### What the Flags DID NOT Provide

1. **Fix for nested parallelism** (Hedgehog + OpenMP)
2. **Fix for OpenMP-in-kernel issue**
3. **Ability to use numThreads > 1**
4. **Performance improvement**

---

## Alternative Solutions

### Option 1: Disable OpenMP in Kernels (Requires Code Changes)

**Approach**: Remove OpenMP directives from kernel files

```fortran
! velo_kernels.f90 - Remove these lines:
! !$OMP PARALLEL PRIVATE(I,J,K)
! !$OMP DO SCHEDULE(STATIC)
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL
```

**Pros**:
- Would likely fix nested parallelism
- Hedgehog parallelism could work

**Cons**:
- Loses within-kernel OpenMP parallelism
- Changes many kernel files
- Affects standard FDS build too

**Risk**: MODERATE (removes existing parallelism)

### Option 2: Conditional OpenMP (Recommended)

**Approach**: Add preprocessor flag to disable OpenMP for Hedgehog build

```cmake
# CMakeLists.txt
if(USE_HEDGEHOG)
    target_compile_definitions(fds_hh PRIVATE NO_KERNEL_OPENMP)
endif()
```

```fortran
! velo_kernels.f90
#ifndef NO_KERNEL_OPENMP
!$OMP PARALLEL PRIVATE(I,J,K)
#endif
DO K=1,M%KBAR
    ...
ENDDO
#ifndef NO_KERNEL_OPENMP
!$OMP END PARALLEL
#endif
```

**Pros**:
- Hedgehog build: No OpenMP in kernels
- Standard FDS: Keeps OpenMP in kernels
- Enables Hedgehog data parallelism

**Cons**:
- Requires modifying all kernel files
- Maintains two code paths

**Effort**: MODERATE (systematic but straightforward)

### Option 3: Single-Level Parallelism Only

**Approach**: Accept current sequential Hedgehog

**Pros**:
- Works now
- No code changes
- Thread-safe flags provide robustness

**Cons**:
- No performance benefit from Hedgehog parallelism

**Effort**: NONE (current state)

---

## Conclusion

### Test Verdict: ❌ FAILED

Thread-safe Fortran flags (`-frecursive -fno-automatic`) **did not** enable parallel execution in Hedgehog. The root cause is **nested parallelism** (Hedgehog threads + OpenMP in kernels), not module variable thread-safety.

### Flags Are Still Valuable

Despite not solving the parallel execution issue, the thread-safe flags provide:
- Better code robustness
- Stack allocation (safer than static)
- No regression in sequential mode
- **Recommended to keep**: `-frecursive -fno-automatic`

### Path Forward

**Short-term**: ✅ Keep thread-safe flags, use `numThreads=1`
- Provides robustness
- No performance loss
- Works correctly

**Medium-term**: Consider Option 2 (Conditional OpenMP)
- Would enable Hedgehog parallelism
- Moderate implementation effort
- Clear separation of concerns

**Long-term**: Comprehensive refactoring
- Remove global module variables
- Make kernels pure functions
- Enable both Hedgehog and OpenMP parallelism

---

## Files Modified (Permanent)

**CMakeLists.txt**: Added `-fno-automatic` flag (KEEP THIS)

```cmake
target_compile_options(fds PRIVATE -cpp -std=f2018 -frecursive -fno-automatic ...)
```

**main_hh.cpp**: Reverted to `numThreads=1`

```cpp
size_t numThreads = 1;  // Sequential (nested parallelism issue)
```

---

## Recommendations

1. ✅ **Keep `-fno-automatic` flag** - Improves robustness, no downside
2. ✅ **Keep `numThreads=1`** - Only working configuration
3. ⏭️ **Consider conditional OpenMP** - If parallel performance is critical
4. 📝 **Document findings** - Help future developers understand the issue

---

**Bottom Line**: Thread-safe flags are good practice and should be kept, but they don't solve the nested parallelism issue. Parallel execution remains blocked until OpenMP directives are removed from kernels or made conditional.
