# Hedgehog Sub-Graph Methodology

## Purpose

This document describes the systematic methodology for converting sequential, single-mesh Hedgehog tasks into parallel, multi-mesh sub-graphs. This pattern enables multi-mesh parallel processing within a single node by replacing orchestration routines that use `POINT_TO_MESH` with sub-graphs that call thread-safe computation kernels directly.

## Core Concept

**Before (Sequential Processing)**:
```
Task A (mesh 1) → Task A (mesh 2) → ... → Task A (mesh N) → Barrier
```
- Each mesh token processed sequentially through the task
- Task calls orchestration routine with `POINT_TO_MESH(NM)`
- Total time: N × T_task

**After (Parallel Processing via Sub-Graph)**:
```
[Orchestrator State] collects all N tokens
    ↓
[Parallel Kernel Task] processes all N meshes concurrently (numThreads=N)
    ↓
[Collector State] gathers results, emits N tokens
```
- All mesh tokens processed in parallel
- Kernel tasks call thread-safe kernels: `KERNEL(MESHES(NM), ...)`
- Total time: T_task + overhead

## Prerequisites

Before creating a sub-graph, verify:

1. **Thread-safe kernel exists**: Check `*_kernels.f90` for routine taking `TYPE(MESH_TYPE), INTENT(INOUT) :: M`
2. **No cross-mesh access in kernel**: Kernel operates only on `M%` arrays, no `OMESH` or `MESHES` array access
3. **Orchestration identified**: Separate orchestration logic (POINT_TO_MESH, CC_IBM, special cases) from kernel call
4. **C bindings available**: Check `fds_fortran_interface.h` and `fds_c_interface.f90`

## Step-by-Step Methodology

### Step 1: Identify the Target Task

**Example**: CorrVelocityTask (corrector_tasks.h:146-160)

```cpp
class CorrVelocityTask : public hh::AbstractTask<1, MeshData, MeshData> {
    void execute(std::shared_ptr<MeshData> data) override {
        fds_velocity_corrector(data->t, data->dt, data->nm);  // Orchestration
        fds_check_divergence(data->nm);                        // Orchestration
        this->addResult(data);
    }
};
```

**Analysis**:
- Inputs: MeshData (single token)
- Outputs: MeshData (single token)
- Calls: Two Fortran orchestration routines
- Current behavior: Sequential (one mesh at a time)

### Step 2: Identify the Kernels

**Location**: Source/velo_kernels.f90, Source/divg_kernels.f90

```fortran
! velo_kernels.f90:176
SUBROUTINE VELOCITY_CORRECTOR_KERNEL(M,DT)
TYPE(MESH_TYPE), INTENT(INOUT) :: M
REAL(EB), INTENT(IN) :: DT
! ... operates only on M% arrays ...
END SUBROUTINE

! divg_kernels.f90:1660
SUBROUTINE CHECK_DIVERGENCE_KERNEL(M)
TYPE(MESH_TYPE), INTENT(INOUT) :: M
! ... operates only on M% arrays ...
END SUBROUTINE
```

**Verification**:
- ✓ Takes `TYPE(MESH_TYPE)` as argument
- ✓ No `POINT_TO_MESH` calls
- ✓ No cross-mesh access

### Step 3: Create Work Token Data Structure

**File**: `Source/hedgehog/data/velocity_corrector_data.h`

```cpp
#ifndef VELOCITY_CORRECTOR_DATA_H
#define VELOCITY_CORRECTOR_DATA_H

#include "mesh_data.h"
#include <memory>

/// Work token for parallel velocity corrector kernel execution
struct VelocityCorrectorWork {
    int nm;          ///< Mesh index
    double t;        ///< Simulation time
    double dt;       ///< Time step

    // Preserve original MeshData for downstream routing
    std::shared_ptr<MeshData> originalMeshData;

    VelocityCorrectorWork(int nm_, double t_, double dt_,
                          std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif // VELOCITY_CORRECTOR_DATA_H
```

**Why separate work token?**
- Decouples kernel execution from graph routing
- Allows kernel-specific parameters
- Preserves original MeshData for downstream passage

### Step 4: Create Orchestrator State

**File**: `Source/hedgehog/state/velocity_corrector_state.h`

**Purpose**: Collect N MeshData tokens, perform sequential pre-processing, emit N work tokens

```cpp
#ifndef VELOCITY_CORRECTOR_STATE_H
#define VELOCITY_CORRECTOR_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/velocity_corrector_data.h"

/// Orchestrator state for velocity corrector sub-graph.
/// Collects all N mesh tokens and dispatches parallel work.
class VelocityCorrectorOrchestrator
    : public hh::AbstractState<1, MeshData, VelocityCorrectorWork> {
public:
    explicit VelocityCorrectorOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, VelocityCorrectorWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // All mesh tokens collected - dispatch parallel work

            // Sequential pre-processing would go here
            // (e.g., CC_IBM setup if needed)

            // Emit work tokens for parallel kernel execution
            for (auto &md : collected_) {
                auto work = std::make_shared<VelocityCorrectorWork>(
                    md->nm, md->t, md->dt, md);
                this->addResult(work);
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // VELOCITY_CORRECTOR_STATE_H
```

**Key points**:
- Collects exactly N tokens before dispatching
- Sequential operations happen here (not in parallel task)
- Emits N work tokens (one per mesh)

### Step 5: Create Parallel Kernel Task

**File**: `Source/hedgehog/task/velocity_corrector_kernel_task.h`

**Purpose**: Execute thread-safe kernels in parallel for multiple meshes

```cpp
#ifndef VELOCITY_CORRECTOR_KERNEL_TASK_H
#define VELOCITY_CORRECTOR_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/velocity_corrector_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls thread-safe velocity corrector kernels.
/// Each thread processes one mesh independently.
class VelocityCorrectorKernelTask
    : public hh::AbstractTask<1, VelocityCorrectorWork, VelocityCorrectorWork> {
public:
    explicit VelocityCorrectorKernelTask(size_t numThreads)
        : hh::AbstractTask<1, VelocityCorrectorWork, VelocityCorrectorWork>(
              "VelocityCorrectorKernel", numThreads) {}

    void execute(std::shared_ptr<VelocityCorrectorWork> work) override {
        // Call thread-safe kernels directly (no POINT_TO_MESH)
        fds_velocity_corrector_kernel(work->nm, work->t, work->dt);
        fds_check_divergence_kernel(work->nm);

        // Pass work token downstream (contains original MeshData)
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, VelocityCorrectorWork, VelocityCorrectorWork>>
    copy() override {
        return std::make_shared<VelocityCorrectorKernelTask>(this->numberThreads());
    }
};

#endif // VELOCITY_CORRECTOR_KERNEL_TASK_H
```

**Key points**:
- `numThreads` parameter controls parallelism (use nmeshes for full parallelism)
- Calls kernel wrappers that bypass orchestration
- No POINT_TO_MESH, no global state access
- Each thread operates on different mesh (indexed by `work->nm`)

### Step 6: Create Collector State

**File**: `Source/hedgehog/state/velocity_corrector_state.h` (add to same file)

**Purpose**: Collect N result tokens, perform sequential post-processing, emit N MeshData tokens

```cpp
/// Collector state for velocity corrector sub-graph.
/// Gathers all N kernel results and emits MeshData tokens downstream.
class VelocityCorrectorCollector
    : public hh::AbstractState<1, VelocityCorrectorWork, MeshData> {
public:
    explicit VelocityCorrectorCollector(int nmeshes)
        : hh::AbstractState<1, VelocityCorrectorWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<VelocityCorrectorWork> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            // All kernel results collected

            // Sequential post-processing would go here
            // (e.g., global diagnostics, reductions)

            // Emit original MeshData tokens to continue graph flow
            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<VelocityCorrectorWork>> results_;
};
```

**Key points**:
- Waits for all N kernels to complete
- Sequential post-processing (reductions, diagnostics)
- Emits original MeshData tokens (preserves graph flow)

### Step 7: Add Kernel Wrapper Functions (Fortran)

**File**: `Source/hedgehog/fds_c_interface.f90`

Add wrappers that call kernels directly without `POINT_TO_MESH`:

```fortran
SUBROUTINE C_FDS_VELOCITY_CORRECTOR_KERNEL(NM, T, DT) BIND(C, NAME="fds_velocity_corrector_kernel")
    USE VELO_KERNELS, ONLY: VELOCITY_CORRECTOR_KERNEL
    USE MESH_VARIABLES, ONLY: MESHES
    REAL(C_DOUBLE), VALUE :: T, DT
    INTEGER(C_INT), VALUE :: NM
    CALL VELOCITY_CORRECTOR_KERNEL(MESHES(NM), DT)
END SUBROUTINE C_FDS_VELOCITY_CORRECTOR_KERNEL

SUBROUTINE C_FDS_CHECK_DIVERGENCE_KERNEL(NM) BIND(C, NAME="fds_check_divergence_kernel")
    USE DIVG_KERNELS, ONLY: CHECK_DIVERGENCE_KERNEL
    USE MESH_VARIABLES, ONLY: MESHES
    INTEGER(C_INT), VALUE :: NM
    CALL CHECK_DIVERGENCE_KERNEL(MESHES(NM))
END SUBROUTINE C_FDS_CHECK_DIVERGENCE_KERNEL
```

**Key differences from orchestration wrappers**:
- No `POINT_TO_MESH` call
- Directly passes `MESHES(NM)` to kernel
- No conditional logic (CC_IBM, cylindrical, etc.)

### Step 8: Declare C Interface (C++)

**File**: `Source/hedgehog/fds_fortran_interface.h`

```cpp
// Thread-safe kernel wrappers (bypass orchestration)
extern "C" {
    void fds_velocity_corrector_kernel(int nm, double t, double dt);
    void fds_check_divergence_kernel(int nm);
}
```

### Step 9: Integrate Sub-Graph into Main Graph

**File**: `Source/hedgehog/graph/fds_graph.h`

Replace the sequential task with the sub-graph pattern:

```cpp
// OLD (lines 218-220):
// graph->edges(corrPressureTask, corrVelocity);
// graph->edges(corrVelocity, collector6bSM);

// NEW:
// Create sub-graph components
auto velCorrOrchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityCorrectorWork>>(
    std::make_shared<VelocityCorrectorOrchestrator>(nmeshes), "VelCorrOrch");
auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(numThreads);
auto velCorrCollectorSM = std::make_shared<hh::StateManager<1, VelocityCorrectorWork, MeshData>>(
    std::make_shared<VelocityCorrectorCollector>(nmeshes), "VelCorrCollector");

// Wire sub-graph
graph->edges(corrPressureTask, velCorrOrchSM);          // Pressure → Orchestrator
graph->edges(velCorrOrchSM, velCorrKernelTask);         // Orchestrator → Kernel
graph->edges(velCorrKernelTask, velCorrCollectorSM);    // Kernel → Collector
graph->edges(velCorrCollectorSM, collector6bSM);        // Collector → MeshExchange(6)
```

**Graph structure**:
```
corrPressureTask (N tokens)
    ↓
[VelCorrOrchestrator] collects N, emits N work tokens
    ↓
[VelocityCorrectorKernelTask] parallel execution (numThreads)
    ↓
[VelCorrCollector] collects N results, emits N MeshData
    ↓
collector6bSM (for MESH_EXCHANGE(6))
```

### Step 10: Include New Headers

**File**: `Source/hedgehog/graph/fds_graph.h`

Add at top:
```cpp
#include "../data/velocity_corrector_data.h"
#include "../state/velocity_corrector_state.h"
#include "../task/velocity_corrector_kernel_task.h"
```

### Step 11: Update CMake (if needed)

If you created new header files, ensure they're accessible. Hedgehog headers are typically header-only, so no CMakeLists changes needed for headers.

## Testing Methodology

### Phase 1: Sequential Verification (numThreads=1)

**Purpose**: Verify byte-identical results with sequential execution

```cpp
size_t numThreads = 1;  // Sequential for correctness check
auto graph = buildFDSGraph(local_nmeshes, t, dt, tEnd, numThreads);
```

**Steps**:
1. Build with sub-graph enabled
2. Run test case (e.g., `dancing_eddies_1mesh_short.fds`)
3. Compare CSV output with baseline:
   ```bash
   diff test_cases/run_1mesh/dancing_eddies_1mesh_short_devc.csv \
        test_cases/saved_results/hh_1mesh/dancing_eddies_1mesh_short_devc.csv
   ```
4. Verify byte-identical (no diff output)

**If results differ**: Kernel is not equivalent to orchestration routine - investigate:
- Are CC_IBM paths handled correctly?
- Are PREDICTOR/CORRECTOR flags set correctly?
- Are pointer aliases resolved correctly?

### Phase 2: Parallel Validation (numThreads=N)

**Purpose**: Verify thread-safety and parallel correctness

```cpp
size_t numThreads = local_nmeshes;  // Parallel execution
auto graph = buildFDSGraph(local_nmeshes, t, dt, tEnd, numThreads);
```

**Steps**:
1. Build with parallel execution
2. Run same test case
3. Compare CSV output with sequential baseline
4. Verify byte-identical

**If results differ**: Thread-safety violation - investigate:
- Are kernels truly independent (no shared state)?
- Are indexed global writes correct (e.g., `DT_NEW(NM)`)?
- Are Fortran I/O statements thread-safe?

### Phase 3: Multi-Mesh Testing

**Purpose**: Verify sub-graph works with multiple meshes

```bash
mpiexec --oversubscribe -n 4 build_hh/Source/hedgehog/fds_hh dancing_eddies_4mesh_short.fds
```

**Verify**:
- All meshes process correctly
- No deadlocks (graph completes)
- Results match baseline

### Phase 4: Performance Measurement

**Purpose**: Quantify parallel speedup

```bash
# Sequential
time mpiexec -n 1 fds_hh test.fds  # numThreads=1

# Parallel
time mpiexec -n 1 fds_hh test.fds  # numThreads=N
```

**Expected**: ~1.5-3× speedup for compute-intensive kernels (diminishing returns beyond core count)

## Common Patterns and Variations

### Pattern 1: Simple Kernel (No Pre/Post Processing)

**Example**: CHECK_DIVERGENCE
- Orchestrator: Collect N, emit N (no logic)
- Kernel Task: Call kernel
- Collector: Collect N, emit N (no logic)

### Pattern 2: Kernel with Pre-Processing

**Example**: VELOCITY_FLUX with CC_IBM
- Orchestrator: Collect N, call `CC_VELOCITY_BC` (sequential), emit N
- Kernel Task: Call `VELOCITY_FLUX_KERNEL`
- Collector: Collect N, call `CC_VELOCITY_FLUX` (sequential), emit N

### Pattern 3: Kernel with Global Reduction

**Example**: CFL time step check
- Orchestrator: Collect N, emit N
- Kernel Task: Call kernel, write to `DT_NEW(NM)` (indexed)
- Collector: Collect N, compute `MINVAL(DT_NEW)` (sequential), emit N

### Pattern 4: Kernel with Retry Loop

**Example**: VELOCITY_PREDICTOR with CFL retry
- Orchestrator: Collect N, emit N work tokens
- Kernel Task: Call kernel, compute DT_NEW(NM)
- Collector: Check global CFL condition
  - If retry: Emit new work tokens → loop back to kernel task
  - Else: Emit N MeshData → continue graph

**Implementation**:
```cpp
class CFLRetryCollector : public hh::AbstractState<...> {
    void execute(...) override {
        results_.push_back(work);
        if (results_.size() == nmeshes_) {
            int needRetry;
            double newDt;
            fds_check_change_time_step(&needRetry, &newDt);

            if (needRetry) {
                // Re-emit work tokens with smaller DT
                for (auto &w : results_) {
                    w->dt = newDt;
                    this->addResult(w);  // Back to kernel task
                }
            } else {
                // Success - emit MeshData to continue
                for (auto &w : results_) {
                    this->addResult(w->originalMeshData);
                }
            }
            results_.clear();
        }
    }
};
```

**Graph wiring**: Collector → Kernel Task (for retry loop)

## Troubleshooting Guide

### Issue: Non-Byte-Identical Results (Sequential)

**Cause**: Kernel not equivalent to orchestration routine

**Debug**:
1. Compare kernel vs orchestration Fortran code line-by-line
2. Check for pointer alias differences (US vs U, RHOS vs RHO)
3. Verify PREDICTOR/CORRECTOR flag is set correctly
4. Check for missing CC_IBM conditional logic

### Issue: Non-Deterministic Results (Parallel)

**Cause**: Thread-safety violation (race condition)

**Debug**:
1. Check for unindexed global writes
2. Look for `POINT_TO_MESH` remnants
3. Verify Fortran I/O is not in kernel
4. Check for module-level SAVE variables

### Issue: Deadlock (Graph Doesn't Terminate)

**Cause**: Collector not emitting tokens

**Debug**:
1. Add debug prints to collector state
2. Verify `nmeshes_` matches actual token count
3. Check if some meshes are filtered (SOLID_PHASE_ONLY)
4. Verify all paths emit tokens (retry loop must emit)

### Issue: Poor Speedup (Parallel < 1.5×)

**Cause**: Kernel is too lightweight or has sequential bottleneck

**Debug**:
1. Profile kernel execution time
2. Check if orchestration pre/post-processing dominates
3. Verify numThreads matches actual parallelism
4. Check for false sharing (adjacent mesh data in cache lines)

## Summary Checklist

For each sub-graph implementation:

- [ ] Verify thread-safe kernel exists
- [ ] Create work token data structure
- [ ] Implement orchestrator state (collect N → emit N work)
- [ ] Implement parallel kernel task (calls kernel, no POINT_TO_MESH)
- [ ] Implement collector state (collect N results → emit N MeshData)
- [ ] Add Fortran kernel wrappers (bypass orchestration)
- [ ] Declare C interface
- [ ] Wire sub-graph into main graph
- [ ] Test sequential (numThreads=1) - byte-identical
- [ ] Test parallel (numThreads=N) - byte-identical
- [ ] Test multi-mesh case
- [ ] Measure performance

## Next Steps

After velocity corrector prototype:
1. Apply pattern to VELOCITY_PREDICTOR (with CFL retry)
2. Apply to DIVERGENCE_PART_1, DIVERGENCE_PART_2
3. Apply to DENSITY, MASS, other kernel-based modules
4. Profile to identify remaining sequential bottlenecks
5. Consider finer-grained sub-graphs (e.g., FVX/FVY/FVZ loops in velocity kernel)
