# Sub-Graph Quick Start Guide

## Overview

This guide provides a streamlined checklist for converting sequential Hedgehog tasks into parallel sub-graphs. Use this as a quick reference when applying the velocity corrector pattern to other FDS modules.

## Prerequisites Checklist

Before starting, verify:

- [ ] Thread-safe kernel exists in `Source/*_kernels.f90`
- [ ] Kernel signature: `SUBROUTINE KERNEL_NAME(M, ...) TYPE(MESH_TYPE), INTENT(INOUT) :: M`
- [ ] Kernel uses NO `POINT_TO_MESH` calls
- [ ] Kernel has NO cross-mesh access (OMESH, MESHES array)
- [ ] Current task calls orchestration routine with `POINT_TO_MESH(NM)`

## Implementation Steps

### 1. Create Work Token (5 minutes)

**File**: `Source/hedgehog/data/<module>_data.h`

```cpp
#ifndef <MODULE>_DATA_H
#define <MODULE>_DATA_H

#include "mesh_data.h"
#include <memory>

struct <Module>Work {
    int nm;
    double t, dt;
    // Add other kernel parameters here
    std::shared_ptr<MeshData> originalMeshData;

    <Module>Work(int nm_, double t_, double dt_, std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};

#endif
```

### 2. Create State Managers (10 minutes)

**File**: `Source/hedgehog/state/<module>_state.h`

```cpp
#ifndef <MODULE>_STATE_H
#define <MODULE>_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/<module>_data.h"

class <Module>Orchestrator : public hh::AbstractState<1, MeshData, <Module>Work> {
public:
    explicit <Module>Orchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, <Module>Work>(), nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing here (if needed)
            for (auto &md : collected_) {
                auto work = std::make_shared<<Module>Work>(md->nm, md->t, md->dt, md);
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

class <Module>Collector : public hh::AbstractState<1, <Module>Work, MeshData> {
public:
    explicit <Module>Collector(int nmeshes)
        : hh::AbstractState<1, <Module>Work, MeshData>(), nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<<Module>Work> work) override {
        results_.push_back(work);
        if (static_cast<int>(results_.size()) == nmeshes_) {
            // Sequential post-processing here (if needed)
            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }
            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<<Module>Work>> results_;
};

#endif
```

### 3. Create Kernel Task (10 minutes)

**File**: `Source/hedgehog/task/<module>_kernel_task.h`

```cpp
#ifndef <MODULE>_KERNEL_TASK_H
#define <MODULE>_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/<module>_data.h"
#include "../fds_fortran_interface.h"

class <Module>KernelTask : public hh::AbstractTask<1, <Module>Work, <Module>Work> {
public:
    explicit <Module>KernelTask(size_t numThreads)
        : hh::AbstractTask<1, <Module>Work, <Module>Work>("<Module>Kernel", numThreads) {}

    void execute(std::shared_ptr<<Module>Work> work) override {
        // Call thread-safe kernel wrappers
        fds_<module>_kernel(work->nm, work->t, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, <Module>Work, <Module>Work>> copy() override {
        return std::make_shared<<Module>KernelTask>(this->numberThreads());
    }
};

#endif
```

### 4. Add Fortran Kernel Wrapper (5 minutes)

**File**: `Source/hedgehog/fds_c_interface.f90`

Add at the bottom of the "Thread-safe kernel wrappers" section:

```fortran
SUBROUTINE C_FDS_<MODULE>_KERNEL(NM, T, DT) BIND(C, NAME="fds_<module>_kernel")
    USE <MODULE>_KERNELS, ONLY: <MODULE>_KERNEL
    USE MESH_VARIABLES, ONLY: MESHES
    INTEGER(C_INT), VALUE :: NM
    REAL(C_DOUBLE), VALUE :: T, DT
    CALL <MODULE>_KERNEL(MESHES(NM), T, DT)
END SUBROUTINE C_FDS_<MODULE>_KERNEL
```

### 5. Add C Interface Declaration (2 minutes)

**File**: `Source/hedgehog/fds_fortran_interface.h`

Add in the "Thread-safe kernel wrappers" section:

```cpp
void fds_<module>_kernel(int nm, double t, double dt);
```

### 6. Integrate into Graph (10 minutes)

**File**: `Source/hedgehog/graph/fds_graph.h`

**Step 6a**: Add includes at top:
```cpp
#include "../data/<module>_data.h"
#include "../task/<module>_kernel_task.h"
#include "../state/<module>_state.h"
```

**Step 6b**: Create components in `buildFDSGraph`:
```cpp
auto <module>OrchSM = std::make_shared<hh::StateManager<1, MeshData, <Module>Work>>(
    std::make_shared<<Module>Orchestrator>(nmeshes), "<Module>Orch");
auto <module>KernelTask = std::make_shared<<Module>KernelTask>(numThreads);
auto <module>CollectorSM = std::make_shared<hh::StateManager<1, <Module>Work, MeshData>>(
    std::make_shared<<Module>Collector>(nmeshes), "<Module>Collector");
```

**Step 6c**: Replace task edges:
```cpp
// OLD:
// graph->edges(prevTask, <module>Task);
// graph->edges(<module>Task, nextTask);

// NEW:
graph->edges(prevTask, <module>OrchSM);
graph->edges(<module>OrchSM, <module>KernelTask);
graph->edges(<module>KernelTask, <module>CollectorSM);
graph->edges(<module>CollectorSM, nextTask);
```

### 7. Build and Test (5 minutes)

```bash
# Build
cd build_hh
cmake --build . --target fds_hh -j$(nproc)

# Test sequential (numThreads=1)
cd ../test_cases/run_1mesh
rm -f *.csv
mpiexec -n 1 ../../build_hh/Source/hedgehog/fds_hh ../dancing_eddies_1mesh_short.fds

# Verify byte-identical
diff dancing_eddies_1mesh_short_devc.csv \
     ../../test_cases/saved_results/hh_1mesh/dancing_eddies_1mesh_short_devc.csv
```

**Expected**: No diff output (byte-identical)

## Common Patterns

### Pattern A: Single Kernel, No Pre/Post Processing

**Example**: CHECK_DIVERGENCE_KERNEL

- Orchestrator: Pure collect → emit
- Task: Single kernel call
- Collector: Pure collect → emit

**Time**: ~30 minutes

### Pattern B: Multiple Kernels

**Example**: VELOCITY_CORRECTOR (velocity kernel + divergence check)

- Task calls multiple kernels in sequence
- All kernels must be thread-safe

**Time**: ~35 minutes

### Pattern C: Kernel with Pre-Processing

**Example**: VELOCITY_FLUX with CC_IBM

- Orchestrator: Collect → call sequential CC_VELOCITY_BC → emit
- Task: Call VELOCITY_FLUX_KERNEL
- Collector: Collect → call sequential CC_VELOCITY_FLUX → emit

**Time**: ~45 minutes

### Pattern D: Kernel with Retry Loop

**Example**: VELOCITY_PREDICTOR with CFL check

- Orchestrator: Collect → emit
- Task: Call kernel, compute DT_NEW(NM)
- Collector: Check global condition
  - If retry: Re-emit work tokens → back to task
  - Else: Emit MeshData → continue

**Time**: ~60 minutes (more complex)

## Troubleshooting

### Build Fails

**Issue**: `undefined reference to fds_<module>_kernel`

**Solution**: Check Fortran wrapper has correct `BIND(C, NAME=...)` and matches C declaration

---

**Issue**: `no matching function for call to make_shared<<Module>Work>`

**Solution**: Check work token constructor signature matches usage

### Test Fails (Non-Byte-Identical)

**Issue**: CSV differs from baseline

**Solution**:
1. Compare kernel vs orchestration Fortran code
2. Verify kernel is truly equivalent (no missing logic)
3. Check PREDICTOR/CORRECTOR flag is set correctly
4. Look for pointer alias issues (US vs U, RHOS vs RHO)

### Graph Deadlocks

**Issue**: Simulation hangs, never completes

**Solution**:
1. Check nmeshes matches actual token count
2. Verify orchestrator/collector emit exactly N tokens
3. Look for missing edges in graph wiring

## Time Estimates

| Module Complexity | Time to Implement | Time to Test |
|-------------------|-------------------|--------------|
| Simple (1 kernel) | 30 minutes        | 5 minutes    |
| Medium (2-3 kernels) | 40 minutes     | 5 minutes    |
| Complex (pre/post) | 60 minutes       | 10 minutes   |
| Advanced (retry) | 90 minutes         | 15 minutes   |

## Priority Order (Recommended)

Based on impact and complexity:

1. ✅ **VELOCITY_CORRECTOR** - DONE (prototype)
2. **DIVERGENCE_PART_2** - Simple, high impact
3. **DIVERGENCE_PART_1** - Simple, high impact
4. **VELOCITY_PREDICTOR** - Complex (retry), highest impact
5. **DENSITY** - Medium, moderate impact
6. **MASS_FINITE_DIFFERENCES** - Simple, moderate impact
7. **VELOCITY_FLUX** - Complex (CC_IBM), moderate impact
8. **Others** - As profiling identifies bottlenecks

## Next Module: Example Walkthrough

To implement DIVERGENCE_PART_2 sub-graph:

1. **Check kernel**: `Source/divg_kernels.f90:1389` - `DIVERGENCE_PART_2_KERNEL(M,DT)`
2. **Create** `divergence_part2_data.h` with `DivergencePart2Work`
3. **Create** `divergence_part2_state.h` with Orchestrator/Collector
4. **Create** `divergence_part2_kernel_task.h`
5. **Add wrapper** in `fds_c_interface.f90`:
   ```fortran
   CALL DIVERGENCE_PART_2_KERNEL(MESHES(NM), DT)
   ```
6. **Declare** in `fds_fortran_interface.h`:
   ```cpp
   void fds_divergence_part_2_kernel(int nm, double dt);
   ```
7. **Wire** in `fds_graph.h` - replace `DivPart2PredTask` and `CorrDivPart2Task`
8. **Build and test** - verify byte-identical

**Time**: ~30 minutes (simple kernel, no pre/post)

## Reference Implementation

See `docs/velocity_corrector_subgraph_implementation.md` for complete working example.

See `docs/hedgehog_subgraph_methodology.md` for detailed methodology and patterns.
