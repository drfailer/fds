# Sub-Graph Extraction Methodology

## Goal

Convert sequential, single-mesh Hedgehog tasks into parallel, multi-mesh sub-graphs. Each sub-graph replaces one sequential task with an orchestrator → parallel kernel → collector pattern, enabling concurrent processing of multiple meshes within a single process.

## Architecture

```
[Orchestrator State]  collects N MeshData tokens, runs sequential pre-processing
        ↓
[Parallel Kernel Task]  processes N meshes concurrently (kernelThreads)
        ↓
[Collector State]  gathers N results, sorts by mesh index, emits N MeshData tokens
```

- **Orchestrator**: always single-threaded (states are not parallelized in Hedgehog)
- **Kernel Task**: `numThreads = kernelThreads` (controlled by `buildFDSGraph` parameter)
- **Collector**: always single-threaded, sorts results by `nm` for deterministic ordering

## Graph Naming Convention

**IMPORTANT**: When creating sub-graphs, use a named graph that matches the Fortran routine being parallelized. This allows tracing graph sections back to original code.

```cpp
// Create named graph for WALL_BC parallelization
auto wallBCGraph = std::make_shared<hh::Graph<...>>("WallBC");

// Add nodes to the named graph
wallBCGraph->input<MeshData>(orchestrator);
wallBCGraph->addEdge(orchestrator, kernelTask);
wallBCGraph->addEdge(kernelTask, collector);
wallBCGraph->output<MeshData>(collector);

// Add the sub-graph to main graph
mainGraph->addEdge(prevNode, wallBCGraph);
mainGraph->addEdge(wallBCGraph, nextNode);
```

**Naming guidelines**:
- Use PascalCase matching Fortran routine name: `"WallBC"`, `"VelocityPredictor"`, `"Radiation"`
- For complex routines with multiple phases, use descriptive names: `"WallBCPhase2"`, `"RadiationSolve"`
- Maintain consistency with existing graphs (see completed sub-graphs in `MEMORY.md`)

**Benefits**:
- Easier debugging: graph name appears in Hedgehog logs
- Code traceability: graph name → Fortran routine
- Documentation: graph structure mirrors Fortran call hierarchy

## Prerequisites

Before creating a sub-graph:

1. **Thread-safe kernel exists** in `*_kernels.f90` with signature `(M, ...)` where `M` is `TYPE(MESH_TYPE)`
2. **No cross-mesh access** in kernel (no `OMESH`, `MESHES()` array, `EXTERNAL_WALL`)
3. **RECURSIVE C wrapper** exists (or will be created) in `fds_c_interface.f90`
4. **C declaration** exists (or will be created) in `fds_fortran_interface.h`

## Patterns

### Pattern A: Pure Kernel (No Pre/Post Processing)

The sequential task is a trivial wrapper around one or more thread-safe kernels. No cross-mesh code.

**Examples**: VelocityCorrector, VelocityPredictor, DivPart2, CorrStep1, DensityPred

```
Orchestrator: collect N → emit N work tokens (no Fortran calls)
Kernel Task:  call fds_*_kernel(nm, ...) per mesh
Collector:    collect N → sort by nm → emit N MeshData
```

### Pattern B: Sequential Pre-Processing + Parallel Kernel

The sequential task contains cross-mesh code (reads `OMESH`) that must stay sequential, followed by a thread-safe kernel.

**Examples**: CorrDivPart1 (COMBUSTION_BC + kernel), DivSetup (VISCOSITY_BC + kernel)

```
Orchestrator: collect N → run cross-mesh code on each → emit N work tokens
Kernel Task:  call fds_*_kernel(nm, ...) per mesh
Collector:    collect N → sort by nm → emit N MeshData
```

The orchestrator calls the cross-mesh routines (e.g., `fds_combustion_bc`, `fds_viscosity_bc`) sequentially for each mesh before dispatching parallel kernel work.

## Step-by-Step Procedure

### Step 1: Analyze the Sequential Task

Read the task's `execute()` method to identify:
- Which Fortran routines it calls
- Which have thread-safe kernels (check `*_kernels.f90`)
- Which have cross-mesh dependencies (check for `OMESH`, `POINT_TO_MESH`)

```cpp
// Example: CorrDivPart1Task
void execute(std::shared_ptr<MeshData> data) override {
    fds_combustion_bc(data->nm);           // Cross-mesh (reads OMESH%Q) → orchestrator
    fds_divergence_part_1(data->t, data->dt, data->nm);  // Has kernel → parallel
    this->addResult(data);
}
```

### Step 2: Create C Kernel Wrapper (if needed)

If no RECURSIVE wrapper exists yet, add one in `fds_c_interface.f90`:

```fortran
RECURSIVE SUBROUTINE C_FDS_DIVERGENCE_PART_1_KERNEL(NM, T, DT) &
    BIND(C, NAME="fds_divergence_part_1_kernel")
    USE DIVG_KERNELS, ONLY: DIVERGENCE_PART_1_KERNEL
    USE MESH_VARIABLES, ONLY: MESHES
    INTEGER(C_INT), VALUE :: NM
    REAL(C_DOUBLE), VALUE :: T, DT
    CALL DIVERGENCE_PART_1_KERNEL(MESHES(NM), T, DT, NM)
END SUBROUTINE C_FDS_DIVERGENCE_PART_1_KERNEL
```

Key points:
- `RECURSIVE` keyword ensures thread-safety (multiple threads can enter simultaneously)
- `BIND(C, NAME="...")` provides the C-callable name
- Passes `MESHES(NM)` directly — NO `POINT_TO_MESH`
- For logical arguments, use a local variable: `EST_FLAG = (ESTIMATED /= 0)` then pass `EST_FLAG` (avoids Fortran keyword argument parsing issues)
- For array workspace (e.g., `GX(0:IBAR_MAX)`), allocate locally in the wrapper

Add the C declaration in `fds_fortran_interface.h`:

```cpp
void fds_divergence_part_1_kernel(int nm, double t, double dt);
```

### Step 3: Create Work Token

`Source/hedgehog/data/<name>_data.h`:

```cpp
struct <Name>Work {
    int nm;
    double t, dt;
    // Add other kernel parameters as needed
    std::shared_ptr<MeshData> originalMeshData;

    <Name>Work(int nm_, double t_, double dt_, std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), originalMeshData(md) {}
};
```

Include only the parameters the kernel needs. The `originalMeshData` pointer preserves the MeshData token for downstream routing.

### Step 4: Create Orchestrator and Collector States

`Source/hedgehog/state/<name>_state.h`:

**Orchestrator** (Pattern A — no pre-processing):
```cpp
class <Name>Orchestrator : public hh::AbstractState<1, MeshData, <Name>Work> {
    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                this->addResult(std::make_shared<<Name>Work>(md->nm, md->t, md->dt, md));
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};
```

**Orchestrator** (Pattern B — with sequential pre-processing):
```cpp
class <Name>Orchestrator : public hh::AbstractState<1, MeshData, <Name>Work> {
    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing (cross-mesh routines)
            for (auto &md : collected_) {
                fds_cross_mesh_routine(md->nm);  // Reads OMESH, must be sequential
            }
            // Dispatch parallel work
            for (auto &md : collected_) {
                this->addResult(std::make_shared<<Name>Work>(md->nm, md->t, md->dt, md));
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};
```

**Collector** (same for both patterns):
```cpp
class <Name>Collector : public hh::AbstractState<1, <Name>Work, MeshData> {
    void execute(std::shared_ptr<<Name>Work> work) override {
        results_.push_back(work);
        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) { return a->nm < b->nm; });
            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }
            results_.clear();
            results_.reserve(nmeshes_);
        }
    }
};
```

**Important**: The collector sorts results by `nm` to ensure deterministic downstream ordering regardless of which thread finishes first.

### Step 5: Create Kernel Task

`Source/hedgehog/task/<name>_kernel_task.h`:

```cpp
class <Name>KernelTask : public hh::AbstractTask<1, <Name>Work, <Name>Work> {
public:
    explicit <Name>KernelTask(size_t numThreads)
        : hh::AbstractTask<1, <Name>Work, <Name>Work>("<Name>Kernel", numThreads) {}

    void execute(std::shared_ptr<<Name>Work> work) override {
        fds_<name>_kernel(work->nm, work->t, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, <Name>Work, <Name>Work>> copy() override {
        return std::make_shared<<Name>KernelTask>(this->numberThreads());
    }
};
```

The `copy()` method is required when `numThreads > 1` — Hedgehog clones the task to create one instance per thread.

### Step 6: Wire into the Graph

In `Source/hedgehog/graph/fds_graph.h`:

**Add includes**:
```cpp
#include "../data/<name>_data.h"
#include "../state/<name>_state.h"
#include "../task/<name>_kernel_task.h"
```

**Remove** the sequential task variable (or comment it out with a NOTE).

**Create** sub-graph components:
```cpp
auto <name>OrchSM = std::make_shared<hh::StateManager<1, MeshData, <Name>Work>>(
    std::make_shared<<Name>Orchestrator>(nmeshes), "<Name>Orch");
auto <name>KernelTask = std::make_shared<<Name>KernelTask>(kernelThreads);
auto <name>CollectorSM = std::make_shared<hh::StateManager<1, <Name>Work, MeshData>>(
    std::make_shared<<Name>Collector>(nmeshes), "<Name>Collector");
```

**Replace** edges:
```cpp
// OLD:
// graph->edges(prevNode, sequentialTask);
// graph->edges(sequentialTask, nextNode);

// NEW:
graph->edges(prevNode, <name>OrchSM);
graph->edges(<name>OrchSM, <name>KernelTask);
graph->edges(<name>KernelTask, <name>CollectorSM);
graph->edges(<name>CollectorSM, nextNode);
```

### Step 7: Build and Test

```bash
cd build_hh
cmake --build . --target fds -j$(nproc)       # Standard FDS (no regression)
cmake --build . --target fds_hh -j$(nproc)    # Hedgehog FDS

# Test all cases
for dir in run_1mesh test_4mesh test_reac test_species; do
    cd ../test_cases/$dir
    rm -f *_devc.csv *_hrr.csv
    mpiexec --oversubscribe -n 1 ../../build_hh/Source/hedgehog/fds_hh *.fds
done

# Compare against baselines
diff run_1mesh/*_devc.csv archive/saved_results/hh_1mesh/*_devc.csv
```

All DEVC outputs must be byte-identical.

## Completed Sub-Graphs

| # | Name | Pattern | Kernel(s) | Pred/Corr |
|---|------|---------|-----------|-----------|
| 1 | VelocityCorrector | A | VELOCITY_CORRECTOR_KERNEL + CHECK_DIVERGENCE_KERNEL | Corr |
| 2 | VelocityPredictor | A | VELOCITY_PREDICTOR_KERNEL + CHECK_STABILITY_KERNEL | Pred |
| 3 | DivergencePart2 | A | DIVERGENCE_PART_2_KERNEL | Both (2 instances) |
| 4 | CorrStep1 | A | COMPUTE_VISCOSITY_KERNEL + MASS_FINITE_DIFFERENCES_NEW_KERNEL + DENSITY_KERNEL | Corr |
| 5 | DensityPred | A | DENSITY_KERNEL | Pred |
| 6 | CorrDivPart1 | B | DIVERGENCE_PART_1_KERNEL (seq: COMBUSTION_BC) | Corr |
| 7 | DivSetup | B | VELOCITY_FLUX_KERNEL (seq: VISCOSITY_BC) | Both (2 instances) |
| 8 | PredFinal | B | VELOCITY_BC_PROCESS_EDGES_KERNEL (seq: MATCH_VELOCITY, SYNTH_TURB, VELOCITY_BC_PREPROCESSING; post: CC_VELOCITY_BC) | Pred |
| 9 | CorrFinal | B | VELOCITY_BC_PROCESS_EDGES_KERNEL (seq: MATCH_VELOCITY, VELOCITY_BC_PREPROCESSING; post: CC_VELOCITY_BC, UPDATE_GLOBAL_OUTPUTS) | Corr |

## Files per Sub-Graph

Each sub-graph adds 3 files:
```
Source/hedgehog/
├── data/<name>_data.h          Work token struct
├── state/<name>_state.h        Orchestrator + Collector states
└── task/<name>_kernel_task.h   Parallel kernel task
```

And modifies 3 files:
```
Source/hedgehog/fds_c_interface.f90     RECURSIVE kernel wrapper (if new kernel)
Source/hedgehog/fds_fortran_interface.h C declaration (if new kernel)
Source/hedgehog/graph/fds_graph.h       Graph wiring
```

## Automation Checklist

For each sub-graph:

1. [ ] Identify sequential task and its Fortran calls
2. [ ] Classify: Pattern A (pure kernel) or Pattern B (sequential pre-processing + kernel)
3. [ ] Create RECURSIVE C wrapper if kernel not yet exposed
4. [ ] Add C declaration in fds_fortran_interface.h
5. [ ] Create work token data structure
6. [ ] Create orchestrator state (with sequential pre-processing for Pattern B)
7. [ ] Create collector state (sort by nm)
8. [ ] Create parallel kernel task (with copy() override)
9. [ ] Add includes to fds_graph.h
10. [ ] Create sub-graph components in buildFDSGraph
11. [ ] Replace sequential task edges with sub-graph edges
12. [ ] Build both fds and fds_hh
13. [ ] Test across all test cases (1/3/4/5-mesh)
14. [ ] Verify byte-identical DEVC output
