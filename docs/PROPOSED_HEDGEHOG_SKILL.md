---
name: hedgehog
description: Details how to use the Hedgehog C++20 dataflow graph library for task-based parallelism.
---

# Hedgehog

Hedgehog is a C++20 library for expressing algorithms as **dataflow graphs**
where **data tokens flow between nodes** (tasks, states, sub-graphs). It is
well-suited for algorithms that exploit **data decomposition** — splitting a
problem into independent blocks that can be processed in parallel by dedicated
threads.

## Core Principles

1. **Data flows, not control flow.** Nodes execute when data arrives, not when
   explicitly called. The graph topology defines execution order.
2. **Tasks do computation.** A task receives a data token, performs work, and
   emits result tokens. Tasks can be multi-threaded (cloned across N threads).
3. **States do orchestration.** A state manages data routing, synchronization,
   and cycle control. States are always single-threaded and thread-safe.
4. **Types define edges.** An edge exists between two nodes when the output type
   of one matches the input type of the other. Hedgehog routes data by type.
5. **Graphs compose.** A graph is itself a node — it can be embedded as a
   sub-graph inside a larger graph.

## API Reference

### Data Types

Data tokens are plain C++ structs or classes wrapped in `std::shared_ptr`.
Design data types to carry everything a downstream node needs to process them.

```cpp
/// Primary token flowing through the graph (one per decomposed block)
struct BlockData {
    int blockId;          // Which block this token represents
    double param1;        // Parameters needed by downstream tasks
    // ... fields the computation needs ...

    BlockData(int id, double p1) : blockId(id), param1(p1) {}
};

/// Work token for parallel kernel dispatch (carries block + context)
struct KernelWork {
    int blockId;
    double param1;
    std::shared_ptr<BlockData> originalData;  // Preserve for downstream routing

    KernelWork(int id, double p1, std::shared_ptr<BlockData> orig)
        : blockId(id), param1(p1), originalData(orig) {}
};

/// Aggregate token for barrier synchronization (collects N blocks)
struct BarrierData {
    std::vector<std::shared_ptr<BlockData>> blocks;
    bool done = false;    // Termination flag for cycle control
};
```

**Design rules:**
- Keep data tokens lightweight — they are copied/moved frequently.
- Include an `originalData` pointer in work tokens to preserve the upstream
  token for downstream routing after parallel processing.
- Use distinct C++ types for different stages of the pipeline. Hedgehog routes
  data by type, so `BlockData` and `KernelWork` are different edge types.

### Tasks

Tasks are the **computation kernels** of the graph. Each task receives data,
performs work, and emits results. Tasks can run on multiple threads.

```cpp
#include <hedgehog/hedgehog.h>

// Template: <NumInputTypes, InputType1, ..., OutputType1, ...>
class ComputeTask : public hh::AbstractTask<1, KernelWork, KernelWork> {
public:
    explicit ComputeTask(size_t numThreads)
        : hh::AbstractTask<1, KernelWork, KernelWork>("ComputeTask", numThreads) {}

    void execute(std::shared_ptr<KernelWork> work) override {
        // Perform computation on work->blockId
        do_computation(work->blockId, work->param1);

        // Emit result downstream
        this->addResult(work);
    }

    // REQUIRED when numThreads > 1: Hedgehog clones the task per thread
    std::shared_ptr<hh::AbstractTask<1, KernelWork, KernelWork>> copy() override {
        return std::make_shared<ComputeTask>(this->numberThreads());
    }
};
```

**Task rules:**
- Tasks are for **computation only**. Do not put orchestration logic,
  synchronization, or data routing decisions in tasks.
- The `execute()` method receives one token at a time. The task processes it
  independently of other tokens — this is what makes parallelism safe.
- When `numThreads > 1`, Hedgehog creates `numThreads` clones via `copy()`.
  Each clone runs on its own thread and pulls tokens from a shared queue.
- **Thread safety**: Each thread gets its own task instance. The `execute()`
  body must not access shared mutable state (global variables, static data,
  shared file handles). All data must come from the work token or be read-only.
- A task can call multiple kernels in sequence within one `execute()` — they
  all run on the same thread for the same block.
- A task can have multiple input types (one `execute()` override per type).

### States

States are the **orchestration nodes** of the graph. They manage data routing,
synchronization, and flow control. States are always single-threaded.

A state is wrapped in a `hh::StateManager` to become a graph node.

```cpp
#include <hedgehog/hedgehog.h>

// Template: <NumInputTypes, InputType1, ..., OutputType1, ...>
class OrchestratorState : public hh::AbstractState<1, BlockData, KernelWork> {
public:
    explicit OrchestratorState(int nblocks)
        : hh::AbstractState<1, BlockData, KernelWork>(),
          nblocks_(nblocks) {
        collected_.reserve(nblocks);
    }

    void execute(std::shared_ptr<BlockData> data) override {
        collected_.push_back(data);

        // Wait until ALL blocks have arrived before dispatching
        if (static_cast<int>(collected_.size()) == nblocks_) {
            // Sequential pre-processing (if needed)
            for (auto &bd : collected_) {
                sequential_setup(bd->blockId);
            }

            // Dispatch work tokens for parallel execution
            for (auto &bd : collected_) {
                this->addResult(std::make_shared<KernelWork>(
                    bd->blockId, bd->param1, bd));
            }

            collected_.clear();
            collected_.reserve(nblocks_);
        }
    }

private:
    int nblocks_;
    std::vector<std::shared_ptr<BlockData>> collected_;
};
```

**State rules:**
- States are for **orchestration only**. Do not put heavy computation in states
  — move it to tasks where it can be parallelized.
- States are single-threaded and mutex-protected by the StateManager. This
  makes them safe for collecting tokens and making routing decisions.
- Use states to implement scatter-gather patterns: collect N tokens, then emit
  N work tokens (or 1 aggregate token).
- A state can have multiple input types and multiple output types.

**Using states as graph nodes:**

States cannot be added directly to a graph. Wrap them in a `hh::StateManager`:

```cpp
auto orchestratorSM = std::make_shared<hh::StateManager<1, BlockData, KernelWork>>(
    std::make_shared<OrchestratorState>(nblocks), "Orchestrator");
```

### Graphs

Graphs define the dataflow topology. Like tasks and states, graphs have typed
inputs and outputs. Graphs can be nested as sub-graphs.

```cpp
#include <hedgehog/hedgehog.h>

// Template: <NumInputTypes, InputType1, ..., OutputType1, ...>
using GraphType = hh::Graph<1, BlockData, BlockData>;
auto graph = std::make_shared<GraphType>("My Graph");

// Create nodes
auto orchestratorSM = std::make_shared<hh::StateManager<...>>(...);
auto computeTask = std::make_shared<ComputeTask>(numThreads);
auto collectorSM = std::make_shared<hh::StateManager<...>>(...);

// Define topology
graph->inputs(orchestratorSM);                    // Graph entry point
graph->edges(orchestratorSM, computeTask);        // Orchestrator -> Task
graph->edges(computeTask, collectorSM);           // Task -> Collector
graph->outputs(collectorSM);                      // Graph exit point
```

**Edge routing:** `graph->edges(A, B)` creates edges for **all output types of
A that match input types of B**. Hedgehog routes data by C++ type — if A emits
`KernelWork` and B accepts `KernelWork`, the edge is created automatically.

**Graph lifecycle:**

```cpp
graph->executeGraph();          // Spawn threads, start processing
graph->pushData(token1);        // Push initial data into graph inputs
graph->pushData(token2);
graph->finishPushingData();     // Signal no more external data
graph->waitForTermination();    // Block until graph completes
graph->createDotFile("stats.dot",
    hh::ColorScheme::EXECUTION,
    hh::StructureOptions::QUEUE);  // Export statistics
```

## Key Patterns

### Pattern 1: Scatter-Gather (Orchestrator -> Parallel Task -> Collector)

The most important pattern for data-parallel execution. Converts N sequential
operations into N parallel operations.

```
[Orchestrator State]  collects N tokens, dispatches N work tokens
        |
[Parallel Task]       numThreads threads process work concurrently
        |
[Collector State]     gathers N results, sorts, emits N tokens downstream
```

**Why it works:** The orchestrator is a synchronization point — it waits for
all N tokens before dispatching. This ensures all blocks are ready before
parallel execution begins. The collector is another synchronization point — it
waits for all N results before emitting downstream, ensuring the next stage
sees a complete set.

**Orchestrator** (collect all, optionally pre-process, dispatch):
```cpp
class ScatterState : public hh::AbstractState<1, BlockData, KernelWork> {
    void execute(std::shared_ptr<BlockData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nblocks_) {
            // Optional: sequential pre-processing on each block
            for (auto &bd : collected_) {
                sequential_setup(bd->blockId);  // e.g., boundary conditions
            }
            // Dispatch parallel work
            for (auto &bd : collected_) {
                this->addResult(std::make_shared<KernelWork>(..., bd));
            }
            collected_.clear();
            collected_.reserve(nblocks_);
        }
    }
};
```

**Collector** (gather all, sort for determinism, emit):
```cpp
class GatherState : public hh::AbstractState<1, KernelWork, BlockData> {
    void execute(std::shared_ptr<KernelWork> work) override {
        results_.push_back(work);
        if (static_cast<int>(results_.size()) == nblocks_) {
            // Sort by block index for deterministic downstream ordering
            std::sort(results_.begin(), results_.end(),
                [](const auto &a, const auto &b) {
                    return a->blockId < b->blockId;
                });
            // Optional: sequential post-processing
            for (auto &w : results_) {
                sequential_finalize(w->blockId);
            }
            // Emit original tokens downstream
            for (auto &w : results_) {
                this->addResult(w->originalData);
            }
            results_.clear();
            results_.reserve(nblocks_);
        }
    }
};
```

**Critical: Sort results in the collector.** When N threads process N blocks,
completion order is non-deterministic. Sorting by block index ensures
downstream nodes always see tokens in the same order, which is essential for
reproducible results.

### Pattern 2: Barrier (Collect N -> Process -> Scatter N)

For global operations that need all blocks before proceeding (e.g., global
reductions, inter-block communication, I/O).

```
[Collector State]    collects N BlockData -> emits 1 BarrierData
        |
[Barrier Task]       processes BarrierData (single-threaded) -> emits N BlockData
```

```cpp
/// Collector: N tokens -> 1 aggregate
class BarrierCollector : public hh::AbstractState<1, BlockData, BarrierData> {
    void execute(std::shared_ptr<BlockData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nblocks_) {
            std::sort(collected_.begin(), collected_.end(),
                [](const auto &a, const auto &b) { return a->blockId < b->blockId; });
            auto bd = std::make_shared<BarrierData>();
            bd->blocks = std::move(collected_);
            collected_ = {};
            collected_.reserve(nblocks_);
            this->addResult(bd);
        }
    }
};

/// Barrier task: 1 aggregate -> N tokens (single-threaded, numThreads=1)
class BarrierTask : public hh::AbstractTask<1, BarrierData, BlockData> {
    void execute(std::shared_ptr<BarrierData> data) override {
        global_operation(data);  // e.g., inter-block exchange
        for (auto &bd : data->blocks) {
            this->addResult(bd);
        }
    }
};
```

### Pattern 3: Cycle Termination with canTerminate()

When a graph has cycles (e.g., iterative time-stepping), Hedgehog cannot
determine termination automatically — nodes in a cycle wait on each other
indefinitely. You must create a **custom StateManager** that overrides
`canTerminate()` to break the cycle.

**This is the most subtle and error-prone part of Hedgehog.** Get it wrong and
the graph either hangs forever or terminates prematurely.

```cpp
/// State that manages the cycle decision
class LoopState : public hh::AbstractState<1, BarrierData, BlockData, BarrierData> {
public:
    void execute(std::shared_ptr<BarrierData> data) override {
        if (data->done) {
            done_ = true;
            // Emit BarrierData on the TERMINATION output (different type)
            this->addResult(data);
            return;
        }
        // Continue: emit BlockData tokens back into the cycle
        for (auto &bd : data->blocks) {
            this->addResult(bd);
        }
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    bool done_ = false;
};

/// Custom StateManager that overrides canTerminate()
class LoopStateManager
    : public hh::StateManager<1, BarrierData, BlockData, BarrierData> {
public:
    LoopStateManager(std::shared_ptr<LoopState> const &state,
                     std::string const &name)
        : hh::StateManager<1, BarrierData, BlockData, BarrierData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        // CRITICAL: Lock the state before accessing it.
        // canTerminate() is called from a different thread than execute().
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<LoopState>(
            this->state())->isDone();
        this->state()->unlock();
        return ret;
    }
};
```

**How it works:**

1. By default, a node can terminate when all its predecessors are done and its
   input queue is empty. In a cycle, predecessors are never "done" because they
   wait on the cycle node itself.
2. `canTerminate()` overrides this check. When it returns `true`, the node
   signals that it is ready to terminate even though predecessors are alive.
3. **Always lock/unlock** the state inside `canTerminate()`. The method is
   called from Hedgehog's internal thread, not the state's execution thread.
4. The cycle-breaking state typically has **two output types**: one for
   continuing the cycle (`BlockData` -> back to computation) and one for
   termination (`BarrierData` -> graph output). This type-based routing
   ensures termination data exits the cycle cleanly.

**Graph wiring for cycles:**

```cpp
// Cycle edge: LoopState -> back to beginning of pipeline
graph->edges(loopSM, firstTaskInCycle);

// Termination edge: LoopState -> graph output (different type)
graph->edges(loopSM, terminationSink);
graph->outputs(terminationSink);
```

**Common mistakes:**
- Forgetting to lock/unlock the state in `canTerminate()` -> race conditions.
- Using the same output type for both cycle and termination -> Hedgehog cannot
  distinguish which edge to route on (both edges get the data).
- Not setting `done_ = true` before emitting termination data -> canTerminate()
  returns false and the graph hangs.

### Pattern 4: Sub-Graphs (Nested Graphs)

Graphs can be used as nodes inside other graphs. This is useful for reusable
pipeline fragments or for cycles that are embedded in a larger linear pipeline.

```cpp
inline auto buildSubGraph(/* params */) {
    using SubGraphType = hh::Graph<1, BarrierData, BlockData>;
    auto subgraph = std::make_shared<SubGraphType>("SubGraph");

    auto entryTask = std::make_shared<EntryTask>();
    auto loopSM = std::make_shared<LoopStateManager>(...);
    auto processTask = std::make_shared<ProcessTask>();

    subgraph->inputs(entryTask);
    subgraph->edges(entryTask, processTask);
    subgraph->edges(processTask, loopSM);
    subgraph->edges(loopSM, processTask);  // Internal cycle
    subgraph->outputs(loopSM);             // Exit via different type

    return subgraph;
}

// Use in parent graph:
auto subgraph = buildSubGraph();
parentGraph->edges(upstreamNode, subgraph);
parentGraph->edges(subgraph, downstreamNode);
```

**Important limitation:** A sub-graph with an internal cycle resets its node
states between invocations when embedded in an outer cycle. This means:
- Internal nodes terminate after each invocation of the sub-graph.
- The `canTerminate()` state must be re-evaluated fresh each time.
- If the sub-graph's cycle state retains a `done_ = true` from a previous
  invocation, it will terminate immediately on the next invocation.

**When sub-graphs with cycles don't work in outer cycles**, replace the
sub-graph's internal cycle with a **while loop inside a single task**:

```cpp
class RetryTask : public hh::AbstractTask<1, BarrierData, BlockData> {
    void execute(std::shared_ptr<BarrierData> data) override {
        while (true) {
            process(data);
            if (!needsRetry(data)) break;
            adjustParameters(data);
        }
        for (auto &bd : data->blocks) {
            this->addResult(bd);
        }
    }
};
```

### Pattern 5: Passthrough Barrier (Synchronization Without Computation)

When adjacent pipeline stages must not run concurrently (e.g., because they
access shared state like global pointers), insert a passthrough barrier that
simply collects all tokens and re-emits them:

```cpp
class PassthroughBarrier : public hh::AbstractState<1, BlockData, BlockData> {
    void execute(std::shared_ptr<BlockData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nblocks_) {
            std::sort(collected_.begin(), collected_.end(),
                [](const auto &a, const auto &b) { return a->blockId < b->blockId; });
            for (auto &bd : collected_) {
                this->addResult(bd);
            }
            collected_.clear();
        }
    }
};
```

This ensures all N tokens from the upstream stage have completed before any
token enters the downstream stage.

## Designing for Parallelism

### Step 1: Identify the Data Decomposition

The first question: **what is the independent unit of work?**

- Matrix algorithms: matrix blocks (tiles)
- Mesh-based simulations: individual meshes or mesh partitions
- Image processing: image tiles
- Monte Carlo: independent samples

Each independent unit becomes a **data token**. The number of tokens determines
the available parallelism.

### Step 2: Separate Orchestration from Computation

For each operation in your algorithm, ask: **can this run independently on each
block, or does it need information from other blocks?**

- **Independent (block-local)** -> put in a parallel **Task** (numThreads > 1)
- **Cross-block (needs neighbor data)** -> put in a sequential **State** or
  single-threaded **Task** (numThreads = 1)

Common cross-block operations that must stay sequential:
- Boundary data exchange between blocks
- Global reductions (min, max, sum across blocks)
- Operations that read neighbor block data (ghost cells, halos)
- I/O operations (file writes, stdout)

### Step 3: Extract Thread-Safe Kernels

A computation kernel is thread-safe when it:
- Takes all data as explicit arguments (no global mutable state)
- Operates only on its own block's data
- Does not read or write other blocks' data
- Does not use thread-unsafe library state (shared FFT plans, etc.)
- Does not perform I/O

If a routine mixes block-local computation with cross-block operations, split
it: move the cross-block part to the orchestrator state, keep the block-local
part in the kernel task.

### Step 4: Choose the Right Number of Threads

```cpp
auto task = std::make_shared<ComputeTask>(numThreads);
```

- `numThreads = 1`: Sequential execution. Use for correctness verification.
- `numThreads = N` (N = number of blocks): Full parallelism. Each block gets
  its own thread.
- Only **kernel tasks** get `numThreads > 1`. All states, orchestrators, and
  collectors are always single-threaded (enforced by Hedgehog for states;
  should be set to 1 for non-kernel tasks).

### Step 5: Profile with Amdahl's Law

After parallelizing, measure the sequential fraction. If S% of execution time
is sequential (barriers, orchestrators, cross-block operations), the maximum
speedup is:

```
max_speedup = 1 / (S + (1-S)/N)
```

With 50% sequential time, even infinite threads give only 2x speedup. To
improve further, you must reduce the sequential fraction by decomposing more
operations into parallel kernels.

Use `graph->createDotFile(...)` to get per-node execution statistics and
identify bottlenecks.

## Code Organization

```
project/
  data/          Data token structs (BlockData, KernelWork, BarrierData)
  task/          Computation tasks (parallel kernels)
  state/         Orchestration states and state managers
  graph/         Graph construction functions
  tool/          Helper functions and utilities
```

- One header file per data type, task, or state group.
- Group related orchestrator + collector states in the same file.
- Use inline builder functions (`buildGraph(...)`, `buildSubGraph(...)`)
  in graph headers.

## Common Pitfalls

### Graph Hangs (Deadlock)

**Symptom:** `waitForTermination()` never returns.

**Causes:**
1. Missing `canTerminate()` override on a cycle node.
2. Collector expects N tokens but fewer arrive (token lost or filtered).
3. Task's `execute()` doesn't call `addResult()` on some code path.
4. Cycle state emits on the cycle-continue type even when done.

**Fix:** Add debug prints in states to trace token counts. Verify N matches
the actual number of tokens in flight.

### Non-Deterministic Results

**Symptom:** Results change between runs with `numThreads > 1`.

**Causes:**
1. Collector doesn't sort results before emitting -> downstream sees tokens
   in random order.
2. Kernel accesses shared mutable state (global variables, static data).
3. Floating-point reduction order depends on thread completion order.

**Fix:** Always sort in collectors. Ensure kernels are pure functions of their
arguments. Use indexed writes (`result[blockId] = value`) instead of
accumulations.

### Premature Termination

**Symptom:** Graph terminates before processing all data.

**Causes:**
1. `canTerminate()` returns true too early (state not properly tracking
   completion).
2. Forgetting to lock/unlock state in `canTerminate()`.
3. Sub-graph internal state carries over `done_ = true` from previous
   invocation.

**Fix:** Always lock the state in `canTerminate()`. Reset state between
invocations if the sub-graph is called multiple times.

### Poor Speedup

**Symptom:** Parallel execution is barely faster than sequential.

**Causes:**
1. Sequential fraction too high (Amdahl's law).
2. Kernel is too lightweight — thread management overhead dominates.
3. False sharing in cache lines between threads.
4. Orchestrator pre-processing serializes most of the work.

**Fix:** Profile with dot file statistics. Move more computation into parallel
kernels. Combine multiple lightweight kernels into one task.

## Advanced: Foreign Function Kernels

When computation kernels are written in another language (e.g., Fortran, C),
wrap them for thread-safe calling:

```cpp
// C++ header declaring the foreign function
extern "C" {
    void compute_kernel(int block_id, double t, double dt);
}

// Task calling the foreign kernel
class ForeignKernelTask : public hh::AbstractTask<1, KernelWork, KernelWork> {
    void execute(std::shared_ptr<KernelWork> work) override {
        compute_kernel(work->blockId, work->param1, work->param2);
        this->addResult(work);
    }
    // ... copy() override ...
};
```

**Thread-safety requirements for foreign kernels:**
- The foreign function must not use global mutable state.
- For Fortran: use `RECURSIVE` keyword on the subroutine, pass all data
  through explicit arguments (not module-level pointer aliases), compile with
  `-frecursive` to make local variables stack-allocated.
- For C: avoid `static` local variables, use thread-local storage if needed.
- The foreign function must be reentrant: multiple threads calling it
  simultaneously with different arguments must produce correct results.

## Source Code

- [Hedgehog API](./ressources/hedgehog-api/)
- [Hedgehog Tutorials](./ressources/hedgehog-Tutorials/)
- [Communicator task implementation](./ressources/communicator_task/)
