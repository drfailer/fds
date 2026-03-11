# Hedgehog Skill Revision

How the Claude Code Hedgehog skill was revised based on practical experience
from the FDS parallelization project.

## Context

The original Hedgehog skill (`~/.claude/skills/hedgehog/SKILL.md`) was written
as a generic API reference with matrix multiplication examples. After building
11 parallel sub-graphs for FDS across multiple sessions — encountering
deadlocks, thread-safety crashes, non-deterministic results, and cycle
termination bugs along the way — the skill was revised to capture all the
hard-won patterns and pitfalls.

## Sources Analyzed

The revision drew from three categories of project artifacts:

### 1. Implementation Files (Source/hedgehog/)

Every `.h` file in the Hedgehog integration was read to extract the actual
patterns used in production:

| Directory | Files Read | Patterns Extracted |
|-----------|-----------|-------------------|
| `data/` | `mesh_data.h`, `barrier_data.h`, `velocity_corrector_data.h`, + 10 work token files | Data token design, work token with `originalData` back-pointer, aggregate barrier token |
| `state/` | `collector_state.h`, `timestep_state.h`, `mesh_barrier_state.h`, `div_setup_state.h`, + 10 orchestrator/collector files | Scatter-gather pattern, collector sorting, cycle termination with `canTerminate()`, passthrough barrier, Pattern A vs Pattern B orchestrators |
| `task/` | `velocity_corrector_kernel_task.h`, `predictor_tasks.h`, `corrector_tasks.h`, + 9 kernel task files | Parallel task with `copy()`, thread count configuration, foreign function kernel wrapping |
| `graph/` | `fds_graph.h`, `change_timestep_subgraph.h` | Full graph wiring with 11 sub-graphs, sub-graph composition, cycle management, data-driven `canTerminate()` |
| Root | `main_hh.cpp`, `fds_c_interface.f90`, `fds_fortran_interface.h` | Graph lifecycle, token injection, Fortran-C++ bridge, `RECURSIVE` wrappers |

### 2. Documentation Files (docs/)

All project documentation was analyzed for lessons learned:

| Document | Key Insights Extracted |
|----------|----------------------|
| `PARALLELIZATION_PROGRESS.md` | Complete inventory of 11 sub-graphs, remaining sequential bottlenecks, Amdahl's law analysis (49% sequential = max 2x speedup), super-linear kernel speedup (9x from 4 threads) |
| `METHOD_SUBGRAPH.md` | Step-by-step sub-graph creation procedure, Pattern A (pure kernel) vs Pattern B (sequential pre-processing + parallel kernel), automation checklist |
| `METHOD_KERNEL_EXTRACTION.md` | Thread-safe kernel requirements, `M%` prefix transformation, `RECURSIVE` keyword, common pitfalls (missing `USE TYPES`, line length, CMake ordering) |
| `hedgehog_subgraph_methodology.md` | Orchestrator/collector design, work token rationale, testing methodology (sequential first, then parallel), troubleshooting guide |
| `subgraph_quick_start.md` | Concise implementation templates, pattern classification (A/B/C/D) |
| `PARALLEL_EXECUTION_FINDINGS.md` | Thread-safety investigation, OpenMP nested parallelism crashes, global module variable issues |
| `FINAL_SUMMARY.md` | Early prototype findings (before parallel execution worked), Fortran thread-safety blockers |
| `PARALLEL_TEST_SUMMARY.md` | Test results across 1-mesh and 4-mesh configurations |
| `THREAD_SAFE_FLAGS_TEST_REPORT.md` | `-frecursive -fno-automatic` flag testing, what works vs what doesn't |
| `RECURSIVE_KEYWORD_TEST_REPORT.md` | `RECURSIVE` keyword effects and limitations |

### 3. Git Commit History

The commit log was reviewed to trace the evolution of patterns:

```
fc4d622cb0  Add CC_IBM integration to parallelized sub-graphs
c63dcdced9  Extract 2 new kernels and add 4 sub-graphs
15e7d34bcb  Add methodology docs for kernel extraction, module split, sub-graph creation
c0d1af60fe  Add CorrDivPart1 and DivSetup sub-graphs
b48a4a45cf  Add density predictor sub-graph
1fafbd57fb  Fix array shape mismatches and sort mesh ordering at sync points
ba808ec249  Replace ChangeTimeStepTask while-loop with data-driven retry subgraph cycle
4e46ef0c92  Fix graph cycle termination so waitForTermination() returns cleanly
6ead0d349f  Add velocity predictor sub-graph
aabbcf5541  Remove OpenMP and enable parallel execution via selective kernel parallelization
```

Key commits that shaped the skill:
- `4e46ef0c92` — Cycle termination debugging led to Pattern 3 (canTerminate)
- `1fafbd57fb` — Collector sorting requirement discovered (deterministic ordering)
- `ba808ec249` — Sub-graph-in-cycle limitation discovered, leading to Pattern 4 caveats
- `aabbcf5541` — OpenMP removal, leading to the Foreign Function Kernels section

## What Changed

### Additions (not in original skill)

| Section | Why Added |
|---------|-----------|
| **Core Principles** | The original skill lacked any conceptual framework. The 5 principles distill the mental model needed to use Hedgehog correctly. |
| **Data Types** | Data token design was never covered. The `originalData` back-pointer pattern was discovered during the velocity corrector implementation and is essential for all scatter-gather sub-graphs. |
| **Pattern 1: Scatter-Gather** | The orchestrator -> kernel -> collector pattern was the primary pattern used across all 11 sub-graphs. Not mentioned in the original skill. |
| **Pattern 2: Barrier** | The N->1->N barrier pattern was used for all MESH_EXCHANGE, PRESSURE_ITERATION, and HVAC synchronization points. |
| **Pattern 3: Cycle Termination** | The original skill showed `canTerminate()` but didn't explain *why* it's needed (nodes in a cycle can't terminate normally), the lock/unlock requirement, or the dual-output-type technique. These were all discovered through debugging deadlocks. |
| **Pattern 4: Sub-Graphs** | Sub-graph composition was not mentioned. The critical limitation (internal cycle state resets) was discovered when the retry sub-graph hung inside the main time-stepping cycle. |
| **Pattern 5: Passthrough Barrier** | Pure synchronization without computation — needed to prevent concurrent Fortran execution of adjacent sequential tasks. |
| **Designing for Parallelism** | 5-step methodology from data decomposition to Amdahl's law profiling. Based on the full parallelization experience. |
| **Common Pitfalls** | 4 failure modes (deadlock, non-determinism, premature termination, poor speedup) with concrete symptoms, causes, and fixes. Each was encountered and debugged during the project. |
| **Foreign Function Kernels** | Fortran/C integration patterns (`RECURSIVE`, `-frecursive`, reentrant requirements). Essential for the FDS use case and any mixed-language project. |
| **Data-driven termination** | `canTerminate()` variant using monotonic data conditions (`t + dt >= tEnd`) instead of simple done flags. Used in the retry loop sub-graph. |

### Improvements to existing content

| Section | What Changed |
|---------|-------------|
| **Tasks** | Fixed constructor bug (original called `AbstractCUDATask` instead of `AbstractTask`). Added thread-safety rules, `copy()` explanation, multi-kernel-per-execute note. |
| **States** | Rewritten to emphasize orchestration-only role. Added StateManager wrapping (was mentioned but not shown). Added scatter-gather collect pattern. |
| **Graphs** | Added edge routing explanation (type-based). Added graph lifecycle (executeGraph, pushData, finishPushingData, waitForTermination, createDotFile). |
| **Code Organization** | Added guidance on file grouping (orchestrator + collector in same file, inline builder functions). |

### Removals

| What | Why Removed |
|------|-------------|
| Matrix multiplication examples dominating each section | Replaced with generic, reusable patterns. The matrix examples were specific to one use case and obscured the general API. |
| "Advanced Usage" section | Was a single vague sentence ("maximize computation capabilities"). Replaced with concrete Designing for Parallelism methodology. |
| `#define` macro pattern for template params | Modern C++ practice avoids macros for type aliases. The new examples use inline template parameters. |

## Methodology

The revision process:

1. **Read every implementation file** — Two research agents explored all
   `Source/hedgehog/` files and all `docs/` files in parallel, producing
   comprehensive inventories of patterns and techniques.

2. **Read the current skill** — Identified gaps between what the skill
   documented and what the project actually used.

3. **Cross-reference with commit history** — Traced which patterns emerged
   from debugging specific problems (cycle termination from deadlock fixes,
   collector sorting from non-determinism fixes, etc.).

4. **Synthesize generic patterns** — Abstracted FDS-specific implementations
   into domain-agnostic Hedgehog patterns that apply to any project.

5. **Write and install** — Produced the new skill at
   `~/.claude/skills/hedgehog/SKILL.md`, preserving the original structure
   where possible while adding all new content.

## Validation

The new skill was validated against the actual implementation:
- Every code pattern in the skill has a corresponding working implementation in
  `Source/hedgehog/`
- Every pitfall described was actually encountered and resolved during the project
- The 5-step parallelization methodology matches the actual steps taken across
  11 sub-graph implementations
- The cycle termination pattern matches the working `TimestepLoopStateManager`
  and `RetryLoopStateManager` implementations
