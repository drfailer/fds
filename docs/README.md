# FDS Hedgehog Parallelization Documentation

## Quick Start

1. Read [PARALLELIZATION_PROGRESS.md](PARALLELIZATION_PROGRESS.md) for current status
2. Review appropriate METHOD files for your task (see below)
3. Check completed examples: [WALL_BC_PARALLELIZATION_PLAN.md](WALL_BC_PARALLELIZATION_PLAN.md)

## Core Methodology Files

### [METHOD_KERNEL_EXTRACTION.md](METHOD_KERNEL_EXTRACTION.md)
How to extract thread-safe kernels from Fortran modules.
- Creating `*_kernels.f90` modules with `TYPE(MESH_TYPE)` arguments
- Removing module-level state (POINT_TO_MESH)
- Patterns: index-based, pointer-based, local pointer alias shadowing

### [METHOD_SUBGRAPH.md](METHOD_SUBGRAPH.md)
Converting sequential tasks into parallel sub-graphs (Pattern A).
- Pure parallel execution (no cross-mesh dependencies)
- Examples: VelocityCorrector, VelocityPredictor, DivPart2

### [METHOD_PATTERN_B_COMPLEX.md](METHOD_PATTERN_B_COMPLEX.md)
Complex routines with cross-mesh dependencies (Pattern B).
- Three-phase: preprocessing → parallel kernel → finalization
- Dedicated sub-graph wrappers (`graph/<routine>_subgraph.h`)
- Example: WallBC, DensityBlock, VelocityFluxBlock

### [METHOD_MODULE_SPLIT.md](METHOD_MODULE_SPLIT.md)
Decomposing large Fortran modules (>5K lines) into sub-modules.

### [METHOD_MESH_BLOCK.md](METHOD_MESH_BLOCK.md)
K-block decomposition for intra-mesh parallelism.
- Splitting mesh computation along K dimension
- Orchestrator/block-kernel/collector pattern
- K-safety analysis (ghost cells, face values, wall loops)

## Progress Tracking

### [PARALLELIZATION_PROGRESS.md](PARALLELIZATION_PROGRESS.md)
**Current status and overall progress.**
- 20 completed sub-graphs (8 with K-block decomposition)
- Phase 4 complete: DivergencePart2, Density block-decomposed; DivPart1 not viable
- Remaining sequential tasks and future work

## Pipelining Research

### [pipelining/README.md](pipelining/README.md)
Data-flow analysis for intra-timestep pipelining parallelism.
- Complete read/write dependency map for all predictor/corrector routines
- Two identified opportunities: predictor two-level pipeline, corrector major pipeline
- Dependency proof tables showing zero data conflicts
- Section-level cost analysis with ops/cell estimates
- Hedgehog graph recommendations (3 tiers by impact/complexity)

### [pipelining/IMPLEMENTATION_PROGRESS.md](pipelining/IMPLEMENTATION_PROGRESS.md)
Phased implementation plan for pipelining changes.
- 6 phases: kernel extraction → scratch arrays → driver → Fork 1 → Fork 2 → predictor
- Current phase tracking and test criteria

## Reference Implementations

### [WALL_BC_PARALLELIZATION_PLAN.md](WALL_BC_PARALLELIZATION_PLAN.md)
Complete Pattern B implementation reference (three-phase decomposition, test results).

### [CHANGE_TIMESTEP_REFACTORING.md](CHANGE_TIMESTEP_REFACTORING.md)
CFL retry loop refactoring into state-managed sub-graph.

### [BLOG_AI_ASSISTED_REWRITE.md](BLOG_AI_ASSISTED_REWRITE.md)
Retrospective on AI-assisted parallelization methodology.

## File Organization

```
docs/
├── README.md                           # This file
├── PARALLELIZATION_PROGRESS.md         # Current status (START HERE)
│
├── METHOD_KERNEL_EXTRACTION.md         # Kernel extraction patterns
├── METHOD_SUBGRAPH.md                  # Pattern A (pure kernel)
├── METHOD_PATTERN_B_COMPLEX.md         # Pattern B (complex routines)
├── METHOD_MODULE_SPLIT.md              # Module decomposition
├── METHOD_MESH_BLOCK.md                # K-block decomposition
│
├── WALL_BC_PARALLELIZATION_PLAN.md     # Reference: Pattern B example
├── CHANGE_TIMESTEP_REFACTORING.md      # Reference: cycle/retry pattern
├── BLOG_AI_ASSISTED_REWRITE.md         # Process retrospective
│
├── pipelining/                         # Pipelining parallelism research
│   ├── README.md                       # Analysis and implementation strategy
│   ├── fds_dataflow.dot / .svg         # Full data-flow dependency graph
│   └── fds_pipeline_opportunities.dot / .svg  # Identified opportunities
│
└── architecture/                       # Codebase analysis and diagrams
```

## Workflow for Parallelizing a New Routine

```
1. Identify target routine
   └─→ Check PARALLELIZATION_PROGRESS.md for remaining tasks

2. Determine pattern
   ├─→ Pure kernel? → METHOD_SUBGRAPH.md (Pattern A)
   ├─→ Cross-mesh deps? → METHOD_PATTERN_B_COMPLEX.md (Pattern B)
   └─→ K-block parallel? → METHOD_MESH_BLOCK.md

3. Extract/convert kernels
   └─→ METHOD_KERNEL_EXTRACTION.md

4. Create sub-graph components
   └─→ Orchestrator, kernel task, collector

5. Test and verify
   └─→ python3 test_cases/run_tests.py -v
   └─→ python3 test_cases/run_verification.py test --no-redundant --max-gold-time 30 --timeout 120 --tolerance 1e-6
```
