# FDS Hedgehog Parallelization Documentation

## Quick Start

1. Read [PARALLELIZATION_PROGRESS.md](PARALLELIZATION_PROGRESS.md) for current status
2. Review appropriate METHOD files for your task (see below)

## Core Methodology Files

### [METHOD_KERNEL_EXTRACTION.md](METHOD_KERNEL_EXTRACTION.md)
How to extract thread-safe kernels from Fortran modules.
- Creating `*_kernels.f90` modules with `TYPE(MESH_TYPE)` arguments
- Removing module-level state (POINT_TO_MESH)
- Patterns: index-based, pointer-based, local pointer alias shadowing

### [METHOD_PATTERN_B_COMPLEX.md](METHOD_PATTERN_B_COMPLEX.md)
Complex routines with cross-mesh dependencies (Pattern B).
- Three-phase: preprocessing → parallel kernel → finalization
- Dedicated sub-graph wrappers (`graph/<routine>_subgraph.h`)
- Example: WallBC, DensityBlock, VelocityFluxBlock

### [METHOD_MODULE_SPLIT.md](METHOD_MODULE_SPLIT.md)
Decomposing large Fortran modules (>5K lines) into sub-modules.

### [METHOD_DEPENDENCY_EXCHANGE.md](METHOD_DEPENDENCY_EXCHANGE.md)
Replacing global exchange barriers with per-mesh dependency tracking.
- Push-then-gate pattern (parallel push task + reusable gate state)
- DynBitset and MeshDependencyGraph infrastructure
- Thread safety analysis (push model avoids read-after-release races)
- Extending to other exchange codes and MPI

## Progress Tracking

### [PARALLELIZATION_PROGRESS.md](PARALLELIZATION_PROGRESS.md)
**Current status and overall progress.**
- 20 completed sub-graphs (8 with K-block decomposition)
- All phases complete (Phases 1-4)
- Remaining work: advanced optimization (hybrid MPI, relaxed barriers, NUMA)

## File Organization

```
docs/
├── README.md                           # This file
├── PARALLELIZATION_PROGRESS.md         # Current status (START HERE)
│
├── METHOD_KERNEL_EXTRACTION.md         # Kernel extraction patterns
├── METHOD_PATTERN_B_COMPLEX.md         # Pattern B (complex routines)
├── METHOD_MODULE_SPLIT.md              # Module decomposition
└── METHOD_DEPENDENCY_EXCHANGE.md      # Dependency-aware exchange
```

## Workflow for Parallelizing a New Routine

```
1. Identify target routine
   └─→ Check PARALLELIZATION_PROGRESS.md for remaining tasks

2. Determine pattern
   ├─→ Cross-mesh deps? → METHOD_PATTERN_B_COMPLEX.md (Pattern B)
   └─→ Exchange barrier? → METHOD_DEPENDENCY_EXCHANGE.md

3. Extract/convert kernels
   └─→ METHOD_KERNEL_EXTRACTION.md

4. Create sub-graph components
   └─→ Orchestrator, kernel task, collector

5. Test and verify
   └─→ python3 test_cases/run_tests.py -v
   └─→ python3 test_cases/run_verification.py test --no-redundant --max-gold-time 30 --timeout 120 --tolerance 1e-6
```
