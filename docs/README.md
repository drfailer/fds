# FDS Hedgehog Parallelization Documentation

This directory contains documentation for the FDS Hedgehog parallelization project.

## Quick Start

**New to the project?** Start here:
1. Read [PARALLELIZATION_PROGRESS.md](PARALLELIZATION_PROGRESS.md) for current status
2. Review appropriate METHOD files for your task (see below)
3. Check completed examples: [WALL_BC_PARALLELIZATION_PLAN.md](WALL_BC_PARALLELIZATION_PLAN.md)

## Core Methodology Files

Step-by-step procedures for parallelizing FDS routines:

### [METHOD_KERNEL_EXTRACTION.md](METHOD_KERNEL_EXTRACTION.md)
How to extract thread-safe kernels from Fortran modules.

**Use when:**
- Creating a new `*_kernels.f90` module
- Converting routines to use `TYPE(MESH_TYPE)` arguments
- Removing module-level state (POINT_TO_MESH)

**Key patterns:**
- Index-based access (small routines)
- Pointer-based access (large routines)
- RECURSIVE keyword for thread safety

### [METHOD_SUBGRAPH.md](METHOD_SUBGRAPH.md)
Converting sequential tasks into parallel sub-graphs (Pattern A).

**Use when:**
- Task is a simple wrapper around thread-safe kernels
- No cross-mesh dependencies (no OMESH access)
- Pure parallel execution possible

**Examples:** VelocityCorrector, VelocityPredictor, DivPart2

### [METHOD_PATTERN_B_COMPLEX.md](METHOD_PATTERN_B_COMPLEX.md)
Complex routines with cross-mesh dependencies (Pattern B).

**Use when:**
- Routine has 150+ lines of orchestration
- Contains OMESH reads/writes
- Only portion is parallelizable (80-90%)
- Multiple phases with different dependencies

**Example:** WallBC (3-phase: preprocessing → parallel kernel → finalization)

**Key concepts:**
- Sequential preprocessing (OMESH reads)
- Parallel kernel execution (local mesh only)
- Sequential finalization (OMESH writes)
- Global parameter computation
- Flag-based cell filtering
- **Dedicated sub-graph wrapper** (`graph/<routine>_subgraph.h`) for traceability

### [METHOD_MODULE_SPLIT.md](METHOD_MODULE_SPLIT.md)
Decomposing large Fortran modules into focused sub-modules.

**Use when:**
- Module exceeds 5000 lines
- Contains multiple independent functional areas
- Before extracting kernels from a large module

## Progress Tracking

### [PARALLELIZATION_PROGRESS.md](PARALLELIZATION_PROGRESS.md)
**Current status and overall progress.**

- Completed sub-graphs (12 total)
- Performance profiling results
- Remaining sequential tasks
- Future work roadmap

**Check this first** to understand what's been done and what's next.

## Reference Implementations

### [WALL_BC_PARALLELIZATION_PLAN.md](WALL_BC_PARALLELIZATION_PLAN.md)
**Complete Pattern B complex routine implementation.**

- Thread-safe callee conversions (4 routines, 1000+ lines)
- Three-phase decomposition (preprocessing, kernel, finalization)
- Hedgehog sub-graph integration
- Test results (byte-identical on 1-5 meshes)

**Best reference for Pattern B implementations.**

## Other Documentation

### [CHANGE_TIMESTEP_REFACTORING.md](CHANGE_TIMESTEP_REFACTORING.md)
Notes on change timestep sub-graph refactoring.

### [PROPOSED_HEDGEHOG_SKILL.md](PROPOSED_HEDGEHOG_SKILL.md)
Proposed Claude skill for Hedgehog parallelization.

### [HEDGEHOG_SKILL_REVISION.md](HEDGEHOG_SKILL_REVISION.md)
Revisions to the Hedgehog skill.

### [BLOG_AI_ASSISTED_REWRITE.md](BLOG_AI_ASSISTED_REWRITE.md)
Blog post about AI-assisted code rewriting.

## Test Results

Test reports are in `../test_cases/`:
- [WALLBC_TEST_REPORT.md](../test_cases/WALLBC_TEST_REPORT.md) - WallBC verification (5 test cases)
- [VELOCITY_CORRECTOR_TEST_REPORT.md](../test_cases/VELOCITY_CORRECTOR_TEST_REPORT.md) - Velocity corrector verification

## Workflow for Parallelizing a New Routine

```
1. Identify target routine
   └─→ Check PARALLELIZATION_PROGRESS.md for remaining tasks

2. Determine pattern
   ├─→ Pure kernel? → METHOD_SUBGRAPH.md (Pattern A)
   └─→ Cross-mesh dependencies? → METHOD_PATTERN_B_COMPLEX.md (Pattern B)

3. Extract/convert kernels
   └─→ METHOD_KERNEL_EXTRACTION.md

4. Create sub-graph components
   └─→ Orchestrator, kernel task, collector
   └─→ See METHOD_SUBGRAPH.md or METHOD_PATTERN_B_COMPLEX.md

5. Test and verify
   └─→ Run test suite: python3 test_cases/run_tests.py -v
   └─→ Verify byte-identical results

6. Document
   └─→ Update PARALLELIZATION_PROGRESS.md
   └─→ Create test report (if complex routine)
```

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
│
├── WALL_BC_PARALLELIZATION_PLAN.md     # Reference: Pattern B example
│
└── [Other docs...]                     # Misc notes and drafts
```

## Key Terminology

**Pattern A**: Pure kernel sub-graph (no preprocessing, just parallel execution)
- Example: VelocityCorrector, DivPart2

**Pattern B**: Sequential pre/post + parallel kernel
- Example: WallBC, CorrDivPart1, PredStep1

**OMESH**: Neighboring mesh data (cross-mesh dependencies)
- Reading OMESH → sequential preprocessing
- Writing OMESH → sequential finalization

**Thread-safe kernel**: Fortran routine marked RECURSIVE with TYPE(MESH_TYPE) argument
- No module-level state access
- No OMESH access
- Can execute in parallel across meshes

**Sub-graph**: Hedgehog graph component replacing a sequential task
- Orchestrator (collect N meshes)
- Kernel task (parallel execution)
- Collector (gather N results)

## Getting Help

1. **Not sure which pattern to use?**
   - Read the routine's code looking for OMESH, POINT_TO_MESH
   - OMESH = Pattern B, otherwise Pattern A

2. **Kernel extraction failing?**
   - Check METHOD_KERNEL_EXTRACTION.md for conversion patterns
   - See WALL_BC_PARALLELIZATION_PLAN.md for complex examples

3. **Tests not byte-identical?**
   - Check collector sorts results by mesh index
   - Verify all OMESH access moved to preprocessing/finalization
   - Review callee thread safety

4. **Performance not improving?**
   - Profile with Hedgehog dot file
   - Check sequential fraction (should be < 20% for good speedup)
   - Verify kernel is doing bulk of work (80-90%)
