# FDS Architecture Documentation

This directory contains a comprehensive analysis of the FDS codebase architecture,
produced to support the Hedgehog integration project (per-mesh parallel execution).

## Contents

| File | Description |
|------|-------------|
| `FDS_ARCHITECTURE.md` | Main analysis document covering all five topics below |
| `graphs/module_dependencies.dot` | Fortran module USE-dependency graph (layered) |
| `graphs/execution_flow.dot` | Main time-stepping loop with predictor/corrector phases |
| `graphs/pressure_iteration.dot` | Pressure solver inner iteration loop detail |
| `graphs/data_flow.dot` | How key 3D/4D field arrays flow between routines |
| `graphs/mpi_communication.dot` | All MPI synchronization points in one time step |
| `graphs/thread_safety.dot` | Thread safety status of every per-mesh routine |

## Rendering the Graphs

All graphs are in Graphviz DOT format. Render with:

```bash
# Single file
dot -Tsvg graphs/module_dependencies.dot -o graphs/module_dependencies.svg

# All at once
for f in graphs/*.dot; do
    dot -Tsvg "$f" -o "${f%.dot}.svg"
done

# PNG alternative (higher DPI for large graphs)
dot -Tpng -Gdpi=150 graphs/execution_flow.dot -o graphs/execution_flow.png
```

Requires `graphviz` package (`sudo apt install graphviz` or `brew install graphviz`).

## Topics Covered

1. **Module Dependencies** - All ~50 Fortran modules, their USE relationships, and the
   orchestration/kernel split pattern used for thread-safe computation.

2. **Execution Flow** - The complete predictor-corrector time-stepping loop, including
   the pressure iteration inner loop and CFL-based time step adaptation.

3. **Data Types & Data Flow** - MESH_TYPE fields, their physical meaning, and how they
   are produced/consumed through one time step. Covers the chain:
   species transport -> EOS -> divergence -> Poisson -> velocity update.

4. **Intra-Node Parallelism** - Thread safety assessment for running multiple meshes
   concurrently within one MPI process. Identifies POINT_TO_MESH, SAVE variables,
   and global accumulators as critical blockers. Maps which routines are already
   safe (kernel modules) vs. which need refactoring.

5. **Inter-Node Parallelism** - All 18 MESH_EXCHANGE codes, their data payloads, when
   they occur, and whether they use persistent or blocking MPI. Also covers global
   reductions (MPI_ALLREDUCE) for pressure zones, CFL sync, and diagnostics.

## How to Use This

- **For understanding the codebase**: Start with `FDS_ARCHITECTURE.md` sections 1-3.
- **For Hedgehog integration**: Focus on sections 4-5 and the `thread_safety.dot` graph.
- **For kernel extraction work**: Check section 4.1 for which modules still need refactoring.
- **For debugging MPI issues**: See section 5 and `mpi_communication.dot`.
