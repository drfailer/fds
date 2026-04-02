# FDS Benchmark Suite

Compares **fds_master** (MPI + OpenMP) vs **fds_hh** (Hedgehog dataflow graph)
across multiple configurations.

## Benchmark Cases

| Case | Category | Meshes | Mesh Size (IJK) | Cells/Mesh | Total Cells | T_END | Pressure Eligible | Notes |
|---|---|---|---|---|---|---|---|---|
| pressure_iteration3d_default | Pressure_Solver | 8 (2x2x2 via MULT) | 16x16x16 | 4,096 | 32,768 | 0.5 | Yes | Largest 3D pressure case |
| shunn3_4mesh_64 | Scalar_Analytical_Solution | 4 (explicit) | 32x1x32 | 1,024 | 4,096 | 1.0 | Yes | MMS analytical solution |
| pressure_iteration2d_default | Pressure_Solver | 8 (4x1x2 via MULT) | 8x1x16 | 128 | 1,024 | 0.5 | Yes | 2D pressure iteration |
| obst_activation_default | Pressure_Solver | 4 (explicit) | 16x1x16 | 256 | 1,024 | 2.0 | Yes | Obstruction activation |

## Run Configurations

Each case is tested with the following configurations:

| Config | Executable | MPI Procs | OMP Threads | Pressure Subgraph | Mesh Re-decomposition |
|---|---|---|---|---|---|
| master_mpi_omp1 | fds_master | nmeshes | 1 | N/A | No |
| master_mpi_omp2 | fds_master | nmeshes | 2 | N/A | No |
| master_mpi_omp4 | fds_master | nmeshes | 4 | N/A | No |
| hh_auto | fds_hh | 1 | 1 | auto | No |
| hh_pres_off | fds_hh | 1 | 1 | off | No |
| hh_pres_on | fds_hh | 1 | 1 | on (forced) | No |
| hh_mpi | fds_hh | nmeshes | 1 | auto | No |
| hh_redec_8 | fds_hh | 1 | 1 | auto | --mesh-dim 8 8 8 |

### Mesh re-decomposition effect (--mesh-dim 8 8 8)

| Case | Original Meshes | Split per Mesh | Resulting Meshes |
|---|---|---|---|
| pressure_iteration3d_default | 8 (MULT) | 2x2x2 = 8 | 64 |
| shunn3_4mesh_64 | 4 | 4x1x4 = 16 | 64 |
| pressure_iteration2d_default | 8 (MULT) | 1x1x2 = 2 | 16 |
| obst_activation_default | 4 | 2x1x2 = 4 | 16 |

## Usage

```bash
# 1. Generate gold files (required once)
python3 run_benchmark.py generate-gold

# 2. Run benchmark
python3 run_benchmark.py run

# 3. View results
python3 run_benchmark.py report

# With custom OMP threads
python3 run_benchmark.py run --omp-threads 2 4 8

# With extra MPI config for fds_hh
python3 run_benchmark.py run --mpi-procs 2

# Skip mesh re-decomposition tests
python3 run_benchmark.py run --no-mesh-dim

# Run specific cases
python3 run_benchmark.py run --case shunn3_4mesh_64

# On SLURM cluster
sbatch benchmark_slurm.sh
```

## Notes

- **fds_master** always runs with MPI (nproc = nmeshes). Pure OpenMP single-process
  is too slow for multi-mesh cases.
- **Pressure subgraph** auto-enables when: single MPI process, FFT or ULMAT solver,
  no tunnel preconditioner, >= 2 local meshes.
- Gold files are generated with fds_master MPI (nproc=nmeshes, OMP=1) for a
  deterministic baseline.
- Timing files (`_cpu.csv`, `_steps.csv`, `_pressit.csv`) are excluded from
  correctness comparison.
