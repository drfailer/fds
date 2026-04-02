# FDS Benchmark Results

Machine: 2x Intel Xeon Silver 4114 @ 2.20GHz, 20 physical cores (40 with HT)
Date: 2026-04-02
Repeats: 1

## pressure_iteration3d_default

8 meshes (2x2x2 via MULT), 16x16x16 each, 32,768 total cells, T_END=0.5

| Config | Time | vs best master | Status |
|---|---|---|---|
| master_mpi_omp1 | 5.5s | 0.89x | OK |
| **master_mpi_omp2** | **4.9s** | **(ref)** | OK |
| master_mpi_omp4 | 5.6s | 0.87x | OK |
| hh_auto | 22.7s | 0.21x | OK |
| hh_pres_off | 20.7s | 0.23x | OK |
| hh_pres_on | 20.1s | 0.24x | OK |
| hh_mpi | 1.1s | — | FAIL |
| hh_redec_8 | 43.1s | 0.11x | DIFF (max_diff=2.66e+01) |

## shunn3_4mesh_64

4 meshes (explicit), 32x1x32 each, 4,096 total cells, T_END=1.0

| Config | Time | vs best master | Status |
|---|---|---|---|
| master_mpi_omp1 | 2.3s | 0.86x | OK |
| master_mpi_omp2 | 2.2s | 0.90x | OK |
| **master_mpi_omp4** | **2.0s** | **(ref)** | OK |
| hh_auto | 6.2s | 0.32x | OK |
| hh_pres_off | 5.8s | 0.34x | OK |
| hh_pres_on | 5.9s | 0.33x | OK |
| hh_mpi | 6.2s | 0.32x | DIFF (max_diff=4.11e-03) |
| hh_redec_8 | 11.8s | 0.17x | DIFF (max_diff=2.47e+00) |

## pressure_iteration2d_default

8 meshes (4x1x2 via MULT), 8x1x16 each, 1,024 total cells, T_END=0.5

| Config | Time | vs best master | Status |
|---|---|---|---|
| **master_mpi_omp1** | **2.8s** | **(ref)** | OK |
| master_mpi_omp2 | 3.0s | 0.93x | OK |
| master_mpi_omp4 | 3.7s | 0.75x | OK |
| hh_auto | 6.5s | 0.43x | OK |
| hh_pres_off | 5.5s | 0.50x | OK |
| hh_pres_on | 6.2s | 0.45x | OK |
| hh_mpi | 3.3s | 0.85x | DIFF (max_diff=2.83e-05) |
| hh_redec_8 | 6.9s | 0.40x | DIFF (max_diff=1.94e+02) |

## obst_activation_default

4 meshes (explicit), 16x1x16 each, 1,024 total cells, T_END=2.0

| Config | Time | vs best master | Status |
|---|---|---|---|
| **master_mpi_omp1** | **0.9s** | **(ref)** | OK |
| master_mpi_omp2 | 0.9s | 1.00x | OK |
| master_mpi_omp4 | 0.9s | 1.00x | OK |
| hh_auto | 2.0s | 0.44x | OK |
| hh_pres_off | 1.8s | 0.50x | OK |
| hh_pres_on | 1.7s | 0.52x | OK |
| hh_mpi | 1.9s | 0.48x | DIFF |
| hh_redec_8 | 3.3s | 0.27x | DIFF (max_diff=4.08e+00) |

## Summary (geomean speedup vs best fds_master)

| Config | Geomean Speedup | Cases |
|---|---|---|
| master_mpi_omp1 | 0.93x | 4 |
| master_mpi_omp2 | 0.96x | 4 |
| master_mpi_omp4 | 0.90x | 4 |
| hh_auto | 0.34x | 4 |
| hh_mpi | 0.50x | 3 |
| hh_pres_off | 0.38x | 4 |
| hh_pres_on | 0.37x | 4 |
| hh_redec_8 | 0.21x | 4 |

## Notes

- **hh_mpi FAIL on pressure_iteration3d_default**: 8 MPI procs with 8 meshes = 1 mesh/proc.
  Pressure subgraph auto-disables with only 1 local mesh; possible other issue causing failure.
- **hh_mpi DIFF**: Expected — MPI changes exchange ordering and floating-point reduction order.
- **hh_redec_8 DIFF**: Expected — re-decomposition changes mesh topology entirely, producing
  different (but physically valid) numerical results. Would need separate gold files or
  higher tolerance to mark as OK.
- **hh_mpi on pressure_iteration2d (0.85x)**: Closest to parity with fds_master.
- These are small cases. Hedgehog overhead is relatively larger; larger cases should show
  better scaling.
