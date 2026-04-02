#!/usr/bin/env python3
"""
FDS Benchmark Suite

Compares fds_master (MPI+OpenMP) vs fds_hh (Hedgehog) across multiple configurations:
  - fds_master with MPI (nproc=nmeshes) + varying OMP_NUM_THREADS
  - fds_hh single process (with and without pressure subgraph)
  - fds_hh with MPI (nproc=nmeshes, same as fds_master)
  - fds_hh with mesh re-decomposition (--mesh-dim 8 8 8)

Usage:
    # Run full benchmark (uses default cases)
    python3 run_benchmark.py run

    # Run with specific OMP thread counts
    python3 run_benchmark.py run --omp-threads 4 8 16 20

    # Run only specific cases
    python3 run_benchmark.py run --case shunn3_4mesh_64 pressure_iteration2d_default

    # Generate gold files first (required before benchmarking)
    python3 run_benchmark.py generate-gold --timeout 1800

    # Show results from a previous run
    python3 run_benchmark.py report

    # Verify correctness of all configurations against gold
    python3 run_benchmark.py verify
"""

import argparse
import copy
import json
import os
import re
import select
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Paths ────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DIR = Path(__file__).resolve().parent
VERIFICATION_DIR = REPO_ROOT / "Verification"
BENCHMARK_DIR = TEST_DIR / "benchmark"
BENCHMARK_RESULTS = BENCHMARK_DIR / "results.json"
COMPARE_SCRIPT = TEST_DIR / "compare_csv.py"

# Binaries
FDS_HH = REPO_ROOT / "build_hh" / "Source" / "hedgehog" / "fds_hh"
FDS_MASTER_PATH = shutil.which('fds_master')
FDS_MASTER = Path(FDS_MASTER_PATH) if FDS_MASTER_PATH else (REPO_ROOT.parent / "fds-master" / "build" / "fds")


# ── Benchmark Case Definitions ──────────────────────────────────────────────

@dataclass
class BenchmarkCase:
    """A benchmark test case."""
    name: str
    category: str
    filename: str         # FDS input filename
    chid: str             # CHID from &HEAD
    nmeshes: int          # Number of meshes in the input file
    total_cells: int      # Total grid cells (for reporting)
    t_end: float
    timeout: int = 1800   # Per-run timeout in seconds
    tolerance: float = 1e-6
    pressure_eligible: bool = False  # Can use pressure subgraph
    mesh_dim: Optional[Tuple[int, int, int]] = None  # For re-decomposition tests


# Cases selected for benchmarking:
# - Multi-mesh with FFT solver (pressure subgraph eligible)
# - No tunnel preconditioner, no CC_IBM
# - Large enough for meaningful timing
BENCHMARK_CASES = [
    BenchmarkCase(
        name="pressure_iteration3d_default",
        category="Pressure_Solver",
        filename="pressure_iteration3d_default.fds",
        chid="pressure_iteration3d_default",
        nmeshes=8,
        total_cells=8 * 16 * 16 * 16,  # 32768 — largest 3D pressure-iteration case
        t_end=0.5,
        timeout=1800,
        pressure_eligible=True,
    ),
    BenchmarkCase(
        name="shunn3_4mesh_64",
        category="Scalar_Analytical_Solution",
        filename="shunn3_4mesh_64.fds",
        chid="shunn3_4mesh_64",
        nmeshes=4,
        total_cells=4 * 32 * 1 * 32,  # 4096 (4 meshes of 32x1x32)
        t_end=1.0,
        timeout=600,
        tolerance=1e-4,  # _mms.csv diffs at ~1.5e-5 from -frecursive FP reordering
        pressure_eligible=True,
    ),
    BenchmarkCase(
        name="pressure_iteration2d_default",
        category="Pressure_Solver",
        filename="pressure_iteration2d_default.fds",
        chid="pressure_iteration2d_default",
        nmeshes=8,
        total_cells=8 * 8 * 1 * 16,  # 1024 (small cells, pressure-iteration heavy)
        t_end=0.5,
        timeout=600,
        pressure_eligible=True,
    ),
    BenchmarkCase(
        name="obst_activation_default",
        category="Pressure_Solver",
        filename="obst_activation_default.fds",
        chid="obst_activation_default",
        nmeshes=4,
        total_cells=4 * 16 * 1 * 16,  # 1024
        t_end=2.0,
        timeout=300,
        pressure_eligible=True,
    ),
    # hallways: 5 meshes, 73K cells, T_END=60 — very long runtime (~30min),
    # uncomment for extended benchmarks or on cluster with more resources
    # BenchmarkCase(
    #     name="hallways",
    #     category="Pressure_Solver",
    #     filename="hallways.fds",
    #     chid="hallways",
    #     nmeshes=5,
    #     total_cells=73728,
    #     t_end=60.0,
    #     timeout=3600,
    #     pressure_eligible=True,
    # ),
]


# ── Run Configurations ──────────────────────────────────────────────────────

@dataclass
class RunConfig:
    """A specific benchmark run configuration."""
    label: str
    exe: str              # 'fds_master' or 'fds_hh'
    omp_threads: int = 1
    mpi_procs: int = 1
    pressure_subgraph: str = "auto"  # "auto", "on", "off"
    mesh_dim: Optional[Tuple[int, int, int]] = None
    description: str = ""


def build_run_configs(omp_thread_counts: List[int],
                      enable_mesh_dim: bool = True,
                      mpi_procs: int = 0) -> List[RunConfig]:
    """Build list of run configurations for benchmarking."""
    configs = []

    # 1. fds_master with MPI + various OMP thread counts
    #    Always use MPI with nproc=nmeshes (pure OMP is too slow for multi-mesh)
    #    OMP=1 is always included as baseline
    omp_set = sorted(set([1] + list(omp_thread_counts)))
    for n in omp_set:
        configs.append(RunConfig(
            label=f"master_mpi_omp{n}",
            exe="fds_master",
            omp_threads=n,
            mpi_procs=-1,  # sentinel: use case.nmeshes
            description=f"fds_master MPI=nmeshes OMP={n}",
        ))

    # 4. fds_hh single process, pressure subgraph auto
    configs.append(RunConfig(
        label="hh_auto",
        exe="fds_hh",
        omp_threads=1,
        mpi_procs=1,
        pressure_subgraph="auto",
        description="fds_hh (pressure=auto)",
    ))

    # 5. fds_hh single process, pressure subgraph off
    configs.append(RunConfig(
        label="hh_pres_off",
        exe="fds_hh",
        omp_threads=1,
        mpi_procs=1,
        pressure_subgraph="off",
        description="fds_hh (pressure=off)",
    ))

    # 6. fds_hh single process, pressure subgraph on (forced)
    configs.append(RunConfig(
        label="hh_pres_on",
        exe="fds_hh",
        omp_threads=1,
        mpi_procs=1,
        pressure_subgraph="on",
        description="fds_hh (pressure=on)",
    ))

    # 7. fds_hh with MPI (nproc=nmeshes, same as fds_master)
    configs.append(RunConfig(
        label="hh_mpi",
        exe="fds_hh",
        omp_threads=1,
        mpi_procs=-1,  # sentinel: use case.nmeshes
        pressure_subgraph="auto",
        description="fds_hh MPI=nmeshes (pressure=auto)",
    ))

    # 8. fds_hh with MPI (custom count, if requested)
    if mpi_procs >= 2:
        configs.append(RunConfig(
            label=f"hh_mpi{mpi_procs}",
            exe="fds_hh",
            omp_threads=1,
            mpi_procs=mpi_procs,
            pressure_subgraph="auto",
            description=f"fds_hh MPI={mpi_procs} (pressure=auto)",
        ))

    # 9. fds_hh with mesh re-decomposition (target 8x8x8 to ensure splitting)
    if enable_mesh_dim:
        configs.append(RunConfig(
            label="hh_redec_8",
            exe="fds_hh",
            omp_threads=1,
            mpi_procs=1,
            pressure_subgraph="auto",
            mesh_dim=(8, 8, 8),
            description="fds_hh --mesh-dim 8 8 8 (pressure=auto)",
        ))

    return configs


# ── FDS Execution ────────────────────────────────────────────────────────────

def resolve_exe(exe_name: str) -> Path:
    """Resolve executable name to path."""
    if exe_name == "fds_master":
        return FDS_MASTER
    elif exe_name == "fds_hh":
        return FDS_HH
    else:
        return Path(exe_name)


def run_fds(exe: Path, input_file: Path, work_dir: Path, chid: str,
            config: RunConfig, timeout: int = 300) -> Tuple[bool, float]:
    """Run FDS with given configuration. Returns (success, elapsed_seconds)."""
    work_input = work_dir / input_file.name
    shutil.copy2(input_file, work_input)

    # Build command
    cmd = ['mpiexec', '--oversubscribe', '-n', str(config.mpi_procs), str(exe)]
    if config.pressure_subgraph != "auto" and "fds_hh" in exe.name:
        cmd.extend(['--pressure-subgraph', config.pressure_subgraph])
    if config.mesh_dim and "fds_hh" in exe.name:
        cmd.extend(['--mesh-dim',
                     str(config.mesh_dim[0]),
                     str(config.mesh_dim[1]),
                     str(config.mesh_dim[2])])
    cmd.append(input_file.name)

    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(config.omp_threads)

    is_fds_hh = "fds_hh" in exe.name

    try:
        start = time.time()

        if is_fds_hh:
            success = _run_fds_hh(cmd, work_dir, timeout, env)
        else:
            result = subprocess.run(
                cmd, cwd=work_dir, capture_output=True, text=True,
                timeout=timeout, env=env)
            success = result.returncode == 0

        elapsed = time.time() - start

        # Verify via .out file
        out_file = work_dir / f"{chid}.out"
        if out_file.exists():
            content = out_file.read_text()
            if "STOP: FDS completed successfully" in content or "Total Time:" in content:
                # Extract CPU time from .out file for more accurate timing
                cpu_match = re.search(r'Total Elapsed Wall Clock Time.*?:\s*([\d.]+)', content)
                if cpu_match:
                    wall_time = float(cpu_match.group(1))
                    return True, wall_time
                return True, elapsed

        if is_fds_hh and success:
            return True, elapsed

        return False, elapsed

    except subprocess.TimeoutExpired:
        return False, float(timeout)
    except Exception as e:
        print(f"    Error: {e}")
        return False, 0.0


def _run_fds_hh(cmd, work_dir, timeout, env=None) -> bool:
    """Run fds_hh with early termination on completion signal."""
    process = subprocess.Popen(
        cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, bufsize=1, env=env)

    start = time.time()
    try:
        while True:
            if process.poll() is not None:
                return process.returncode == 0

            ready, _, _ = select.select([process.stdout], [], [], 0.5)
            if ready:
                line = process.stdout.readline()
                if not line:
                    continue
                if "Graph dot file written" in line:
                    time.sleep(0.5)
                    process.terminate()
                    try:
                        process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    return True

            if time.time() - start > timeout:
                process.kill()
                process.wait()
                return False
    except Exception:
        process.kill()
        process.wait()
        return False


def clean_work_dir(work_dir: Path, chid: str):
    """Remove binary output files, keeping CSVs and .out."""
    for pattern in ['*.sf', '*.sf.bnd', '*.smv', '*.fds', '*_git.txt',
                    '*.s3d', '*.s3d.bnd', '*.xyz', '*.be', '*.ge', '*.iso',
                    '*.prt5', '*.restart', '*.q', '*.szz',
                    '*.fed', '*.end']:
        for f in work_dir.glob(pattern):
            f.unlink()


# ── Comparison ───────────────────────────────────────────────────────────────

def compare_with_gold(chid: str, run_dir: Path, gold_dir: Path,
                      tolerance: float = 1e-6) -> Tuple[bool, float]:
    """Compare CSV outputs against gold. Returns (passed, max_diff)."""
    skip_suffixes = {'_cpu.csv', '_steps.csv', '_devc_ctrl_log.csv', '_pressit.csv'}
    gold_csvs = sorted(gold_dir.glob(f"{chid}_*.csv"))

    max_diff = 0.0
    all_passed = True

    for gold_file in gold_csvs:
        suffix = gold_file.name[len(chid):]
        if suffix in skip_suffixes:
            continue

        run_file = run_dir / gold_file.name
        if not run_file.exists():
            all_passed = False
            continue

        try:
            cmd = [sys.executable, str(COMPARE_SCRIPT),
                   str(run_file), str(gold_file),
                   '--tolerance', str(tolerance), '--verbose']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode != 0:
                all_passed = False
            for line in result.stdout.split('\n'):
                if 'Max difference' in line:
                    try:
                        d = float(line.split(':')[-1].strip())
                        max_diff = max(max_diff, d)
                    except ValueError:
                        pass
        except Exception:
            all_passed = False

    return all_passed, max_diff


# ── Results I/O ──────────────────────────────────────────────────────────────

def load_results() -> Dict:
    if BENCHMARK_RESULTS.exists():
        with open(BENCHMARK_RESULTS) as f:
            return json.load(f)
    return {}


def save_results(results: Dict):
    BENCHMARK_RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with open(BENCHMARK_RESULTS, 'w') as f:
        json.dump(results, f, indent=2)


# ── Commands ─────────────────────────────────────────────────────────────────

def cmd_generate_gold(args, cases: List[BenchmarkCase]):
    """Generate gold files for benchmark cases using fds_master."""
    exe = resolve_exe("fds_master")
    if not exe.exists():
        print(f"Error: fds_master not found at {exe}")
        return 1

    gold_dir = BENCHMARK_DIR / "gold"
    gold_dir.mkdir(parents=True, exist_ok=True)

    print(f"Generating gold with: {exe}")
    print(f"Output: {gold_dir}")
    print(f"Cases: {len(cases)}, Timeout: {args.timeout}s\n")

    results = load_results()

    for i, case in enumerate(cases, 1):
        case_gold = gold_dir / case.name
        has_gold = (case_gold.exists() and
                    any(case_gold.glob(f"{case.chid}_*.csv")))
        if has_gold and not args.force:
            print(f"[{i}/{len(cases)}] {case.name:50s} SKIP (exists)")
            continue

        case_gold.mkdir(parents=True, exist_ok=True)
        input_file = VERIFICATION_DIR / case.category / case.filename

        # Gold uses MPI with 1 process per mesh (standard FDS parallelism)
        gold_config = RunConfig(
            label="gold", exe="fds_master",
            omp_threads=1, mpi_procs=case.nmeshes)

        print(f"[{i}/{len(cases)}] {case.name:50s} (MPI={case.nmeshes}) ", end='', flush=True)
        success, elapsed = run_fds(exe, input_file, case_gold, case.chid,
                                   gold_config, timeout=args.timeout)
        if success:
            clean_work_dir(case_gold, case.chid)
            print(f"OK ({elapsed:.1f}s)")
        else:
            print(f"FAIL ({elapsed:.1f}s)")

        results.setdefault('gold', {})[case.name] = {
            'success': success,
            'elapsed': round(elapsed, 2),
        }
        save_results(results)

    return 0


def cmd_run(args, cases: List[BenchmarkCase]):
    """Run the full benchmark suite."""
    configs = build_run_configs(
        omp_thread_counts=args.omp_threads,
        enable_mesh_dim=not args.no_mesh_dim,
        mpi_procs=args.mpi_procs,
    )

    # Validate executables
    exes_needed = set(c.exe for c in configs)
    for exe_name in exes_needed:
        exe = resolve_exe(exe_name)
        if not exe.exists():
            print(f"Error: {exe_name} not found at {exe}")
            return 1

    # Check gold files exist
    gold_dir = BENCHMARK_DIR / "gold"
    cases_with_gold = []
    for case in cases:
        case_gold = gold_dir / case.name
        if not case_gold.exists() or not any(case_gold.glob(f"{case.chid}_*.csv")):
            print(f"Warning: No gold for {case.name}. Run 'generate-gold' first.")
        else:
            cases_with_gold.append(case)

    if not cases_with_gold:
        print("No cases with gold files. Run 'generate-gold' first.")
        return 1

    runs_dir = BENCHMARK_DIR / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    results = load_results()
    n_repeats = args.repeats

    total_runs = len(cases_with_gold) * len(configs) * n_repeats
    run_num = 0

    print(f"\nBenchmark Configuration:")
    print(f"  Cases: {len(cases_with_gold)}")
    print(f"  Configs: {len(configs)}")
    print(f"  Repeats: {n_repeats}")
    print(f"  Total runs: {total_runs}")
    print(f"  Configs:")
    for c in configs:
        print(f"    {c.label:25s} {c.description}")
    print()

    for case in cases_with_gold:
        input_file = VERIFICATION_DIR / case.category / case.filename
        case_gold = gold_dir / case.name

        print(f"\n{'='*70}")
        print(f"Case: {case.name} ({case.nmeshes} meshes, {case.total_cells:,} cells, T_END={case.t_end})")
        print(f"{'='*70}")

        for config in configs:
            # Skip mesh_dim for cases with irregular meshes
            if config.mesh_dim and case.name in ('hallways',):
                print(f"  {config.label:25s} SKIP (irregular meshes)")
                continue

            # Resolve sentinel values (-1) to case-specific values
            run_config = config
            if config.mpi_procs == -1 or config.omp_threads == -1:
                import multiprocessing
                hw_threads = multiprocessing.cpu_count()
                run_config = copy.copy(config)
                if config.mpi_procs == -1:
                    run_config.mpi_procs = case.nmeshes
                if config.omp_threads == -1:
                    run_config.omp_threads = max(1, hw_threads // run_config.mpi_procs)
                run_config.description = (config.description
                    .replace('nmeshes', str(run_config.mpi_procs))
                    .replace('hw/nmeshes', str(run_config.omp_threads)))

            # Skip MPI if case has fewer meshes than MPI procs
            if run_config.mpi_procs > case.nmeshes:
                print(f"  {config.label:25s} SKIP (nmeshes={case.nmeshes} < mpi_procs={run_config.mpi_procs})")
                continue

            timings = []
            for rep in range(n_repeats):
                run_num += 1
                work_dir = runs_dir / case.name / f"{config.label}_r{rep}"
                if work_dir.exists():
                    shutil.rmtree(work_dir)
                work_dir.mkdir(parents=True)

                label_str = f"{config.label} [{rep+1}/{n_repeats}]" if n_repeats > 1 else config.label
                print(f"  [{run_num}/{total_runs}] {label_str:30s} ", end='', flush=True)

                success, elapsed = run_fds(
                    resolve_exe(run_config.exe), input_file, work_dir, case.chid,
                    run_config, timeout=case.timeout)

                if success:
                    timings.append(elapsed)
                    # Verify against gold
                    passed, max_diff = compare_with_gold(
                        case.chid, work_dir, case_gold, case.tolerance)
                    status = "OK" if passed else "DIFF"
                    diff_str = f" (max_diff={max_diff:.2e})" if max_diff > 0 else ""
                    print(f"{elapsed:7.1f}s  {status}{diff_str}")
                else:
                    print(f"{elapsed:7.1f}s  FAIL")

                clean_work_dir(work_dir, case.chid)

            # Store results
            if timings:
                key = f"{case.name}/{config.label}"
                results.setdefault('benchmark', {})[key] = {
                    'case': case.name,
                    'config': config.label,
                    'description': config.description,
                    'nmeshes': case.nmeshes,
                    'total_cells': case.total_cells,
                    'timings': timings,
                    'min': round(min(timings), 2),
                    'max': round(max(timings), 2),
                    'mean': round(sum(timings) / len(timings), 2),
                    'repeats': len(timings),
                }
                save_results(results)

    print(f"\n{'='*70}")
    print(f"Benchmark complete. Results in {BENCHMARK_RESULTS}")
    print(f"Run 'python3 {__file__} report' to see summary.")
    return 0


def cmd_verify(args, cases: List[BenchmarkCase]):
    """Verify all fds_hh configurations produce correct results."""
    gold_dir = BENCHMARK_DIR / "gold"
    runs_dir = BENCHMARK_DIR / "runs"

    configs = [
        RunConfig(label="hh_auto", exe="fds_hh", pressure_subgraph="auto",
                  description="fds_hh (pressure=auto)"),
        RunConfig(label="hh_pres_off", exe="fds_hh", pressure_subgraph="off",
                  description="fds_hh (pressure=off)"),
        RunConfig(label="hh_pres_on", exe="fds_hh", pressure_subgraph="on",
                  description="fds_hh (pressure=on)"),
    ]

    exe = resolve_exe("fds_hh")
    if not exe.exists():
        print(f"Error: fds_hh not found at {exe}")
        return 1

    passed = 0
    failed = 0

    for case in cases:
        case_gold = gold_dir / case.name
        if not case_gold.exists():
            print(f"  {case.name:50s} NO GOLD")
            continue

        input_file = VERIFICATION_DIR / case.category / case.filename

        for config in configs:
            work_dir = runs_dir / case.name / f"verify_{config.label}"
            if work_dir.exists():
                shutil.rmtree(work_dir)
            work_dir.mkdir(parents=True)

            print(f"  {case.name}/{config.label:20s} ", end='', flush=True)

            success, elapsed = run_fds(exe, input_file, work_dir, case.chid,
                                       config, timeout=case.timeout)

            if not success:
                print(f"RUN_FAIL ({elapsed:.1f}s)")
                failed += 1
                continue

            ok, max_diff = compare_with_gold(
                case.chid, work_dir, case_gold, case.tolerance)

            if ok:
                print(f"PASS ({elapsed:.1f}s)")
                passed += 1
            else:
                print(f"FAIL (max_diff={max_diff:.2e}, {elapsed:.1f}s)")
                failed += 1

            clean_work_dir(work_dir, case.chid)

    print(f"\nVerification: {passed} passed, {failed} failed")
    return 0 if failed == 0 else 1


def cmd_report(args, cases: List[BenchmarkCase]):
    """Print benchmark results report."""
    results = load_results()
    benchmark = results.get('benchmark', {})

    if not benchmark:
        print("No benchmark results. Run 'run' first.")
        return 1

    # Group by case
    by_case: Dict[str, Dict[str, Dict]] = {}
    for key, data in benchmark.items():
        case_name = data['case']
        config_label = data['config']
        by_case.setdefault(case_name, {})[config_label] = data

    print(f"\n{'='*90}")
    print(f"FDS BENCHMARK RESULTS")
    print(f"{'='*90}")

    for case in cases:
        if case.name not in by_case:
            continue

        case_data = by_case[case.name]
        print(f"\n{'─'*90}")
        print(f"Case: {case.name}")
        print(f"  Meshes: {case.nmeshes}, Cells: {case.total_cells:,}, T_END: {case.t_end}")
        print(f"  Pressure subgraph eligible: {case.pressure_eligible}")
        print(f"{'─'*90}")

        # Find best fds_master time for speedup calculation
        master_times = {}
        for label, data in case_data.items():
            if label.startswith('master_'):
                master_times[label] = data['min']

        best_master_label = min(master_times, key=master_times.get) if master_times else None
        best_master_time = master_times.get(best_master_label, 0) if best_master_label else 0

        print(f"\n  {'Config':<30s} {'Min':>8s} {'Mean':>8s} {'Max':>8s} {'Reps':>5s} {'Speedup':>8s}")
        print(f"  {'─'*30} {'─'*8} {'─'*8} {'─'*8} {'─'*5} {'─'*8}")

        # Print fds_master configs first, then fds_hh
        for label in sorted(case_data.keys(), key=lambda l: (0 if l.startswith('master') else 1, l)):
            data = case_data[label]
            speedup = ""
            if best_master_time > 0 and data['min'] > 0:
                ratio = best_master_time / data['min']
                speedup = f"{ratio:.2f}x"
                if label == best_master_label:
                    speedup = "(ref)"

            desc = data.get('description', label)
            print(f"  {desc:<30s} {data['min']:7.1f}s {data['mean']:7.1f}s {data['max']:7.1f}s"
                  f" {data['repeats']:>5d} {speedup:>8s}")

        if best_master_label:
            print(f"\n  Reference: {best_master_label} ({best_master_time:.1f}s)")

    # Overall summary
    print(f"\n{'='*90}")
    print(f"SUMMARY")
    print(f"{'='*90}")

    # For each config, compute geometric mean of speedups across cases
    config_speedups: Dict[str, List[float]] = {}
    for case in cases:
        if case.name not in by_case:
            continue
        case_data = by_case[case.name]

        master_times = {l: d['min'] for l, d in case_data.items() if l.startswith('master_')}
        if not master_times:
            continue
        best_master = min(master_times.values())

        for label, data in case_data.items():
            if data['min'] > 0 and best_master > 0:
                ratio = best_master / data['min']
                config_speedups.setdefault(label, []).append(ratio)

    print(f"\n  {'Config':<30s} {'Geomean Speedup':>15s} {'Cases':>6s}")
    print(f"  {'─'*30} {'─'*15} {'─'*6}")
    for label in sorted(config_speedups.keys(), key=lambda l: (0 if l.startswith('master') else 1, l)):
        speedups = config_speedups[label]
        geomean = 1.0
        for s in speedups:
            geomean *= s
        geomean = geomean ** (1.0 / len(speedups))
        print(f"  {label:<30s} {geomean:14.2f}x {len(speedups):>6d}")

    return 0


# ── Case Filtering ───────────────────────────────────────────────────────────

def get_cases(args) -> List[BenchmarkCase]:
    """Get benchmark cases, optionally filtered by name."""
    cases = BENCHMARK_CASES
    if hasattr(args, 'case') and args.case:
        cases = [c for c in cases if c.name in args.case]
    return cases


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='FDS Benchmark Suite',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # generate-gold
    p_gold = subparsers.add_parser('generate-gold', help='Generate gold files')
    p_gold.add_argument('--timeout', type=int, default=1800)
    p_gold.add_argument('--force', action='store_true')
    p_gold.add_argument('--case', nargs='+', help='Run specific cases only')

    # run
    p_run = subparsers.add_parser('run', help='Run benchmark')
    p_run.add_argument('--omp-threads', type=int, nargs='+', default=[2, 4],
                       help='OMP thread counts for fds_master (default: 2 4). '
                            'OMP=1 is always included. Avoid oversubscription: '
                            'nmeshes * omp_threads <= hw_threads')
    p_run.add_argument('--repeats', type=int, default=3,
                       help='Number of repeat runs per config (default: 3)')
    p_run.add_argument('--mpi-procs', type=int, default=0,
                       help='Number of MPI processes for fds_hh MPI test (0=skip)')
    p_run.add_argument('--no-mesh-dim', action='store_true',
                       help='Skip mesh re-decomposition tests')
    p_run.add_argument('--case', nargs='+', help='Run specific cases only')

    # verify
    p_verify = subparsers.add_parser('verify', help='Verify correctness')
    p_verify.add_argument('--case', nargs='+', help='Verify specific cases only')

    # report
    p_report = subparsers.add_parser('report', help='Show results')
    p_report.add_argument('--case', nargs='+', help='Report specific cases only')

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        return 1

    cases = get_cases(args)
    if not cases:
        print("No matching cases.")
        return 1

    if args.command == 'generate-gold':
        return cmd_generate_gold(args, cases)
    elif args.command == 'run':
        return cmd_run(args, cases)
    elif args.command == 'verify':
        return cmd_verify(args, cases)
    elif args.command == 'report':
        return cmd_report(args, cases)

    return 0


if __name__ == '__main__':
    sys.exit(main())
