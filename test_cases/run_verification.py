#!/usr/bin/env python3
"""
FDS Verification Suite Runner

Runs verification cases from Verification/FDS_Cases.sh,
comparing fds_hh (Hedgehog) output against fds_master (ground truth).

Usage:
    # Discover all cases
    python run_verification.py discover [--category Aerosols]

    # Generate gold files with fds_master (ground truth)
    python run_verification.py generate-gold --exe /path/to/fds_master [--category Aerosols] [--timeout 300] [--jobs 4]

    # Run tests with fds_hh and compare against gold
    python run_verification.py test [--category Aerosols] [--timeout 300] [--jobs 4]

    # Show results report (correctness + performance)
    python run_verification.py report
"""

import os
import sys
import subprocess
import argparse
import json
import re
import select
import shutil
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

# ── Paths ──────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent
VERIFICATION_DIR = REPO_ROOT / "Verification"
FDS_CASES_SCRIPT = VERIFICATION_DIR / "FDS_Cases.sh"
TEST_DIR = Path(__file__).resolve().parent
GOLD_DIR = TEST_DIR / "verification_gold"
RUN_DIR = TEST_DIR / "verification_run"
RESULTS_FILE = TEST_DIR / "verification_results.json"
BUILD_DIR = REPO_ROOT / "build_hh"
FDS_HH = BUILD_DIR / "Source" / "hedgehog" / "fds_hh"
FDS_MASTER_PATH = shutil.which('fds_master')
FDS_MASTER = Path(FDS_MASTER_PATH) if FDS_MASTER_PATH else None
COMPARE_SCRIPT = TEST_DIR / "compare_csv.py"


# ── FDS Case Parsing ──────────────────────────────────────────────────────────

def parse_fds_cases() -> List[Dict]:
    """Parse FDS_Cases.sh to extract all active test case definitions."""
    cases = []
    with open(FDS_CASES_SCRIPT) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('$QFDS') or line.startswith('#'):
                continue

            # Parse: $QFDS [-p N] -d DIR FILE [other flags...]
            m = re.match(r'\$QFDS\s+(?:-p\s+(\d+)\s+)?-d\s+(\S+)\s+(\S+)', line)
            if not m:
                continue

            nproc = int(m.group(1)) if m.group(1) else 1
            category = m.group(2)
            filename = m.group(3)

            fds_file = VERIFICATION_DIR / category / filename
            if not fds_file.exists():
                continue

            chid = extract_chid(fds_file)
            nmeshes = count_meshes(fds_file)
            t_end = extract_t_end(fds_file)

            cases.append({
                'category': category,
                'filename': filename,
                'chid': chid,
                'nproc': nproc,
                'nmeshes': nmeshes,
                't_end': t_end,
                'id': f"{category}/{chid}",
            })

    return cases


def extract_chid(fds_file: Path) -> str:
    """Extract CHID from &HEAD line."""
    with open(fds_file) as f:
        content = f.read()
    m = re.search(r"CHID\s*=\s*'([^']+)'", content)
    return m.group(1) if m else fds_file.stem


def count_meshes(fds_file: Path) -> int:
    """Count meshes including MULT expansion."""
    with open(fds_file) as f:
        content = f.read()

    # Remove comment lines for parsing
    lines = [l for l in content.split('\n') if not l.strip().startswith('!')]
    content = '\n'.join(lines)

    mesh_lines = re.findall(r'&MESH\b(.*?)/', content, re.DOTALL)
    total = 0

    for mesh_body in mesh_lines:
        mult_match = re.search(r"MULT_ID\s*=\s*'([^']+)'", mesh_body)
        if mult_match:
            mult_id = mult_match.group(1)
            mult_pattern = rf"&MULT\s[^/]*ID\s*=\s*'{re.escape(mult_id)}'(.*?)/"
            mult_m = re.search(mult_pattern, content, re.DOTALL)
            if mult_m:
                mult_def = mult_m.group(1)
                ni = _get_upper(mult_def, 'I_UPPER') + 1
                nj = _get_upper(mult_def, 'J_UPPER') + 1
                nk = _get_upper(mult_def, 'K_UPPER') + 1
                total += ni * nj * nk
            else:
                total += 1
        else:
            total += 1

    return max(total, 1)


def _get_upper(text: str, key: str) -> int:
    m = re.search(rf'{key}\s*=\s*(\d+)', text)
    return int(m.group(1)) if m else 0


def extract_t_end(fds_file: Path) -> float:
    """Extract T_END from &TIME line."""
    with open(fds_file) as f:
        content = f.read()
    m = re.search(r'T_END\s*=\s*([\d.eE+-]+)', content)
    return float(m.group(1)) if m else 0.0


# ── FDS Execution ─────────────────────────────────────────────────────────────

def run_fds(exe: Path, input_file: Path, work_dir: Path, chid: str,
            timeout: int = 300, is_fds_hh: bool = False,
            nproc: int = 1, omp_threads: int = 1) -> Tuple[bool, float]:
    """
    Run an FDS executable on the given input file.

    Args:
        nproc: Number of MPI processes (from -p N in FDS_Cases.sh).
               For fds_hh, always pass 1 (Hedgehog handles all meshes internally).
        omp_threads: Number of OpenMP threads (for fds_master performance).
               For fds_hh, always pass 1 (uses Hedgehog threading instead).

    Returns (success, elapsed_seconds).
    """
    # Copy input file to work directory
    work_input = work_dir / input_file.name
    shutil.copy2(input_file, work_input)

    cmd = ['mpiexec', '--oversubscribe', '-n', str(nproc), str(exe), input_file.name]

    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(omp_threads)

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

        # Verify completion via .out file
        out_file = work_dir / f"{chid}.out"
        if out_file.exists():
            content = out_file.read_text()
            if "STOP: FDS completed successfully" in content or "Total Time:" in content:
                return True, elapsed

        # fds_hh: completion detected via stdout signal
        if is_fds_hh and success:
            return True, elapsed

        # Check if we hit the timeout
        if elapsed >= timeout - 1:
            return False, elapsed  # Timeout

        return False, elapsed

    except subprocess.TimeoutExpired:
        return False, float(timeout)
    except Exception as e:
        return False, 0.0


def _run_fds_hh(cmd, work_dir, timeout, env=None) -> bool:
    """Run fds_hh with early termination on 'Graph dot file written'."""
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
    """Remove non-essential output files, keeping only CSVs."""
    for pattern in ['*.sf', '*.sf.bnd', '*.smv', '*.out', '*.fds', '*_git.txt',
                    '*.s3d', '*.s3d.bnd', '*.xyz', '*.be', '*.ge', '*.iso',
                    '*.prt5', '*.restart', '*.q', '*.szz',
                    '*.fed', '*.end']:
        for f in work_dir.glob(pattern):
            f.unlink()


# ── Comparison ────────────────────────────────────────────────────────────────

def compare_outputs(chid: str, run_dir: Path, gold_dir: Path,
                    tolerance: float = 1e-10) -> Dict:
    """Compare all CSV output files between run and gold directories.

    Discovers gold CSV files dynamically (devc, hrr, and any others).
    Skips timing-only files (cpu, steps).
    """
    results = {}

    # Discover all gold CSV files for this CHID
    skip_suffixes = {'_cpu.csv', '_steps.csv', '_devc_ctrl_log.csv'}
    gold_csvs = sorted(gold_dir.glob(f"{chid}_*.csv"))

    for gold_file in gold_csvs:
        suffix = gold_file.name[len(chid):]  # e.g., "_devc.csv"
        if suffix in skip_suffixes:
            continue

        run_file = run_dir / gold_file.name

        if not run_file.exists():
            results[suffix] = {'status': 'missing', 'passed': False}
            continue

        try:
            cmd = [sys.executable, str(COMPARE_SCRIPT),
                   str(run_file), str(gold_file),
                   '--tolerance', str(tolerance), '--verbose']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            passed = result.returncode == 0
            # Extract max diff from output
            max_diff = 0.0
            for line in result.stdout.split('\n'):
                if 'Max difference' in line:
                    try:
                        max_diff = float(line.split(':')[-1].strip())
                    except ValueError:
                        pass
            results[suffix] = {
                'status': 'compared',
                'passed': passed,
                'max_diff': max_diff,
                'output': result.stdout[-500:] if not passed else '',
            }
        except Exception as e:
            results[suffix] = {'status': 'error', 'passed': False, 'output': str(e)}

    return results


# ── Results I/O ───────────────────────────────────────────────────────────────

def load_results() -> Dict:
    """Load results from JSON file."""
    if RESULTS_FILE.exists():
        with open(RESULTS_FILE) as f:
            return json.load(f)
    return {}


def save_results(results: Dict):
    """Save results to JSON file."""
    with open(RESULTS_FILE, 'w') as f:
        json.dump(results, f, indent=2)


# ── Commands ──────────────────────────────────────────────────────────────────

def cmd_discover(args, cases: List[Dict]):
    """Print a summary of all discovered cases."""
    categories = {}
    for c in cases:
        cat = c['category']
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(c)

    total = len(cases)
    multi = sum(1 for c in cases if c['nmeshes'] > 1)

    print(f"Discovered {total} cases ({multi} multi-mesh) in {len(categories)} categories\n")

    for cat in sorted(categories.keys()):
        cat_cases = categories[cat]
        print(f"  {cat:40s} {len(cat_cases):3d} cases")

        if args.verbose:
            for c in cat_cases:
                mesh_str = f"{c['nmeshes']}m" if c['nmeshes'] > 1 else "  "
                print(f"    {c['chid']:50s} {mesh_str:>4s}  T_END={c['t_end']}")

    print(f"\n  {'TOTAL':40s} {total:3d} cases")
    print(f"  Multi-mesh (>1):  {multi}")
    print(f"  Single-mesh:      {total - multi}")


def _run_single_case(case: Dict, exe: Path, output_dir: Path,
                     timeout: int, is_fds_hh: bool,
                     omp_threads: int = 1) -> Dict:
    """Run a single case (used by both generate-gold and test). Returns result dict."""
    case_id = case['id']
    chid = case['chid']
    category = case['category']
    input_file = VERIFICATION_DIR / category / case['filename']

    work_dir = output_dir / category / chid
    work_dir.mkdir(parents=True, exist_ok=True)

    # fds_hh always runs with 1 MPI process (Hedgehog handles all meshes internally)
    # fds_master uses the nproc from FDS_Cases.sh for multi-mesh parallelism
    nproc = 1 if is_fds_hh else case.get('nproc', 1)
    # fds_hh uses Hedgehog threading, not OpenMP
    omp = 1 if is_fds_hh else omp_threads

    success, elapsed = run_fds(
        exe, input_file, work_dir, chid,
        timeout=timeout, is_fds_hh=is_fds_hh,
        nproc=nproc, omp_threads=omp)

    # Clean up binary files
    if success:
        clean_work_dir(work_dir, chid)

    return {
        'case_id': case_id,
        'chid': chid,
        'category': category,
        'nmeshes': case['nmeshes'],
        'nproc': nproc,
        'omp_threads': omp,
        'success': success,
        'elapsed': round(elapsed, 2),
    }


def cmd_generate_gold(args, cases: List[Dict]):
    """Generate gold files using fds_master."""
    exe = Path(args.exe) if args.exe else FDS_MASTER
    if exe is None:
        print("Error: fds_master not found in PATH. Use --exe /path/to/fds_master")
        return 1
    if not exe.exists():
        print(f"Error: {exe} does not exist")
        return 1

    print(f"Generating gold files with: {exe}")
    print(f"Output directory: {GOLD_DIR}")
    print(f"Cases: {len(cases)}, Timeout: {args.timeout}s")
    print(f"OMP_NUM_THREADS: {args.omp_threads}, MPI: per-case nproc from FDS_Cases.sh\n")

    GOLD_DIR.mkdir(parents=True, exist_ok=True)
    results = load_results()

    # Skip cases that already have successful gold (unless --force)
    to_run = []
    for case in cases:
        gold_case_dir = GOLD_DIR / case['category'] / case['chid']
        has_gold = (gold_case_dir.exists() and
                    any(gold_case_dir.glob(f"{case['chid']}_*.csv")))
        # Also check if the gold run succeeded (don't skip failed/timed-out gold)
        gold_result = results.get('gold', {}).get(case['id'], {})
        gold_succeeded = gold_result.get('success', False) if gold_result else False
        if has_gold and gold_succeeded and not args.force:
            continue
        to_run.append(case)

    skipped = len(cases) - len(to_run)
    if skipped:
        print(f"Skipping {skipped} cases with existing gold (use --force to regenerate)\n")

    if not to_run:
        print("Nothing to do.")
        return 0

    passed = 0
    failed = 0
    failed_cases = []

    if args.jobs <= 1:
        for i, case in enumerate(to_run, 1):
            nproc_str = f" -p {case['nproc']}" if case['nproc'] > 1 else ""
            print(f"[{i}/{len(to_run)}] {case['id']:50s}{nproc_str:>5s} ", end='', flush=True)
            result = _run_single_case(case, exe, GOLD_DIR, args.timeout, False,
                                      omp_threads=args.omp_threads)
            if result['success']:
                passed += 1
                print(f"OK  ({result['elapsed']:.1f}s)")
            elif result['elapsed'] >= args.timeout - 1:
                failed += 1
                failed_cases.append(f"{case['id']} (TIMEOUT)")
                print(f"TIMEOUT ({result['elapsed']:.0f}s)")
            else:
                failed += 1
                failed_cases.append(case['id'])
                print(f"FAIL ({result['elapsed']:.1f}s)")

            results.setdefault('gold', {})[case['id']] = result
            save_results(results)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {}
            for case in to_run:
                fut = executor.submit(
                    _run_single_case, case, exe, GOLD_DIR, args.timeout, False,
                    omp_threads=args.omp_threads)
                futures[fut] = case

            for i, fut in enumerate(as_completed(futures), 1):
                case = futures[fut]
                try:
                    result = fut.result()
                except Exception as e:
                    result = {'case_id': case['id'], 'success': False, 'elapsed': 0}
                if result['success']:
                    passed += 1
                    status = "OK"
                elif result.get('elapsed', 0) >= args.timeout - 1:
                    failed += 1
                    failed_cases.append(f"{case['id']} (TIMEOUT)")
                    status = "TIMEOUT"
                else:
                    failed += 1
                    failed_cases.append(case['id'])
                    status = "FAIL"
                print(f"[{i}/{len(to_run)}] {case['id']:50s} {status:4s} ({result.get('elapsed', 0):.1f}s)")

                results.setdefault('gold', {})[case['id']] = result
                save_results(results)

    print(f"\nGold generation: {passed} passed, {failed} failed out of {len(to_run)}")
    if failed_cases:
        print("Failed cases:")
        for c in failed_cases:
            print(f"  - {c}")
    return 0 if failed == 0 else 1


def cmd_test(args, cases: List[Dict]):
    """Run tests with fds_hh and compare against gold."""
    if not FDS_HH.exists():
        print(f"Error: fds_hh not found at {FDS_HH}")
        print(f"Build with: cd build_hh && cmake --build . --target fds_hh -j$(nproc)")
        return 1

    # fds_hh uses Hedgehog internal threading; always run 1 at a time
    if args.jobs > 1:
        print("Note: fds_hh uses internal threading; forcing jobs=1")
        args.jobs = 1

    print(f"Testing with: {FDS_HH}")
    print(f"Gold directory: {GOLD_DIR}")
    print(f"Cases: {len(cases)}, Timeout: {args.timeout}s")
    print(f"Tolerance: {args.tolerance:.0e}")
    print(f"fds_hh: MPI=1, OMP=1 (uses Hedgehog threading)\n")

    RUN_DIR.mkdir(parents=True, exist_ok=True)
    results = load_results()

    # Filter to cases that have complete gold files
    to_run = []
    no_gold = 0
    incomplete_gold = 0
    for case in cases:
        gold_case_dir = GOLD_DIR / case['category'] / case['chid']
        has_gold = (gold_case_dir.exists() and
                    any(gold_case_dir.glob(f"{case['chid']}_*.csv")))
        if not has_gold:
            no_gold += 1
            continue
        # Check if gold run completed successfully
        gold_result = results.get('gold', {}).get(case['id'], {})
        if gold_result and not gold_result.get('success', True):
            incomplete_gold += 1
            continue
        to_run.append(case)

    if no_gold:
        print(f"Skipping {no_gold} cases without gold files "
              f"(run 'generate-gold' first)")
    if incomplete_gold:
        print(f"Skipping {incomplete_gold} cases with incomplete gold "
              f"(gold run timed out or failed)")
    if no_gold or incomplete_gold:
        print()

    if not to_run:
        print("No cases to test.")
        return 0

    passed = 0
    failed = 0
    errored = 0
    test_results = []

    for i, case in enumerate(to_run, 1):
        print(f"[{i}/{len(to_run)}] {case['id']:50s} ", end='', flush=True)

        # Run fds_hh
        result = _run_single_case(case, FDS_HH, RUN_DIR, args.timeout, True)

        if not result['success']:
            errored += 1
            result['status'] = 'run_failed'
            print(f"RUN_FAIL ({result['elapsed']:.1f}s)")
        else:
            # Compare with gold
            run_case_dir = RUN_DIR / case['category'] / case['chid']
            gold_case_dir = GOLD_DIR / case['category'] / case['chid']
            comparison = compare_outputs(
                case['chid'], run_case_dir, gold_case_dir, args.tolerance)

            result['comparison'] = comparison

            if not comparison:
                # No comparable files found
                result['status'] = 'no_output'
                errored += 1
                print(f"NO_OUTPUT ({result['elapsed']:.1f}s)")
            elif all(c['passed'] for c in comparison.values()):
                result['status'] = 'passed'
                passed += 1
                # Get gold timing for comparison
                gold_result = results.get('gold', {}).get(case['id'], {})
                gold_time = gold_result.get('elapsed', 0)
                if gold_time > 0:
                    speedup = gold_time / result['elapsed'] if result['elapsed'] > 0 else 0
                    print(f"PASS ({result['elapsed']:.1f}s, "
                          f"fds_master={gold_time:.1f}s, "
                          f"ratio={speedup:.2f}x)")
                else:
                    print(f"PASS ({result['elapsed']:.1f}s)")
            else:
                result['status'] = 'failed'
                failed += 1
                # Show which files failed
                fail_files = [k for k, v in comparison.items() if not v['passed']]
                print(f"FAIL ({result['elapsed']:.1f}s) [{', '.join(fail_files)}]")

        results.setdefault('test', {})[case['id']] = result
        save_results(results)
        test_results.append(result)

    # Summary
    total = len(to_run)
    print(f"\n{'=' * 70}")
    print(f"TEST SUMMARY: {passed} passed, {failed} failed, "
          f"{errored} errors out of {total}")
    print(f"{'=' * 70}")

    if failed > 0:
        print("\nFailed cases:")
        for r in test_results:
            if r.get('status') == 'failed':
                print(f"  - {r['case_id']}")

    if errored > 0:
        print("\nErrored cases:")
        for r in test_results:
            if r.get('status') in ('run_failed', 'no_output'):
                print(f"  - {r['case_id']} ({r['status']})")

    return 0 if (failed == 0 and errored == 0) else 1


def cmd_report(args, cases: List[Dict]):
    """Show results report with correctness and performance comparison."""
    results = load_results()
    gold_results = results.get('gold', {})
    test_results = results.get('test', {})

    if not test_results:
        print("No test results found. Run 'test' first.")
        return 1

    # Collect performance data
    perf_data = []
    status_counts = {'passed': 0, 'failed': 0, 'run_failed': 0, 'no_output': 0}

    for case_id, tr in sorted(test_results.items()):
        status = tr.get('status', 'unknown')
        status_counts[status] = status_counts.get(status, 0) + 1

        gr = gold_results.get(case_id, {})
        gold_time = gr.get('elapsed', 0)
        hh_time = tr.get('elapsed', 0)

        if status == 'passed' and gold_time > 0 and hh_time > 0:
            perf_data.append({
                'case_id': case_id,
                'nmeshes': tr.get('nmeshes', 1),
                'nproc': gr.get('nproc', 1),
                'fds_master_time': gold_time,
                'fds_hh_time': hh_time,
                'ratio': gold_time / hh_time if hh_time > 0 else 0,
            })

    # Print correctness summary
    total = sum(status_counts.values())
    print(f"{'=' * 80}")
    print(f"VERIFICATION REPORT")
    print(f"{'=' * 80}")
    print(f"\nCorrectness: {status_counts.get('passed', 0)} passed, "
          f"{status_counts.get('failed', 0)} failed, "
          f"{status_counts.get('run_failed', 0)} run errors, "
          f"{status_counts.get('no_output', 0)} no output "
          f"(total: {total})")

    # Print failed cases with max_diff info
    failed = [(k, v) for k, v in sorted(test_results.items())
              if v.get('status') == 'failed']
    if failed:
        print(f"\nFailed cases ({len(failed)}):")
        for case_id, tr in failed:
            comp = tr.get('comparison', {})
            fail_files = []
            max_case_diff = 0.0
            for k, v in comp.items():
                if not v.get('passed', True):
                    fail_files.append(k)
                    max_case_diff = max(max_case_diff, v.get('max_diff', 0.0))
            diff_str = f" max_diff={max_case_diff:.2e}" if max_case_diff > 0 else ""
            print(f"  {case_id:50s} [{', '.join(fail_files)}]{diff_str}")

    # Print performance summary
    if perf_data:
        print(f"\n{'=' * 80}")
        print(f"PERFORMANCE COMPARISON (fds_master vs fds_hh)")
        print(f"{'=' * 80}")

        # Separate single-mesh and multi-mesh
        single = [p for p in perf_data if p['nmeshes'] <= 1]
        multi = [p for p in perf_data if p['nmeshes'] > 1]

        for label, data in [("Single-mesh", single), ("Multi-mesh", multi)]:
            if not data:
                continue

            total_fds_master = sum(p['fds_master_time'] for p in data)
            total_hh = sum(p['fds_hh_time'] for p in data)
            avg_ratio = total_fds_master / total_hh if total_hh > 0 else 0

            print(f"\n{label} cases ({len(data)}):")
            print(f"  Total fds_master time:   {total_fds_master:8.1f}s")
            print(f"  Total fds_hh time: {total_hh:8.1f}s")
            print(f"  Aggregate ratio:   {avg_ratio:8.2f}x")

            # Show top 20 cases by fds_hh time
            data_sorted = sorted(data, key=lambda p: p['fds_hh_time'], reverse=True)
            print(f"\n  {'Case':50s} {'fds_master':>8s} {'fds_hh':>8s} {'ratio':>7s} {'meshes':>6s} {'nproc':>5s}")
            print(f"  {'-'*50} {'-'*8} {'-'*8} {'-'*7} {'-'*6} {'-'*5}")
            for p in data_sorted[:20]:
                nproc_str = str(p.get('nproc', '?'))
                print(f"  {p['case_id']:50s} {p['fds_master_time']:7.1f}s {p['fds_hh_time']:7.1f}s "
                      f"{p['ratio']:6.2f}x {p['nmeshes']:5d} {nproc_str:>5s}")

        # Overall
        total_fds_master = sum(p['fds_master_time'] for p in perf_data)
        total_hh = sum(p['fds_hh_time'] for p in perf_data)
        avg_ratio = total_fds_master / total_hh if total_hh > 0 else 0

        print(f"\n{'=' * 80}")
        print(f"Overall ({len(perf_data)} cases): "
              f"fds_master={total_fds_master:.1f}s, fds_hh={total_hh:.1f}s, "
              f"ratio={avg_ratio:.2f}x")
        print(f"{'=' * 80}")

    return 0


# ── Case Filtering ────────────────────────────────────────────────────────────

def filter_cases(cases: List[Dict], args) -> List[Dict]:
    """Filter cases based on command-line arguments."""
    filtered = cases

    if args.category:
        cats = [c.strip() for c in args.category.split(',')]
        filtered = [c for c in filtered if c['category'] in cats]

    if hasattr(args, 'max_meshes') and args.max_meshes:
        filtered = [c for c in filtered if c['nmeshes'] <= args.max_meshes]

    if hasattr(args, 'case') and args.case:
        patterns = args.case
        filtered = [c for c in filtered
                    if any(p in c['id'] or p in c['chid'] for p in patterns)]

    if hasattr(args, 'exclude_category') and args.exclude_category:
        excl = [c.strip() for c in args.exclude_category.split(',')]
        filtered = [c for c in filtered if c['category'] not in excl]

    if hasattr(args, 'min_tend') and args.min_tend is not None:
        filtered = [c for c in filtered if c['t_end'] > args.min_tend]

    if hasattr(args, 'no_redundant') and args.no_redundant:
        filtered = _remove_redundant(filtered)

    if hasattr(args, 'max_gold_time') and args.max_gold_time is not None:
        results = load_results()
        gold = results.get('gold', {})
        before = len(filtered)
        filtered = [c for c in filtered
                    if gold.get(c['id'], {}).get('success', False)
                    and gold[c['id']].get('elapsed', 0) <= args.max_gold_time]
        removed = before - len(filtered)
        if removed:
            print(f"Filtered {removed} cases exceeding gold time {args.max_gold_time}s")

    return filtered


def _remove_redundant(cases: List[Dict]) -> List[Dict]:
    """Remove redundant grid-refinement cases.

    For families of cases that differ only by a grid size number,
    keep only the smallest single-mesh version and all multi-mesh versions.

    Detects patterns like: ns2d_8/ns2d_16/ns2d_64,
    ns2d_8_nupt1/ns2d_16_nupt1/ns2d_64_nupt1,
    plate_view_factor_2D_30/60/100, vort2d_40/80/160/320.
    """
    # Group by (category, pattern) where pattern replaces number sequences with '#'
    families: Dict[tuple, List[Dict]] = {}
    for c in cases:
        pattern = re.sub(r'\d+', '#', c['chid'])
        key = (c['category'], pattern)
        families.setdefault(key, []).append(c)

    result = []
    for (cat, pattern), members in families.items():
        if len(members) <= 1:
            result.extend(members)
            continue

        single_mesh = [m for m in members if m['nmeshes'] <= 1]
        multi_mesh = [m for m in members if m['nmeshes'] > 1]

        if len(single_mesh) <= 1:
            # At most one single-mesh - no redundancy to remove
            result.extend(single_mesh)
        else:
            # Multiple single-mesh cases with the same pattern.
            # Find the varying number position(s) and keep the smallest.
            # Extract all numbers from each CHID
            nums_per_case = []
            for m in single_mesh:
                nums = [int(x) for x in re.findall(r'\d+', m['chid'])]
                nums_per_case.append(nums)

            # Find which number position varies
            if nums_per_case and all(len(n) == len(nums_per_case[0]) for n in nums_per_case):
                # Find varying positions
                varying_pos = []
                for i in range(len(nums_per_case[0])):
                    vals = set(n[i] for n in nums_per_case)
                    if len(vals) > 1:
                        varying_pos.append(i)

                if len(varying_pos) == 1:
                    # Single varying number - sort by it, keep smallest
                    pos = varying_pos[0]
                    sorted_cases = sorted(single_mesh,
                                          key=lambda m: [int(x) for x in
                                                         re.findall(r'\d+', m['chid'])][pos])
                    result.append(sorted_cases[0])
                else:
                    # Multiple varying positions - keep all
                    result.extend(single_mesh)
            else:
                result.extend(single_mesh)

        # Always keep multi-mesh variants
        result.extend(multi_mesh)

    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='FDS Verification Suite Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(dest='command', help='Command to run')

    # discover
    p_disc = subparsers.add_parser('discover', help='List all verification cases')
    p_disc.add_argument('-v', '--verbose', action='store_true')
    p_disc.add_argument('--category', help='Filter by category (comma-separated)')
    p_disc.add_argument('--max-meshes', type=int, help='Max mesh count')
    p_disc.add_argument('--case', action='append', help='Filter by case name/CHID')
    p_disc.add_argument('--exclude-category', help='Exclude categories (comma-separated)')
    p_disc.add_argument('--min-tend', type=float, help='Skip cases with T_END <= value (e.g., 0 to skip viz-only)')
    p_disc.add_argument('--no-redundant', action='store_true', help='Keep only smallest grid size per family + multi-mesh')
    p_disc.add_argument('--max-gold-time', type=float, help='Skip cases where gold took longer than N seconds')

    # generate-gold
    p_gold = subparsers.add_parser('generate-gold', help='Generate gold files with fds_master')
    p_gold.add_argument('--exe', help='Path to fds_master binary (default: fds_master on PATH)')
    p_gold.add_argument('--category', help='Filter by category (comma-separated)')
    p_gold.add_argument('--max-meshes', type=int, help='Max mesh count')
    p_gold.add_argument('--case', action='append', help='Filter by case name/CHID')
    p_gold.add_argument('--exclude-category', help='Exclude categories (comma-separated)')
    p_gold.add_argument('--min-tend', type=float, help='Skip cases with T_END <= value (e.g., 0 to skip viz-only)')
    p_gold.add_argument('--timeout', type=int, default=300, help='Per-case timeout (default: 300s)')
    p_gold.add_argument('--jobs', '-j', type=int, default=1, help='Parallel jobs (default: 1)')
    p_gold.add_argument('--force', action='store_true', help='Regenerate existing gold files')
    p_gold.add_argument('--omp-threads', type=int, default=4, help='OpenMP threads for fds_master (default: 4)')
    p_gold.add_argument('--no-redundant', action='store_true', help='Keep only smallest grid size per family + multi-mesh')
    p_gold.add_argument('--max-gold-time', type=float, help='Skip cases where gold took longer than N seconds')

    # test
    p_test = subparsers.add_parser('test', help='Run fds_hh tests and compare with gold')
    p_test.add_argument('--category', help='Filter by category (comma-separated)')
    p_test.add_argument('--max-meshes', type=int, help='Max mesh count')
    p_test.add_argument('--case', action='append', help='Filter by case name/CHID')
    p_test.add_argument('--exclude-category', help='Exclude categories (comma-separated)')
    p_test.add_argument('--min-tend', type=float, help='Skip cases with T_END <= value (e.g., 0 to skip viz-only)')
    p_test.add_argument('--timeout', type=int, default=300, help='Per-case timeout (default: 300s)')
    p_test.add_argument('--jobs', '-j', type=int, default=1, help='Parallel jobs (default: 1)')
    p_test.add_argument('--tolerance', type=float, default=1e-10, help='Comparison tolerance (default: 1e-10)')
    p_test.add_argument('--no-redundant', action='store_true', help='Keep only smallest grid size per family + multi-mesh')
    p_test.add_argument('--max-gold-time', type=float, help='Skip cases where gold took longer than N seconds')

    # report
    p_report = subparsers.add_parser('report', help='Show results report')
    p_report.add_argument('--category', help='Filter by category (comma-separated)')
    p_report.add_argument('--max-meshes', type=int, help='Max mesh count')
    p_report.add_argument('--case', action='append', help='Filter by case name/CHID')
    p_report.add_argument('--exclude-category', help='Exclude categories (comma-separated)')
    p_report.add_argument('--min-tend', type=float, help='Skip cases with T_END <= value')
    p_report.add_argument('--no-redundant', action='store_true', help='Keep only smallest grid size per family + multi-mesh')
    p_report.add_argument('--max-gold-time', type=float, help='Skip cases where gold took longer than N seconds')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    # Parse all cases
    print("Parsing FDS_Cases.sh...", end=' ', flush=True)
    all_cases = parse_fds_cases()
    print(f"{len(all_cases)} cases found.")

    cases = filter_cases(all_cases, args)
    if len(cases) != len(all_cases):
        print(f"After filtering: {len(cases)} cases selected.")

    # Dispatch command
    if args.command == 'discover':
        return cmd_discover(args, cases)
    elif args.command == 'generate-gold':
        return cmd_generate_gold(args, cases)
    elif args.command == 'test':
        return cmd_test(args, cases)
    elif args.command == 'report':
        return cmd_report(args, cases)

    return 0


if __name__ == '__main__':
    sys.exit(main())
