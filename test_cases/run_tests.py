#!/usr/bin/env python3
"""
FDS Hedgehog Test Runner

Runs test cases and compares outputs with gold files.
Can generate gold files from master branch if needed.
"""

import os
import sys
import subprocess
import argparse
import json
import select
from pathlib import Path
from typing import Dict, List, Tuple
import time

# Paths
REPO_ROOT = Path(__file__).parent.parent
TEST_DIR = Path(__file__).parent
INPUTS_DIR = TEST_DIR / "inputs"
GOLD_DIR = TEST_DIR / "gold"
RUN_DIR = TEST_DIR / "run"
BUILD_DIR = REPO_ROOT / "build_hh"
FDS_HH = BUILD_DIR / "Source" / "hedgehog" / "fds_hh"
FDS_ORIG = REPO_ROOT / "Build" / "ompi_gnu_linux_db" / "fds_ompi_gnu_linux_db"
COMPARE_SCRIPT = TEST_DIR / "compare_csv.py"

# Test cases configuration
TEST_CASES = {
    'dancing_eddies_1mesh': {
        'input': 'dancing_eddies_1mesh_short.fds',
        'meshes': 1,
        'description': '1-mesh Dancing Eddies',
        'compare_files': ['_devc.csv', '_hrr.csv']
    },
    'dancing_eddies_2mesh': {
        'input': 'dancing_eddies_2mesh.fds',
        'meshes': 2,
        'description': '2-mesh Dancing Eddies (embedded)',
        'compare_files': ['_devc.csv']
    },
    'multiple_reac_3mesh': {
        'input': 'multiple_reac_3mesh.fds',
        'meshes': 3,
        'description': '3-mesh Multiple Reactions',
        'compare_files': ['_devc.csv']
    },
    'dancing_eddies_4mesh': {
        'input': 'dancing_eddies_4mesh_short.fds',
        'meshes': 4,
        'description': '4-mesh Dancing Eddies',
        'compare_files': ['_devc.csv', '_hrr.csv']
    },
    'species_props_5mesh': {
        'input': 'species_props_5mesh.fds',
        'meshes': 5,
        'description': '5-mesh Species Properties',
        'compare_files': ['_devc.csv']
    },
}

class TestRunner:
    def __init__(self, verbose: bool = False, tolerance: float = 1e-10):
        self.verbose = verbose
        self.tolerance = tolerance
        self.results = {}

    def log(self, message: str, level: str = "INFO"):
        """Print log message."""
        prefix = {
            "INFO": "  ",
            "PASS": "✓ ",
            "FAIL": "✗ ",
            "WARN": "⚠ ",
            "RUN": "→ "
        }.get(level, "  ")
        print(f"{prefix}{message}")

    def ensure_directories(self):
        """Create necessary directories."""
        RUN_DIR.mkdir(exist_ok=True)
        GOLD_DIR.mkdir(exist_ok=True)
        self.log(f"Work directory: {RUN_DIR}")
        self.log(f"Gold directory: {GOLD_DIR}")

    def check_executable(self, exe_path: Path, name: str) -> bool:
        """Check if executable exists."""
        if not exe_path.exists():
            self.log(f"{name} not found at {exe_path}", "FAIL")
            return False
        if not os.access(exe_path, os.X_OK):
            self.log(f"{name} is not executable: {exe_path}", "FAIL")
            return False
        return True

    def run_fds(self, input_file: Path, exe: Path, work_dir: Path, timeout: int = 60) -> Tuple[bool, float]:
        """
        Run FDS simulation.

        For fds_hh: Monitors output and terminates process as soon as dot file is written
        (which happens right before the waitForTermination() hang).

        Returns:
            (success, elapsed_time)
        """
        chid = input_file.stem
        input_path = INPUTS_DIR / input_file

        if not input_path.exists():
            self.log(f"Input file not found: {input_path}", "FAIL")
            return False, 0.0

        # Copy input to work directory
        work_input = work_dir / input_file.name
        subprocess.run(['cp', str(input_path), str(work_input)], check=True)

        # Run FDS
        cmd = ['mpiexec', '--oversubscribe', '-n', '1', str(exe), input_file.name]

        try:
            start_time = time.time()

            # Use Popen to monitor output in real-time and kill when complete
            process = subprocess.Popen(
                cmd,
                cwd=work_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1
            )

            stdout_lines = []
            completion_detected = False

            # Monitor stdout for completion signal
            while True:
                # Check if process is still running
                if process.poll() is not None:
                    break

                # Read available output with timeout
                ready, _, _ = select.select([process.stdout], [], [], 0.1)
                if ready:
                    line = process.stdout.readline()
                    if line:
                        stdout_lines.append(line)
                        if self.verbose:
                            print(line, end='')

                        # Detect completion signal (dot file written = simulation complete)
                        if "Graph dot file written" in line:
                            completion_detected = True
                            # Give it a moment to finish writing, then kill
                            time.sleep(0.5)
                            process.terminate()
                            try:
                                process.wait(timeout=2)
                            except subprocess.TimeoutExpired:
                                process.kill()
                                process.wait()
                            break

                # Safety timeout
                elapsed = time.time() - start_time
                if elapsed > timeout:
                    process.kill()
                    process.wait()
                    break

            elapsed = time.time() - start_time
            stdout = ''.join(stdout_lines)

            # Check for success
            out_file = work_dir / f"{chid}.out"
            if completion_detected or out_file.exists():
                with open(out_file, 'r') as f:
                    content = f.read()
                    if "Total Time:" in content:
                        return True, elapsed

            self.log(f"FDS did not complete successfully", "FAIL")
            return False, elapsed

        except Exception as e:
            self.log(f"Error running FDS: {e}", "FAIL")
            import traceback
            if self.verbose:
                traceback.print_exc()
            return False, 0.0

    def compare_files(self, chid: str, work_dir: Path, gold_dir: Path, compare_files: List[str]) -> Tuple[bool, Dict]:
        """
        Compare output files with gold files.

        Returns:
            (all_passed, comparison_results)
        """
        results = {}
        all_passed = True

        for suffix in compare_files:
            test_file = work_dir / f"{chid}{suffix}"
            gold_file = gold_dir / f"{chid}{suffix}"

            if not test_file.exists():
                self.log(f"Output file missing: {test_file.name}", "WARN")
                results[suffix] = {'status': 'missing', 'passed': False}
                all_passed = False
                continue

            if not gold_file.exists():
                self.log(f"Gold file missing: {gold_file.name}", "WARN")
                results[suffix] = {'status': 'no_gold', 'passed': False}
                all_passed = False
                continue

            # Run comparison
            cmd = [
                sys.executable,
                str(COMPARE_SCRIPT),
                str(test_file),
                str(gold_file),
                '--tolerance', str(self.tolerance)
            ]

            try:
                result = subprocess.run(cmd, capture_output=True, text=True)
                passed = (result.returncode == 0)

                results[suffix] = {
                    'status': 'compared',
                    'passed': passed,
                    'output': result.stdout
                }

                if not passed:
                    all_passed = False
                    if self.verbose:
                        self.log(result.stdout)

            except Exception as e:
                self.log(f"Comparison error for {suffix}: {e}", "FAIL")
                results[suffix] = {'status': 'error', 'passed': False}
                all_passed = False

        return all_passed, results

    def run_test(self, test_name: str, test_config: Dict) -> Dict:
        """Run a single test case."""
        self.log(f"\n{'='*60}")
        self.log(f"Test: {test_name}")
        self.log(f"Description: {test_config['description']}")
        self.log(f"Meshes: {test_config['meshes']}")
        self.log(f"{'='*60}")

        input_file = Path(test_config['input'])
        chid = input_file.stem
        work_dir = RUN_DIR / test_name
        work_dir.mkdir(exist_ok=True)

        result = {
            'test_name': test_name,
            'description': test_config['description'],
            'meshes': test_config['meshes'],
            'passed': False,
            'run_time': 0.0,
            'comparison': {}
        }

        # Run FDS
        self.log(f"Running FDS...", "RUN")
        success, elapsed = self.run_fds(input_file, FDS_HH, work_dir)
        result['run_time'] = elapsed

        if not success:
            self.log(f"FDS run failed ({elapsed:.2f}s)", "FAIL")
            return result

        self.log(f"FDS completed ({elapsed:.2f}s)", "PASS")

        # Compare with gold
        self.log(f"Comparing outputs...", "RUN")
        gold_subdir = GOLD_DIR / test_name
        all_passed, comparison = self.compare_files(chid, work_dir, gold_subdir, test_config['compare_files'])
        result['comparison'] = comparison

        if all_passed:
            result['passed'] = True
            self.log(f"All comparisons passed", "PASS")
        else:
            self.log(f"Some comparisons failed", "FAIL")

        return result

    def generate_gold(self, test_name: str, test_config: Dict, use_original: bool = True) -> bool:
        """Generate gold files for a test case."""
        self.log(f"\n{'='*60}")
        self.log(f"Generating gold files for: {test_name}")
        self.log(f"{'='*60}")

        input_file = Path(test_config['input'])
        chid = input_file.stem
        gold_subdir = GOLD_DIR / test_name
        gold_subdir.mkdir(exist_ok=True)

        # Choose which FDS to use
        if use_original:
            if not self.check_executable(FDS_ORIG, "FDS original"):
                self.log("Original FDS not available, using fds_hh", "WARN")
                exe = FDS_HH
            else:
                exe = FDS_ORIG
        else:
            exe = FDS_HH

        self.log(f"Using: {exe.name}", "INFO")

        # Run FDS in gold directory
        self.log(f"Running FDS...", "RUN")
        success, elapsed = self.run_fds(input_file, exe, gold_subdir)

        if not success:
            self.log(f"Failed to generate gold files", "FAIL")
            return False

        self.log(f"Gold files generated ({elapsed:.2f}s)", "PASS")

        # Clean up non-essential files
        for pattern in ['*.sf', '*.sf.bnd', '*.smv', '*.out', '*.fds', '*_git.txt']:
            for f in gold_subdir.glob(pattern):
                f.unlink()
                if self.verbose:
                    self.log(f"Cleaned up: {f.name}")

        return True

    def print_summary(self):
        """Print test summary."""
        print(f"\n{'='*60}")
        print("TEST SUMMARY")
        print(f"{'='*60}")

        total = len(self.results)
        passed = sum(1 for r in self.results.values() if r['passed'])
        failed = total - passed

        for test_name, result in sorted(self.results.items()):
            status = "PASS" if result['passed'] else "FAIL"
            status_symbol = "✓" if result['passed'] else "✗"
            print(f"{status_symbol} {test_name:30s} {status:6s} ({result['run_time']:.2f}s)")

        print(f"\n{'='*60}")
        print(f"Total: {total}  Passed: {passed}  Failed: {failed}")
        print(f"{'='*60}")

        return failed == 0

def main():
    parser = argparse.ArgumentParser(description='Run FDS Hedgehog tests')
    parser.add_argument('--generate-gold', action='store_true',
                        help='Generate gold files instead of running tests')
    parser.add_argument('--use-hedgehog-gold', action='store_true',
                        help='Use fds_hh instead of original FDS for gold generation')
    parser.add_argument('--test', '-t', action='append',
                        help='Run specific test(s) only')
    parser.add_argument('--tolerance', type=float, default=1e-10,
                        help='Comparison tolerance (default: 1e-10)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Verbose output')

    args = parser.parse_args()

    runner = TestRunner(verbose=args.verbose, tolerance=args.tolerance)
    runner.ensure_directories()

    # Select tests to run
    if args.test:
        tests = {k: v for k, v in TEST_CASES.items() if k in args.test}
        if not tests:
            print(f"No matching tests found. Available tests: {', '.join(TEST_CASES.keys())}")
            return 1
    else:
        tests = TEST_CASES

    if args.generate_gold:
        # Generate gold files
        print("Generating gold files...")
        use_original = not args.use_hedgehog_gold

        for test_name, test_config in tests.items():
            success = runner.generate_gold(test_name, test_config, use_original)
            if not success:
                print(f"Failed to generate gold for {test_name}")
                return 1

        print("\nGold file generation complete!")
        return 0
    else:
        # Run tests
        print("Running tests...")

        # Check executable
        if not runner.check_executable(FDS_HH, "fds_hh"):
            print("\nPlease build fds_hh first:")
            print("  cd build_hh && cmake --build . --target fds_hh")
            return 1

        # Run each test
        for test_name, test_config in tests.items():
            result = runner.run_test(test_name, test_config)
            runner.results[test_name] = result

        # Print summary
        all_passed = runner.print_summary()

        return 0 if all_passed else 1

if __name__ == '__main__':
    sys.exit(main())
