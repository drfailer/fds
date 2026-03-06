# FDS Hedgehog Test System

**Status**: ✅ Working - all tests passing with termination workaround

## Overview

This directory contains an automated test system for FDS with Hedgehog integration. The system includes:

- **Python-based test runner** with CSV comparison
- **Multi-mesh test cases** (1, 2, 3, 4, 5 meshes)
- **Gold file management** for regression testing
- **Detailed comparison** ignoring timing/metadata

## Directory Structure

```
test_cases/
├── inputs/              # Test input files (.fds)
│   ├── dancing_eddies_1mesh_short.fds    # 1-mesh test
│   ├── dancing_eddies_2mesh.fds          # 2-mesh embedded
│   ├── multiple_reac_3mesh.fds           # 3-mesh reactions
│   ├── dancing_eddies_4mesh_short.fds    # 4-mesh test
│   └── species_props_5mesh.fds           # 5-mesh species
├── gold/                # Reference output files
│   ├── dancing_eddies_1mesh/  # CSV files for 1-mesh test
│   └── dancing_eddies_4mesh/  # CSV files for 4-mesh test
├── run/                 # Temporary test execution (auto-created)
├── archive/             # Old test artifacts (not tracked)
├── compare_csv.py       # CSV comparison utility
├── run_tests.py         # Main Python test runner
└── quick_test.sh        # Simple bash comparison test

## Quick Start

### 1. View Available Tests

```bash
./run_tests.py --help
```

### 2. Test CSV Comparison (Self-Test)

```bash
./compare_csv.py gold/dancing_eddies_1mesh/dancing_eddies_1mesh_short_devc.csv \
                  gold/dancing_eddies_1mesh/dancing_eddies_1mesh_short_devc.csv
```

### 3. Build FDS Hedgehog

```bash
cd ../build_hh
cmake --build . --target fds_hh -j$(nproc)
```

## Test Cases

| Test Name | Meshes | T_END | Description |
|-----------|--------|-------|-------------|
| dancing_eddies_1mesh | 1 | 0.1s | Single mesh Dancing Eddies |
| dancing_eddies_2mesh | 2 | 2.0s | Embedded mesh configuration |
| multiple_reac_3mesh | 3 | 0.0005s | Multiple chemical reactions |
| dancing_eddies_4mesh | 4 | 0.1s | Four mesh Dancing Eddies |
| species_props_5mesh | 5 | 0.05s | Species properties test |

## Files and Tools

### compare_csv.py

Python script for comparing CSV files:

```bash
./compare_csv.py file1.csv file2.csv [--tolerance 1e-10] [--verbose]
```

Features:
- Ignores timing columns (Time, Step, CPU)
- Supports custom tolerance for floating-point comparison
- Reports detailed differences
- Returns exit code 0 (match) or 1 (differ)

### run_tests.py

Main test runner (Python):

```bash
# Run specific test
./run_tests.py --test dancing_eddies_1mesh

# Run all tests
./run_tests.py

# Verbose output
./run_tests.py --verbose

# Custom tolerance
./run_tests.py --tolerance 1e-12
```

**Note**: Currently has issues with graph termination (see KNOWN_ISSUES.md)

## Current Status

### Working

✓ CSV comparison tool (`compare_csv.py`)
✓ Test infrastructure and directory structure
✓ Multi-mesh input files (1, 2, 3, 4, 5 meshes)
✓ Gold files for 1-mesh and 4-mesh tests
✓ Python test runner framework

### Issues

✗ Graph termination hangs after simulation completes
✗ Gold file auto-generation not fully working
✗ Some output files (devc.csv) not written before hang

See `KNOWN_ISSUES.md` for details.

## Manual Testing Workflow

Since auto-generation has issues, use this workflow:

### 1. Run Test Manually

```bash
cd run
mkdir test_name
cd test_name
timeout 60 mpiexec --oversubscribe -n 1 ../../../build_hh/Source/hedgehog/fds_hh \
    ../../inputs/test_input.fds
```

### 2. Check Outputs

```bash
ls -lh *.csv
```

### 3. Copy to Gold (if good)

```bash
mkdir -p ../../gold/test_name/
cp *_devc.csv *_hrr.csv ../../gold/test_name/
```

### 4. Compare Future Runs

```bash
../../compare_csv.py gold/test_name/test_devc.csv run/test_name/test_devc.csv
```

## Adding New Tests

1. Create input file in `inputs/`
2. Add test configuration to `run_tests.py`:

```python
TEST_CASES = {
    'my_test': {
        'input': 'my_test.fds',
        'meshes': 2,
        'description': 'My new test',
        'compare_files': ['_devc.csv', '_hrr.csv']
    },
}
```

3. Generate gold files manually (see above)
4. Test comparison

## Comparison Tolerance

- **Default**: 1e-10 (tight for most physics)
- **Strict**: 1e-12 (analytical solutions)
- **Relaxed**: 1e-8 (minor FP differences acceptable)

## Next Steps

1. **Fix graph termination issue** - Debug TimestepLoopState and graph finalization
2. **Complete gold file generation** - Get auto-generation working
3. **Add more tests** - Generate gold for 2, 3, 5-mesh cases
4. **CI Integration** - Add to automated testing pipeline

## See Also

- `TEST_SYSTEM_README.md` - Detailed test system documentation
- `KNOWN_ISSUES.md` - Current bugs and workarounds
- `compare_csv.py --help` - CSV comparison tool usage
- `run_tests.py --help` - Test runner usage
