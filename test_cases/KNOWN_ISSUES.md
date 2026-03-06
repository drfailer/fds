# Known Issues with Test System

## Graph Termination Hang (March 6, 2026)

**Issue**: The Hedgehog graph sometimes hangs after completing all timesteps and doesn't properly terminate.

**Symptoms**:
- FDS completes all timesteps (e.g., Time Step: 27, Simulation Time: 0.1000000 s)
- Output files are partially written (HRR complete, DEVC empty or partial)
- Process never prints "[FDS-HH] Graph terminated" or "STOP: FDS completed successfully"
- Process must be killed with timeout or manual intervention

**Current Status**:
- This affects the `run_tests.py --generate-gold` functionality
- Tests run manually but hang at completion
- Workaround: Use existing known-good outputs as gold files

**Root Cause** (suspected):
- TimestepLoopState may not be properly emitting final output tokens
- Graph termination detection may be waiting for something that never arrives
- Possible issue with finishPushingData() or waitForTermination() logic

**Workaround**:
1. Use gold files generated from known-good test runs
2. For manual testing, kill process after simulation completes
3. CSV output files are valid even if process hangs

**TODO**:
- Debug TimestepLoopState termination logic
- Check if done flag is being set correctly in BarrierData
- Verify all graph outputs are properly wired
- Add timeout-based success detection as fallback

## Impact

- Gold file generation requires manual intervention
- Automated testing works (comparison only, not generation)
- Manual test execution works with timeout

## Testing Without Auto-Generation

Current workaround for testing:
1. Copy known-good CSV outputs to `gold/test_name/` directory
2. Run `./run_tests.py --test test_name` for comparison
3. Comparison works correctly even though generation doesn't
