# Test Suite

This folder contains the public test runners, regression data, unit tests,
and **component** benchmark material for the public `icemodel` repo.

## Layout

1. `artifacts/`
   - compare-run outputs grouped by batch run under
     `test/artifacts/<yyyymmdd-HHMMSS>/`
2. `baselines/`
   - mutable rolling baselines and frozen release baselines for perf and
     regression
3. `references/`
   - static external reference data such as `runoff_reference.mat`
4. `regression/`
   - software-level regression classes, including performance regression
5. `unit/`
   - ordinary unit tests intended for default discovery
6. `benchmarks/`
   - component benchmarks and selected exploratory microbenchmarks
   - top-level benchmark files are the core kernel benchmarks run by default
   - opt-in microbenchmarks can live in subfolders such as `benchmarks/micro/`
   - representative core benchmark files now live in:
     `SebKernelPerfTest.m`, `ColumnKernelPerfTest.m`,
     and `SpectralKernelPerfTest.m`
7. `tools/`
   - explicit build/snapshot utilities
8. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/+icemodel/+test/+helpers/`
   - shared helper functions under the `icemodel.test.helpers.*` namespace

## Files

### References

1. `references/runoff_reference.mat`

This file stores external runoff comparison context and does not represent
accepted model output.

### Perf baselines

1. `baselines/perf_baseline_<simyear>_rolling_<smbmodel>.mat`
2. `baselines/perf_baseline_<simyear>_<version>_<smbmodel>.mat`

### Regression baselines

1. `baselines/regression_baseline_rolling_<smbmodel>.mat`
2. `baselines/regression_baseline_<version>_<smbmodel>.mat`

## Regression Matrix

Default software-level regression coverage:

1. `icemodel`, `kanm`, `2016`, `solver = 1, 2, 3`
2. `icemodel`, `kanl`, `2016`, `solver = 1, 2, 3`
3. `skinmodel`, `kanm`, `2016`, `solver = 1`
4. `skinmodel`, `kanl`, `2016`, `solver = 1`
5. self-forced station runs only (`sitename == forcings`)

Programmatic regression helpers:

1. `icemodel.test.helpers.getPerfCaseMatrix(...)`
   - canonical performance-regression case matrix
2. `icemodel.test.helpers.getRegressionCaseMatrix(...)`
   - canonical numerical-regression case matrix
3. `icemodel.test.helpers.setModelOptsForCase(...)`
   - canonical builder for the model `opts` used by one regression case
4. `icemodel.test.helpers.getFormalTestSuiteCases()`
   - canonical ordered regression/bootstrap cases used by
     `run_test_bootstrap(...)`

## Which tool to use

1. `run_test_bootstrap(...)`
   - Use for first-time setup or full refreshes.
   - This is the only tool that owns cleanup/backups of `test/artifacts` and
     managed perf/regression baseline files.
2. `build_regression_baseline(...)`
   - Use to accept new rolling or versioned model-output baselines.
   - This writes baseline files and does not produce compare artifacts.
3. `build_perf_baseline(...)`
   - Use to accept new rolling or versioned runtime baselines.
   - This writes baseline files and does not produce compare artifacts.
4. `snapshot_regression_baseline(...)`
   - Use to freeze a release regression baseline from the current rolling
     regression baseline.
5. `snapshot_perf_baseline(...)`
   - Use to freeze a release perf baseline from the current rolling perf
     baseline.
6. `run_regression_suite(...)`
   - Use for normal compare runs against rolling or release regression
     baselines.
   - This writes artifacts under `test/artifacts/<yyyymmdd-HHMMSS>/`.
7. `run_perf_suite(...)`
   - Use for normal compare runs against rolling or release perf baselines.
   - This writes artifacts under `test/artifacts/<yyyymmdd-HHMMSS>/`.
   - By default it also runs the managed core benchmark suite and stores the
     benchmark timing comparison alongside the formal perf artifact.
8. `run_unit_suite(...)`
   - Use for folder-based unit-test discovery under `test/unit/`.
   - Supports `debug=true` to stop on first failure for inspection.
9. `run_benchmark_suite(...)`
   - Use for the formal benchmark suite under `test/benchmarks/`.
   - This remains the standalone component-benchmark runner.
   - By default it runs only the top-level benchmark files.
   - Use `include_subfolders=true` to opt into nested microbenchmarks.
   - Use `sampling_profile="fast|default|strict"` for common sampling
     budgets, or override the numeric runner controls directly when needed.
   - The benchmark suite is intended to explain where runtime is spent, not
     just compare alternative implementations in isolation.
   - Sampling-error warnings from the MATLAB perf framework are not test
     failures; they indicate that a microbenchmark stayed noisy at the
     current sampling budget even though the benchmark itself remained valid.
   - Benchmark-specific interpretation notes should live with the benchmark
     file itself when the timing result motivated a code choice.
   - The rename/round history is reconciled into `RenameRoundTest.m` rather
     than split across separate manual scripts.
10. `build_runoff_reference_from_runoff(...)`
   - Use to refresh the static runoff reference data in `test/references/`.
   - This is separate from baseline management and requires the sibling
     `runoff` project.
11. `validate_test_suite(...)`
   - Use to exercise the public test-suite surface end to end without
     mutating managed baselines.
   - This validates signatures, Code Analyzer cleanliness, runner selector
     variants, per-file discovery, and build/snapshot tools against
     temporary outputs.

## Execution policy

Formal suites run from `icemodel` only and read these local files:

1. baselines from `test/baselines/`
2. references from `test/references/`

They do not require `runoff` on path at execution time.

Rolling baseline rebuilds automatically archive the prior managed MAT file
and any saved profiler artifacts under `test/baselines/archive/` before the
new rolling baseline is written.
