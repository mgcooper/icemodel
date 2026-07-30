# Test Suite

This folder contains the public test runners, regression data, unit tests,
and **component** benchmark material for the public `icemodel` repo.

Operator-facing usage notes for the public runners and study tools live in:

- `/Users/mattcooper/MATLAB/projects/icemodel/test/TOOL_REFERENCE.md`

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

Formal runtime contract for the default regression/perf matrices:

1. the case matrix carries one canonical retained `simyear`
2. the runtime contract expands that to `[simyear - 1, simyear]`
3. `n_spinup_years = 1`
4. `output_years = simyear`

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
   - First-time setup or full refresh / pre-release orchestration entry point.
   - The only tool that owns cleanup/backups of `test/artifacts` and
     managed perf/regression baseline files.
2. `build_regression_baseline(...)` and `build_perf_baseline(...)`
   - Rebuild/accept new rolling or versioned baselines.
   - Writes baseline files; does not produce compare artifacts.
3. `snapshot_regression_baseline(...)` and `snapshot_perf_baseline(...)`
   - Freeze a release regression/perf baseline from the current rolling
     regression/perf baseline.
4. `run_regression_suite(...)` and `run_perf_suite(...)`
   - Compare against existing rolling or release baselines
   - Does not mutate baselines
    - Writes artifacts under `test/artifacts/<yyyymmdd-HHMMSS>/`.
    - Renders a self-contained Quarto HTML report with plots and a compact CSV
      in the same directory; use `build_report=false` only for an artifact-only
      diagnostic run.
   - `run_perf_suite` also runs the core benchmark suite and stores
     benchmark timing comparison alongside the formal perf artifact.
   - Formal performance runs force the MATLAB profiler off. Whole-model
     runtimes must stay inside the accepted two-sided tolerance band so an
     unexplained speedup cannot hide an inflated or incomplete reference.
5. `run_unit_suite(...)`
   - Use for folder-based unit-test discovery under `test/unit/`.
   - Use `debug=true` to stop on first failure for inspection.
6. `run_benchmark_suite(...)`
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
7.  `build_runoff_reference_from_runoff(...)`
   - Refresh the static runoff reference data in `test/references/`.
   - This is separate from baseline management and requires the sibling
     `runoff` project.
8.  `validate_test_suite(...)`
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

The SUMup canonical identity-union regression is deliberately excluded from an
ordinary regression pass by an assumption because it restages 47 observation
cases from three multi-million-row NetCDF files. Before a canonical SUMup
replacement, opt in explicitly and run only that file:

```matlab
addpath("icemodel")
setenv("ICEMODEL_RUN_SUMUP_IDENTITY_UNION", "1")
results = runtests("test/regression/test_sumup_identity_union.m");
setenv("ICEMODEL_RUN_SUMUP_IDENTITY_UNION", "")
assertSuccess(results)
```

The procedure uses a temporary root, pins all three source SHA-256 values,
checks the 47-case plus external-Humphrey identity union and per-case exclusive
contribution, verifies per-artifact duplicate counts, and runs artifact QA. It
does not replace the canonical tree.

Rolling baseline rebuilds automatically archive the prior managed MAT file
and any saved profiler artifacts under `test/baselines/archive/` before the
new rolling baseline is written.
