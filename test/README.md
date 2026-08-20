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
   - end-to-end tests of import pipelines, staged-data audits, report and
     figure builders, demo scripts, and tests that run the real model or
     need locally installed scientific archives
   - these run only when explicitly targeted (`runtests` on a file or on
     `test/regression`); no default runner discovers this folder
5. `unit/`
   - ordinary unit tests intended for default discovery
   - fast, function-level contract tests only: kernels, helpers,
     validators, and schema gates that everyday development must rerun
   - a test belongs in `regression/`, not here, when it stages datasets
     end to end, renders a report or figures, runs the model, or needs
     local archive data
6. `benchmarks/`
   - component benchmarks and selected exploratory microbenchmarks
   - top-level benchmark files are the core kernel benchmarks run by default
   - opt-in microbenchmarks can live in subfolders such as `benchmarks/micro/`
   - the core benchmark files are `SebKernelPerfTest.m`,
     `ColumnKernelPerfTest.m`, and `SpectralKernelPerfTest.m`
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
5. rolling/default comparisons use the official gap-filled PROMICE forcing
   (`forcings = "promice_filled"`), with each case's `sitename` selecting the
   station-specific artifact
6. explicit `baseline = "v1.1"` comparisons retain the frozen release's
   historical station forcing (`kanm` or `kanl`) under the existing case ids;
   any other release tag must register its forcing identity before use

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
   - both matrix helpers accept `baseline`; blank/rolling selects
     `promice_filled`, while registered frozen releases retain their own
     forcing identity
3. `icemodel.test.helpers.setModelOptsForCase(...)`
   - canonical builder for the model `opts` used by one regression case
4. `icemodel.test.helpers.runModelCase(...)`
   - shared case setup, model dispatch, and canonical postprocessing path used
     by regression and verification investigations
5. `icemodel.test.helpers.getFormalTestSuiteCases()`
   - canonical ordered regression/bootstrap cases used by
     `run_test_bootstrap(...)`

## Which tool to use

1. `run_test_bootstrap(...)`
   - First-time setup or full refresh / pre-release orchestration entry point.
   - The only tool that cleans up and backs up `test/artifacts` and the
     managed perf/regression baseline files.
   - Cleanup removes mutable rolling MAT files only. Registered immutable
     releases such as v1.1 are loaded for release checks, not snapshotted from
     a rolling baseline with a different forcing identity.
2. `build_regression_baseline(...)` and `build_perf_baseline(...)`
   - Rebuild/accept new rolling or versioned baselines.
   - Writes baseline files; does not produce compare artifacts.
   - Direct versioned builds cannot replace an existing release file.
3. `snapshot_regression_baseline(...)` and `snapshot_perf_baseline(...)`
   - Freeze a release regression/perf baseline from the current rolling
     regression/perf baseline.
   - Existing release files are immutable, including when `overwrite=true` is
     supplied. The rolling source's saved forcing must match the registered
     current rolling identity and the registered release identity before a new
     snapshot is written.
4. `run_regression_suite(...)` and `run_perf_suite(...)`
   - Compare against existing rolling or release baselines
   - Does not mutate baselines
   - Blank `data_root` selects the baseline registration's canonical tree:
     rolling uses the verification tree, while frozen v1.1 uses historical
     `test/data`. An explicit root remains authoritative.
   - The baseline selector also selects the matching formal forcing identity:
     rolling uses `promice_filled`; frozen v1.1 uses its historical station
     aliases. Unknown release identities fail before model dispatch.
   - A rolling file that predates `promice_filled` stops with
     `rollingBaselineForcingAcceptanceRequired`; accept the authorized rolling
     baseline before using it for comparison.
   - Writes artifacts under `test/artifacts/<yyyymmdd-HHMMSS>/`.
   - Renders a self-contained Quarto HTML report with plots and a compact CSV
     in the same directory; use `build_report=false` only for an artifact-only
     diagnostic run.
   - `run_perf_suite` also runs the core benchmark suite and stores
     benchmark timing comparison alongside the formal perf artifact.
   - Formal performance runs force the MATLAB profiler off. Whole-model
     runtimes must stay inside the accepted two-sided tolerance band so an
     unexplained speedup cannot hide an inflated or incomplete reference.
   - Measurement protocol. Formal timings are sensitive to session
     history (JIT state, persistents, heap layout), so `run_perf_suite`
     enforces four controls:
     1. `isolation="process"` (the default) runs every case in a fresh
        `matlab -batch` subprocess. Use this mode for every formal
        accept/reject verdict on a refactor. The per-case spec and result
        MAT files land in the run's artifact folder. The opt-in
        `isolation="session"` times cases in the current session (with
        `clear functions` hygiene before each case) and is for quick
        diagnostics only.
     2. Every suite runner records itself in the session-activity record
        (`icemodel.test.helpers.markTestSessionDirty`). An in-session
        formal perf run REFUSES to start in a session that already ran
        another suite (`icemodel:test:perf:contaminatedSession`); use
        process isolation or a fresh session. `matlab -batch` one-shot
        runs always start clean.
     3. Case order is randomized (the seed rides the artifact), and each
        case's samples pass a dispersion validity gate
        (`max/median <= 1.5`). An invalid sample set re-measures once,
        then fails as "measurement invalid" — never a phantom verdict.
     4. An ambient anchor re-measures the first executed case at the end
        of the run. The dispersion gate cannot see load or scheduling
        shifts that are steady within each case but different across
        cases; an anchor drift above 15 percent marks every verdict in
        the run ambient-invalid (`meta.ambient_stable = false`).
   - A/A acceptance for the protocol: two consecutive
     `isolation="process"` runs of the same commit must pass the
     tolerance band against each other on all rows.
   - Isolation joins the baseline-compatibility check: timings compare
     against a baseline only when both used the same isolation protocol
     (a baseline without the field counts as "session"). Measured on
     this host, process-isolated cases run a systematic ~25 percent
     slower than session-shared-state cases, so cross-protocol
     comparison would produce phantom verdicts. Refactor gating in
     process mode therefore compares two isolated runs (before vs
     after) rather than the session-built rolling baseline.
5. `run_unit_suite(...)`
   - Use for folder-based unit-test discovery under `test/unit/`.
   - Use `debug=true` to stop on first failure for inspection.
   - The suite runs one test file at a time. A timestamped progress line
     prints before and after each file (to stdout in a desktop session,
     to stderr otherwise so a redirected `matlab -batch` run stays
     observable), and a per-file wall-clock table prints at the end,
     slowest first.
   - Use `progress_log="<path>"` to also append each progress line to
     that file with a per-line open/write/close. A hung run then leaves
     a log whose last `...` line names the file that never finished.
   - Expected full-suite wall-clock: roughly 8 minutes for the 63-file
     `test/unit` tree (measured 2026-08-17: 474.8 s, 1155 tests, under
     `matlab -nodisplay -nosplash -batch`; the small runner-contract
     test file landed after that measurement). No single file takes more
     than about 45 s; check the progress log before assuming a hang.
     The end-to-end staging, report, and model-run tests that used to
     dominate the runtime live in `test/regression` and run only when
     explicitly targeted.
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
   - The rename/round benchmarks live in `RenameRoundTest.m`.
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
9.  `run_promice_ablation_evaluation(...)`
   - Audits all canonical PROMICE case-years and runs only an explicit
     `case_ids` selection (or `"all"`) against the pinned `promice_filled`
     artifact.
   - Initializes each selected year on January 1 and saves the June 1--October
     1 snow-aware comparison payload: observed 600--900 kg m^-3 broad porous-
      weathering-crust sensitivity band, legacy melt, the diagnosed runoff
      proxy, refreezing, net physical solid loss from modeled phase and vapor
      terms, closure, and quantized top-layer deletion. That solid-loss series
      excludes remeshing and domain exchange. The 600 kg m^-3 density endpoint
      is not intact glacier-ice density, and the numeric band edges are
      ordered pointwise for signed lowering.
   - Chooses one longest summer interval with a 0.05 m trace-snow continuity
     threshold and a 0.01 m exposed-ice threshold for paired values. Censored
     rows stay explicit: the interval is not split, and snow-covered rows are
     not scored.
   - An empty `case_ids` selection is readiness-only and writes no run
     artifacts. To persist the ledger, summaries, and saved result bundle,
     name the cases (or pass `case_ids="all"`) and set
     `write_artifacts=true`. Setting `write_artifacts=true` with an empty
     selection is an error.
   - Render the saved bundle with
     `icemodel.verification.report.buildAblationEvaluationReport(...)`; the
     report layer never reruns the model or rereads canonical science inputs.
     Its cross-site figure is limited to completed selected site-years, while
     operational readiness and run accounting remain available as appendix
     prose and downloadable CSV files.
   - On macOS 26.5.2 (Mac16,12) with MATLAB R2025b, this benchmark measured
     elapsed time from artifact-folder creation through `results.mat`. One
     site-year took about two minutes. Nine one-site runs gave a median of 96
     seconds and a range of 71 to 107 seconds. The two-site run took 169 seconds.
     Full-cohort runs took 3,050
     seconds and 4,851 seconds, or 51 to 81 minutes, for 118 admitted rows.
     Allow 90 minutes for a conservative full-cohort estimate; these are
     operational estimates, not performance gates.

## Execution policy

Formal suites run from `icemodel` only and read these local files:

1. baselines from `test/baselines/`
2. references from `test/references/`

They do not require `runoff` on path at execution time.

The SUMup canonical identity-union regression is excluded from an ordinary
regression pass by an assumption because it restages 47 observation
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
