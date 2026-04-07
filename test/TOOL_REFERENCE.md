# Test Tool Reference

This file documents the public test and study entrypoints under:

- `/Users/mattcooper/MATLAB/projects/icemodel/test/`
- `/Users/mattcooper/MATLAB/projects/icemodel/test/tools/`

Use this as the operator-facing reference for what to run, when to run it,
and which options matter.

## Main Test Runners

### `run_unit_suite`

Purpose:

- Run correctness/unit tests under `test/unit/`

Default use:

```matlab
clear functions
results = run_unit_suite();
```

Common use:

```matlab
clear functions
results = run_unit_suite(debug=true, stop_on_failure=true);
```

Use when:

- you changed core kernels/helpers and want quick correctness feedback

Key options:

- `selector`
  - optional file/folder/name filter
- `debug`
  - stop on first failure for interactive inspection
- `stop_on_failure`
  - stop after first failing test
- `verbosity`
  - command-window detail level

### `run_benchmark_suite`

Purpose:

- Run the standalone component benchmark suite under `test/benchmarks/`

Default use:

```matlab
clear functions
results = run_benchmark_suite();
```

Common use:

```matlab
clear functions
results = run_benchmark_suite( ...
    sampling_profile="fast", ...
    show_summary=true);
```

Use when:

- you want component/kernel timings
- you want to compare implementations outside the full model

Key options:

- `testname`
  - benchmark file/folder/name filter
- `include_subfolders`
  - opt into nested microbenchmarks
- `sampling_profile`
  - `fast`, `default`, `strict`
- `show_summary`
  - print compact timing table

Notes:

- This is not the formal end-to-end perf regression tool.
- Sampling-error warnings do not mean the benchmark is invalid.

### `run_regression_suite`

Purpose:

- Compare formal model outputs against accepted regression baselines

Default use:

```matlab
clear functions
results = run_regression_suite();
```

Common use:

```matlab
clear functions
results = run_regression_suite( ...
    tier="smoke", ...
    smbmodel="icemodel", ...
    solver=2, ...
    baseline="rolling");
```

Use when:

- you want to know whether modeled outputs changed relative to baseline

Key options:

- `tier`
  - `smoke`, `full`, `all`
- `smbmodel`
  - `icemodel`, `skinmodel`, `all`
- `solver`
  - optional subset for `icemodel`
- `simyear`
  - canonical retained output year
- `baseline`
  - usually `rolling`
- `smoke_sites`, `full_sites`
  - advanced site overrides

Important runtime contract:

- formal regression cases currently run with one leading spinup year
- for `simyear=2016`, the formal runtime contract is `[2015 2016]` with
  `n_spinup_years = 1`

### `run_perf_suite`

Purpose:

- Compare formal model runtime against accepted perf baselines

Default use:

```matlab
clear functions
results = run_perf_suite();
```

Common use:

```matlab
clear functions
results = run_perf_suite( ...
    tier="smoke", ...
    smbmodel="icemodel", ...
    solver=2, ...
    baseline="rolling", ...
    n_runs=1, ...
    include_benchmarks=false);
```

Use when:

- you want formal model-call timing vs accepted perf baselines

Key options:

- same case-selector options as `run_regression_suite`
- `n_runs`
  - fixed sample count for the perf harness
- `tol_perf`
  - allowed fractional slowdown gate when building baselines
- `include_benchmarks`
  - also run managed component benchmarks
- `benchmark_sampling_profile`
  - sampling budget for those managed benchmarks

Important note:

- the measured region in the formal perf class is the model call only
- it does not include report formatting, baseline loading, or runner overhead
- formal perf cases currently use the same canonical runtime contract as
  regression: for `simyear=2016`, the runtime contract is `[2015 2016]` with
  `n_spinup_years = 1`
- whole-model perf gating is skipped when the accepted perf baseline was built
  under a different MATLAB version/platform than the current run

### `run_test_bootstrap`

Purpose:

- Canonical full-suite orchestration entrypoint

Use when:

- you want the top-level full refresh / pre-release workflow

This is broader than the individual runners above.

## Baseline Management

### `build_regression_baseline`

Purpose:

- accept new regression outputs as the rolling or versioned regression
  baseline

Default use:

```matlab
clear functions
RegressionBaseline = build_regression_baseline(baseline="rolling");
```

Use when:

- a regression change is accepted

Notes:

- writes baseline files
- does not compare against an older baseline
- rolling builds archive the previous managed baseline first
- by default the rebuilt baselines use the formal 2-year contract:
  retained year plus one leading spinup year

### `build_perf_baseline`

Purpose:

- accept new perf measurements as the rolling or versioned perf baseline

Default use:

```matlab
clear functions
PerfBaseline = build_perf_baseline(baseline="rolling");
```

Use when:

- a runtime change is accepted

Notes:

- writes baseline files
- also stores managed benchmark baselines in the same perf MAT file
- by default the rebuilt baselines use the formal 2-year contract:
  retained year plus one leading spinup year

### `snapshot_regression_baseline`

Purpose:

- freeze a release regression baseline from the current rolling baseline

### `snapshot_perf_baseline`

Purpose:

- freeze a release perf baseline from the current rolling baseline

Use snapshots when:

- you want an immutable release baseline

## Spectral / Postprocess Study Tools

These are exploratory benchmark/study tools, not the core suite runners.

### `summarize_spectral_perf`

Purpose:

- produce one spectral comparison report

It reports:

1. kernel timings (inlined, exact, lookup — calls functions directly)
2. direct whole-model timings (exact vs lookup via `opts.lookup_k_bulk`)
3. output-agreement metrics

Default use:

```matlab
clear functions
report = summarize_spectral_perf( ...
    simyear=2016, ...
    n_direct_runs=1);
```

Kernel-only use:

```matlab
clear functions
report = summarize_spectral_perf( ...
    simyear=2016, ...
    include_full_model=false, ...
    include_direct_model=false);
```

Key options:

- `include_full_model`
  - master switch for the whole-model sections
- `include_direct_model`
  - include direct `tic/toc` whole-model timings
- `n_direct_runs`
  - number of repeated direct whole-model runs
- `output_file`
  - save one MAT report

Use when:

- you are comparing exact vs lookup spectral paths
- you want a focused spectral report instead of the full bootstrap workflow

### `run_spectral_study_bootstrap`

Purpose:

- orchestrate the whole spectral study workflow

This is the spectral-study analogue of `run_test_bootstrap(...)`, not of
`run_perf_suite(...)`.

Default use:

```matlab
clear functions
results = run_spectral_study_bootstrap();
```

Recommended diagnostic use:

```matlab
clear functions
results = run_spectral_study_bootstrap( ...
    include_spectral_perf=true, ...
    include_figures=true, ...
    include_density_floor=true);
```

Acceptance-only use:

```matlab
clear functions
results = run_spectral_study_bootstrap( ...
    include_spectral_perf=false, ...
    include_figures=false, ...
    include_density_floor=false, ...
    include_regression_acceptance=true, ...
    include_perf_acceptance=true);
```

What `include_spectral_perf=false` means:

- skip the `summarize_spectral_perf(...)` report entirely
- still allow bootstrap to run figures, density-floor audit, acceptance, or
  baseline rebuild steps

Why that option exists:

- sometimes you only want the orchestration around acceptance/baselines
- the direct spectral report can be the most expensive part of the workflow

What gets saved:

- the returned `results` struct is also saved to:
  - `results.bootstrap_file`
- subordinate reports are also saved separately when requested

### `summarize_spectral_density_floor`

Purpose:

- quantify how often the `rho >= 300` floor is active and how large the
  induced bulk-coefficient error can get

Use when:

- deciding whether the density floor is still justified

### `plot_spectral_variant_profiles`

Purpose:

- save comparison figures for bulk extinction coefficients, spectral net flux,
  spectral divergence, and thermal source term

Use when:

- you want to inspect shape differences, not just scalar errors

### `summarize_postprocess_perf`

Purpose:

- compare the fixed-step hourly aggregation path against legacy `retime(...)`

Use when:

- evaluating postprocess runtime only

This is intentionally separate from the spectral bootstrap.

### `audit_formal_substep_failures`

Purpose:

- run the formal regression cases and detect `dt_min` / `maxsubstep`
  fallback without mutating baselines

Default use:

```matlab
clear functions
report = audit_formal_substep_failures();
```

Common use:

```matlab
clear functions
report = audit_formal_substep_failures( ...
    tier="full", ...
    smbmodel="icemodel", ...
    stop_on_first_issue=true);
```

Use when:

- a baseline build or regression run prints `check_substep` fallback messages
- you want to know which exact formal case triggered them before rebuilding
  baselines

## Recommended Order

### Ordinary correctness work

1. `run_unit_suite(...)`
2. `run_benchmark_suite(...)` if the change is performance-sensitive
3. `run_regression_suite(...)` / `run_perf_suite(...)` if the change touches
   formal outputs or runtime

### Spectral study workflow

1. `run_spectral_study_bootstrap(... include_spectral_perf=true, ... )`
2. inspect:
   - kernel timings
   - direct-model timings
   - direct-model agreement
   - figures
   - density-floor audit
3. if any 2-year formal case emits `check_substep` fallback messages, run:
   - `audit_formal_substep_failures(...)`
4. `run_spectral_study_bootstrap(... include_regression_acceptance=true, include_perf_acceptance=true)`
5. only after that, rebuild rolling baselines if the change is accepted

### Baseline refresh / release

1. `build_regression_baseline(...)`
2. `build_perf_baseline(...)`
3. optional:
   - `snapshot_regression_baseline(...)`
   - `snapshot_perf_baseline(...)`

## Spectral Variant Interface

Production interface:

- `opts.lookup_k_bulk = true` (default) — lookup-table bulk extinction
- `opts.lookup_k_bulk = false` — exact bulk-extinction transform

Kernel benchmarks and study tools call the spectral functions directly and
retain all three historical variants (inlined, exact, lookup) for comparison.
The `inlined` variant lives in `test/legacy/SPECTRALSOURCETERM_INLINE.m`.
