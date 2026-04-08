# Spectral Performance Study

## Scope

This study compares three spectral source-term paths:

1. `inlined`
   - historical fully inlined implementation preserved in
     `/Users/mattcooper/MATLAB/projects/icemodel/test/legacy/SPECTRALSOURCETERM_INLINE.m`
2. `functions`
   - organized exact production path in
     `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/icemodel.column.shortwave_source_term.m`
3. `lookup`
   - same organized production path, but with bulk extinction coefficients
     gathered from a precomputed density lookup table

The study answers two separate questions:

1. Is organized helper-call code slower than the historical inlined code?
2. How much speed is available from replacing the exact bulk-extinction
   transform with a lookup table?

## Production Code Shape

`opts.lookup_k_bulk` (default `true`) controls whether the production path
uses the lookup-table approximation or the exact bulk-extinction transform.
`icemodel.radiation.initialize_spectral_model` builds the lookup table when `opts.lookup_k_bulk` is true and
passes it as `k_bulk_lookup` to `icemodel.column.shortwave_source_term`, which dispatches on
`isempty(k_bulk_lookup)`.

1. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/icemodel.radiation.initialize_spectral_model.m`
   - initializes spectral geometry and optical coefficients
   - builds `k_bulk_lookup` when `opts.lookup_k_bulk` is true
2. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/icemodel.column.shortwave_source_term.m`
   - main entrypoint called by `icemodel.m`
   - remaps thermal density to the spectral grid
   - chooses exact or lookup bulk extinction coefficients via `isempty(k_bulk_lookup)`
   - solves the two-stream system
   - reconstructs net spectral flux
   - collapses that flux to the thermal-grid source term and `chi`
3. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/icemodel.radiation.bulk_extinction_coefficients.m`
   - exact bulk-extinction transform
4. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/icemodel.radiation.bulk_extinction_coefficients_lookup.m`
   - lookup bulk-extinction transform
5. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/icemodel.radiation.solvetwostream.m`
   - two-stream solve on the spectral grid
6. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/SPECTRALNETFLUX.m`
   - reconstructs net spectral flux from the up/down solution
7. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/icemodel.radiation.updateextcoefs.m`
   - rebuilds `tau_N`, `tau_S`, `k_ext` (and optionally `k_bulk_lookup`) for
     one optical grain-radius index; designed for future grain-size evolution
8. `/Users/mattcooper/MATLAB/projects/icemodel/test/legacy/SPECTRALSOURCETERM_INLINE.m`
   - preserved historical inline reference used only by kernel benchmarks

## How To Reproduce

### Bootstrap

```matlab
clear functions
results = run_spectral_study_bootstrap();
```

### Kernel + whole-model spectral summary

```matlab
clear functions
report = summarize_spectral_perf( ...
    simyear=2016, ...
    include_formal_perf=false, ...
    n_direct_runs=1);
```

### Profile figures

```matlab
clear functions
report = plot_spectral_variant_profiles();
```

Figures are written under:

- `/Users/mattcooper/MATLAB/projects/icemodel/test/benchmarks/figures/spectral/`

### Density-floor audit

```matlab
clear functions
report = summarize_spectral_density_floor();
```

## Interpretation Notes

- `bulk_*` timings and errors are only for the bulk-extinction transform.
- `source_*` timings and errors are for the full spectral source-term path.
- `bulk_*` speedups are relative to `bulk_exact`.
- `source_*` speedups are relative to `source_inlined`.
- The exact `functions` path and the historical `inlined` path should match to
  roundoff. Any real speed difference between them is code-construction cost,
  not a physics difference.
- The `lookup` path changes only the bulk-extinction coefficient calculation.
  The two-stream solve and source-term collapse are still shared.
- Direct whole-model timings use the formal 2-year smoke case contract:
  retained year `simyear` with one leading spinup year.

## Current Decision Logic

- Keep `inlined` only as a preserved legacy reference under `test/legacy/`.
- Keep `functions` as the exact organized reference path.
- Use `lookup` as the accepted production default path.

## Latest Validated Results (`R2024b`)

Kernel summary from:

```matlab
clear functions
report = summarize_spectral_perf( ...
    simyear=2016, ...
    include_formal_perf=false, ...
    n_direct_runs=1);
```

| variant | seconds/call | reference | speedup |
| --- | ---: | --- | ---: |
| `bulk_exact` | `0.0079394` | `bulk_exact` | `1.0000` |
| `bulk_lookup` | `0.00087535` | `bulk_exact` | `9.0699` |
| `source_inlined` | `0.016555` | `source_inlined` | `1.0000` |
| `source_functions` | `0.016450` | `source_inlined` | `1.0064` |
| `source_lookup` | `0.0015460` | `source_inlined` | `10.708` |

Kernel agreement against `inlined`:

| variant | max_rel_err | chi_abs_err |
| --- | ---: | ---: |
| `bulk_lookup` | `2.6863e-05` | `NaN` |
| `source_functions` | `4.4551e-13` | `0` |
| `source_lookup` | `2.6758e-04` | `3.3669e-06` |

Direct whole-model summary for the formal 2-year smoke case:

| variant | median_wall_s | speedup_vs_inlined |
| --- | ---: | ---: |
| `inlined` | `249.55` | `1.0000` |
| `functions` | `129.86` | `1.9217` |
| `lookup` | `40.055` | `6.2301` |

Direct whole-model output agreement:

| variant | pct_dif_mean_Tice_numiter | pct_dif_seb_rmse | max_Tice_numiter | n_not_converged |
| --- | ---: | ---: | ---: | ---: |
| `inlined` | `0` | `0` | `4` | `0` |
| `functions` | `0` | `0` | `4` | `0` |
| `lookup` | `-0.14508` | `-0.10225` | `4` | `0` |

Formal smoke regression acceptance against the accepted `functions` rolling
baseline:

- `spectral_variant="functions"`: passed
- `spectral_variant="lookup"`: passed

Formal smoke perf acceptance under `R2024b`:

- the whole-model perf gate is currently skipped because the accepted rolling
  perf baselines were built under `R2025b`
- the runner now records that incompatibility explicitly instead of reporting a
  false failure
- measured smoke times under `R2024b` were:
  - `functions`: `129.86 s` direct whole-model summary
  - `lookup`: `77.638 s` formal perf `icemodel` case
  - `lookup`: `26.633 s` formal perf `skinmodel` case

Interpretation:

1. `functions` now matches `inlined` to roundoff and is not slower in the
   current validated setup.
2. `lookup` remains the clear production speed path.
3. The solver/control fixes did not introduce any formal `dt_min` or
   `maxsubstep` failures in the `icemodel` full matrix.

## Resolution

`lookup` is accepted as the production default. `opts.lookup_k_bulk` controls
whether the production path uses the lookup-table approximation (`true`, default)
or the exact bulk-extinction transform (`false`).

The `inlined` variant is preserved as a kernel-only reference in `test/legacy/`.
The `functions` path (exact organized helpers) is retained as the analytical
reference and accessible via `opts.lookup_k_bulk = false`.

## Direct-Model Timing Note

The direct whole-model timings above were measured during a long MATLAB session and
appear inflated relative to expected steady-state performance. Expected 1-year
runtimes in a fresh session are approximately:

- `icemodel`, solver 1 or 2: ~12 s
- `icemodel`, solver 3: ~17 s
- `skinmodel`: ~10 s

The formal test suite runs 2-year cases (1 spinup + 1 retained output year), so
formal runtimes are roughly 2x the 1-year values. These timings will be
regenerated in a fresh MATLAB session and updated here before final acceptance.
