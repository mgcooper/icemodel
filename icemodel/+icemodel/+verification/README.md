# Verification Namespace

`icemodel.verification` contains the snow-verification workflow used to inspect
staged validation data and compare future model outputs against those targets.
It is separate from the formal regression and performance suites.

## Normal Workflow

Use the top-level functions for ordinary verification runs:

- `icemodel.verification.listcases` lists staged cases from committed manifests.
- `icemodel.verification.loadmanifest` resolves one case manifest and its artifact paths.
- `icemodel.verification.comparecase` compares staged targets with a candidate or smoke reference.
- `icemodel.verification.plotcase` plots staged targets, references, or candidate comparisons.

These functions read staged data under `demo/data/eval/snow` by default. They do
not import raw upstream data or mutate the staged dataset tree.

## Candidate Contract

Candidates are model simulations to be compared with verification targets.
`comparecase` accepts an in-memory `candidate` bundle or a `candidate_file`.
`run_snow_verification_suite` can also receive a `candidate_provider` function
handle. The provider receives one resolved case manifest row from `listcases`
and must return a candidate bundle with the same format as that case's staged
reference:

- ESM-SnowMIP site cases use `format="timeseries"` with a `data` timetable.
- Laugh-Tests process cases use `format="experiment_bundle"` with named
  experiment timetables.

This lets a developer run a synthetic or real snow-model adapter against the
committed targets without changing the staged data.

Use `run_snow_verification_suite(run_icemodel=true)` when the candidate should
come from `icemodel(opts)`. Until production snow physics exists, that route
uses `icemodel.verification.runIcemodelSnowCandidate`, which activates an
explicit verification-only synthetic snow hook inside `icemodel`. Developers
should keep the runner and comparison functions unchanged unless development
is required, replace the synthetic stand-in with real snow-model outputs, and
map those `ice1`/`ice2` outputs through
`icemodel.verification.candidateFromIcemodelOutput`.

Metrics and scatter plots compare only finite paired samples after target and
candidate timetables are synchronized on common timestamps. Missing target data
are not treated as zero and do not contribute to bias, RMSE, correlation, peak
timing, or fitted-line diagnostics.

## Plotting Contract

Plotting uses three independent controls:

- `make_plots` creates comparison figures (`true` / `false`).
- `save_plots` exports comparison PNGs to the run artifact folder (`true` / `false`).
- `plot_visible` controls MATLAB figure visibility (`on` or `off`).

The runner defaults are command-line safe: plots are generated, saved, hidden, and
then closed. For live (interactive) visual review, run:

```matlab
results = run_snow_verification_suite(plot_visible="on");
```

Comparison figures include the time-series overlay plus target-versus-candidate
scatter figures with a 1:1 reference line and fitted linear trend. Scatter
figures are separate from the time-series figure and are only produced for
time-series site cases, not the current Colbeck experiment bundle. The scatter
are generated with `icemodel.plot.scatterplot`. The timeseries plots are generated
with `icemodel.plot.timeseries`.

## Setup Workflow

Setup and refresh tooling lives under `icemodel.verification.setup`:

- `fetchEsmSnowmip` resolves / verifies the local ESM-SnowMIP source cache.
- `fetchLaughTests` resolves / verifies the local Laugh-Tests checkout.
- `importEsmSnowmip` stages curated ESM-SnowMIP site data for select time windows.
- `importLaughTests` stages selected Laugh-Tests synthetic process cases.
- `prepareCaseRoot`, `writeManifest`, `makeFamilyManifest`, `makeCaseManifestEntry`,
   and `metadataStruct` are setup helpers used by the importers.

The setup functions are intentionally separate because they create or overwrite
MAT artifacts and manifests. Use `overwrite=true` only when deliberately
refreshing staged setup data.

### Source-cache layout

Generated / staged smoke artifacts (committed):

```
demo/data/eval/snow/<dataset_family>/<case_id>/{forcing,evaluation,reference}.mat
demo/data/eval/snow/<dataset_family>/manifest.json
```

Local raw source cache (gitignored under `data/verification/**`):

```
data/verification/snow/esm_snowmip/    # cdp + wfj NetCDFs
data/verification/snow/laugh_tests/    # Laugh-Tests checkout
```

Caller pattern:

```matlab
src = icemodel.verification.setup.fetchEsmSnowmip();
icemodel.verification.setup.importEsmSnowmip(src, overwrite=true);

src = icemodel.verification.setup.fetchLaughTests();
icemodel.verification.setup.importLaughTests(src, overwrite=true);
```

When the source cache is incomplete, the fetch helpers print actionable
retrieval instructions (DOI, URL, expected filenames, target paths) and
either error with a stable error id or return the partial path.

## Support Namespaces

- `helpers` contains normal workflow helpers for path discovery, manifest reads,
  artifact loading, candidate resolution, and metric schema definition.
- `namelists` contains canonical selector lists for dataset families, case ids,
  case types, and tiers.
- `validators` contains argument-block validators that consume the namelists.

## Data Contract

Each dataset family has one `manifest.json` under:

`demo/data/eval/snow/<dataset_family>/manifest.json`

Each case folder stores:

- `forcing.mat`
- `evaluation.mat`
- `reference.mat`

Manifests keep case paths relative to the dataset-family folder. Normal workflow
functions resolve those paths to absolute paths at read time.

### Target schema variants

Two evaluation.mat shapes are supported:

1. **Single-bundle** (default for ESM-SnowMIP cdp / wfj):
   ```
   targets.format       = "timeseries" | "experiment_bundle"
   targets.data         (timeseries case)
   targets.experiments  (experiment_bundle case)
   ```

2. **Multi-source** (Colbeck 1976 case): the same evaluation.mat carries
   two reference bundles keyed by source:
   ```
   targets.numerical_summa.experiments.exp{1,2,3}     (frozen SUMMA)
   targets.analytical_clark2017.experiments.exp{1,2,3} (Clark 2017)
   ```
   Generic `comparecase` and `plotcase` callers auto-pick `numerical_summa`
   when the loaded targets struct has no top-level `format` field. The
   case-specific 4-way driver is `icemodel.verification.colbeck.compareSolutions`.

## Variable Mapping Contract

`candidateFromIcemodelOutput(ice1, ice2, opts, manifest)` adapts the icemodel
output (ICE1 / ICE2 timetables / structs) into the candidate bundle consumed by
`comparecase`. Currently supported mappings:

| Verification variable     | Source field    | Derivation                              |
|---------------------------|-----------------|-----------------------------------------|
| `snow_depth_m`            | `ice1.snow_depth` | direct                                  |
| `swe_kg_m2`               | derived         | `snow_depth_m * snow_density_kg_m3`     |
| `surface_temp_C`          | derived         | `Tsfc - Tf` (Tf from physicalConstant)  |
| `bottom_outflow_mps`      | derived         | runoff/outflow proxy from ice2          |
| `snow_liquid_water_storage_m` | derived     | column-integrated f_liq*dz over snow    |

Future snow-model agents who need additional verification variables (cold
content, density profile, f_ice/f_liq snapshots) should extend the adapter and
update this table; do not bury new mappings inside individual cases.

Until production snow physics exists, the suite uses
`verification_synthetic_snow=true` which routes to
`icemodel.verification.syntheticSnowModelRun` and applies hard-coded
perturbations (snow_depth +0.02 m, swe x 1.05, surface_temp +0.25 K,
liquid_water x 1.05) to the staged targets to prove the end-to-end
adapter and comparison path. **The synthetic candidate is NOT a real
model output**; the +5 % storage bias visible in `run_icemodel=true`
metrics is the synthetic perturbation, not a model error. Retirement
of this hook is tracked under `icemodel-tk6.7`.

## Metrics Contract

`comparecase` produces one row per case x experiment x variable pair and
computes the following metrics on aligned finite pairs (`isfinite(target) &
isfinite(candidate)`):

| Metric                    | Variable types        | Description                            |
|---------------------------|-----------------------|----------------------------------------|
| `bias`                    | continuous, sparse    | `mean(candidate - target)`             |
| `rmse`                    | continuous, sparse    | `sqrt(mean((candidate - target).^2))`  |
| `correlation`             | continuous            | Pearson correlation; `NaN` when std=0  |
| `peak_target`             | continuous            | `max(target)` over the comparison window |
| `peak_candidate`          | continuous            | `max(candidate)` over the same window  |
| `peak_error`              | continuous            | `peak_candidate - peak_target`         |
| `peak_time_error_hours`   | continuous            | offset between candidate and target peak times |
| `melt_out_time_error_hours` | snow_depth / swe    | offset between candidate and target return-to-near-zero times |

`status` is `"ok"` when at least one finite pair exists, `"not_applicable"` when
no finite pairs are available (e.g. the candidate omits a variable, or all
observations are missing for the window). `status` is the right column for
filtering before computing summaries.

For the Colbeck multi-source case, `compareSolutions` produces a long-format
table with these same metrics plus `axis_role` (`"formal"` or `"diagnostic"`)
and `target_source` / `candidate_source` columns identifying which pair the
row evaluates. Per-variable RMSE tolerances drive the formal PASS/FAIL summary
(default storage 5 mm, outflow 5e-7 m/s).
