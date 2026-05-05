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

Plotting and artifact writing are opt-in for `run_snow_verification_suite`:

- `make_plots` creates comparison figures (default `false`).
- `save_plots` exports comparison PNGs to the run artifact folder (default `false`).
- `plot_visible` controls MATLAB figure visibility (default `"off"`).
- `write_artifacts` writes `summary.csv` / `summary.mat` / `report.md` and
  per-case `metrics.csv` / `result.mat` (default `false`).

In the default (interactive-development) mode the runner writes nothing to
disk, returns the result struct, and does not even create the run-artifact
directory. The runner promotes flags transitively: `plot_visible="on"` implies
`make_plots=true`; `save_plots=true` implies `make_plots=true` and
`write_artifacts=true`. Pass `write_artifacts=true` to opt in to the
persisted-snapshot workflow regardless of plotting:

```matlab
% Interactive figure, no on-disk artifacts:
results = run_snow_verification_suite(plot_visible="on", save_plots=false);

% Persisted snapshot:
results = run_snow_verification_suite(write_artifacts=true, save_plots=true);
```

Comparison figures include the time-series overlay plus target-versus-candidate
scatter figures with a 1:1 reference line and fitted linear trend. Scatter
figures are separate from the time-series figure and are only produced for
time-series site cases, not the current Colbeck experiment bundle. The scatter
are generated with `icemodel.plot.scatterplot`. The timeseries plots are generated
with `icemodel.plot.timeseries`.

## Time-window policy

Selection is by `case_ids` and optional `startdate` / `enddate`. The
tier abstraction (smoke / full / custom) is gone; selection is
explicit and the same vocabulary applies at both staging and runtime:

- **Staging** (`importEsmSnowmip`): with no explicit window, each
  requested site is staged using its
  `icemodel.verification.helpers.default_smoke_window` (one snow
  water year, second insitu year, `[YYYY-10-01, YYYY+1-09-30]` UTC).
  Pass `startdate` / `enddate` kwargs to stage a different window.
- **Runtime** (`run_snow_verification_suite`): with no explicit
  window and a single ESM-SnowMIP case, the runner narrows to that
  site's `default_smoke_window`. Pass `startdate` / `enddate` to
  override; the staged window remains the upper bound.
- **Default suite call**: `run_snow_verification_suite()` (no args)
  runs Col de Porte (cdp) over its smoke window. CDP is the most
  canonical / widely-cited ESM-SnowMIP snow verification site
  (Menard 2019 ESSD); single-site / single-year keeps interactive
  runs fast.

When the staged window is wider than the runtime window, comparecase
subsets the staged target on the fly via `opts.startdate` /
`opts.enddate` — no re-staging required. See the open architectural
assessment "Custom-date staging vs on-the-fly subsetting" below for
the reasoning behind this contract.

## ESM-SnowMIP sites

10 reference sites are staged from the upstream PANGAEA bundle
(`https://doi.org/10.1594/PANGAEA.897575`):

| code | name                  | location              | insitu years |
|------|-----------------------|-----------------------|--------------|
| cdp  | Col de Porte          | France                | 1994-2014    |
| oas  | Old Aspen             | Saskatchewan, Canada  | 1997-2010    |
| obs  | Old Black Spruce      | Saskatchewan, Canada  | 1997-2010    |
| ojp  | Old Jack Pine         | Saskatchewan, Canada  | 1997-2010    |
| rme  | Reynolds Mountain East| Idaho, USA            | 1988-2008    |
| sap  | Sapporo               | Hokkaido, Japan       | 2005-2015    |
| snb  | Senator Beck          | Colorado, USA         | 2005-2015    |
| sod  | Sodankylä             | Finland               | 2007-2014    |
| swa  | Swamp Angel           | Colorado, USA         | 2005-2015    |
| wfj  | Weissfluhjoch         | Switzerland           | 1996-2016    |

Sites differ in available observation channels: boreal forest sites
(`oas`, `obs`, `ojp`) report `snd_gap_auto` rather than `snd_auto`;
some sites lack `tsl` (no soil-temp obs) or `albs` (no observed
albedo). The builders auto-detect these and the importer derives a
per-site `comparison_variables` list so `comparecase` does not
emit `not_applicable` rows for variables the site never observed.

## Setup Workflow

Setup and refresh tooling lives under `icemodel.verification.setup`:

- `fetchEsmSnowmip` resolves / verifies the local ESM-SnowMIP source cache.
- `fetchLaughTests` resolves / verifies the local Laugh-Tests checkout.
- `buildEsmSnowmipForcing(site, ...)` converts ESM-SnowMIP NetCDF
  files to icemodel-native forcing (timetable). Reusable for
  staging *and* for any future on-the-fly icemodel run.
- `buildEsmSnowmipObservations(site, ...)` converts ESM-SnowMIP obs
  NetCDF files to verification-target observation timetables.
- `importEsmSnowmip` stages all 10 ESM-SnowMIP site cases via the
  builders. With no explicit window, each site uses its
  `default_smoke_window`; pass `startdate` / `enddate` to stage a
  different window.
- `importLaughTests` stages selected Laugh-Tests synthetic process cases.
- `prepareCaseRoot`, `writeManifest`, `makeFamilyManifest`, `makeCaseManifestEntry`,
   and `metadataStruct` are setup helpers used by the importers.

The setup functions are intentionally separate because they create or overwrite
MAT artifacts and manifests. Use `overwrite=true` only when deliberately
refreshing staged setup data.

### Source-cache layout

Generated / staged smoke artifacts (committed):

```sh
# Forcing follows the standard icemodel input layout so configureRun +
# createMetFileNames + loadmet resolve it without verification-only branches.
demo/data/input/met/met_<case_id>_<case_id>_<year>_1hr.mat

# Observation targets and reference bundles stay in the eval tree.
demo/data/eval/snow/<dataset_family>/<case_id>/{evaluation,reference}.mat
demo/data/eval/snow/<dataset_family>/manifest.json
```

Local raw source cache (gitignored under `data/verification/**`):

```sh
data/verification/snow/esm_snowmip/    # 10-site PANGAEA NetCDFs
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

- `helpers` contains normal workflow helpers for path discovery (`evaluationDataRoot`,
  `inputDataRoot`, `snowDataRoot`), manifest reads, artifact loading, candidate
  resolution, metric schema definition, the per-run markdown report writer
  (`writeRunReport`), the per-site default window (`default_smoke_window`), and
  the standard-contract opts builder (`caseSetopts`) used by
  `runIcemodelSnowCandidate`.
- `namelists` contains canonical selector lists for dataset families, case ids,
  case types, the ESM-SnowMIP site-name namelist (`snowmipsite`), the richer
  ESM-SnowMIP catalog (`snowmipcatalog`), and the Laugh-Tests case-id namelist
  (`laughtests`). `caseid` dispatches uniformly across families using these
  per-family namelists.
- `validators` contains argument-block validators that consume the namelists,
  including `mustBeSnowmipSite` for per-site builders.

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

   ```text
   targets.format       = "timeseries" | "experiment_bundle"
   targets.data         (timeseries case)
   targets.experiments  (experiment_bundle case)
   ```

2. **Multi-source** (Colbeck 1976 case): the same evaluation.mat carries
   two reference bundles keyed by source:

   ```matlab
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

| Verification variable         | Source field      | Derivation                              |
|-------------------------------|-------------------|-----------------------------------------|
| `snow_depth_m`                | `ice1.snow_depth` | direct                                  |
| `swe_kg_m2`                   | derived           | `snow_depth_m * snow_density_kg_m3`     |
| `surface_temp_C`              | derived           | `Tsfc - Tf` (Tf from physicalConstant)  |
| `bottom_outflow_mps`          | derived           | runoff/outflow proxy from ice2          |
| `snow_liquid_water_storage_m` | derived           | column-integrated f_liq*dz over snow    |

The forcing side of the verification adapter runs in the opposite
direction: `buildEsmSnowmipForcing(site, ...)` converts ESM-SnowMIP
NetCDF channels (Tair, SWdown, LWdown, Wind, Psurf, Qair, Rainf,
Snowf, plus obs sdepth/albs) into icemodel's native forcing
timetable (tair, swd, lwd, albedo, wspd, rh, psfc, ppt,
snow_depth). The conversion uses
`icemodel.vapor.relative_humidity_from_specific_humidity` for
humidity and `icemodel.physicalConstant('ro_liq')` for the
mass-flux to volumetric-flux conversion of Rainf+Snowf, so all
quantity conversions go through canonical icemodel kernels.

Future snow-model developers who need additional verification variables (cold
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

| Metric                      | Variable types        | Description                              |
|-----------------------------|-----------------------|------------------------------------------|
| `bias`                      | continuous, sparse    | `mean(candidate - target)`               |
| `rmse`                      | continuous, sparse    | `sqrt(mean((candidate - target).^2))`    |
| `correlation`               | continuous            | Pearson correlation; `NaN` when std=0    |
| `peak_target`               | continuous            | `max(target)` over the comparison window |
| `peak_candidate`            | continuous            | `max(candidate)` over the same window    |
| `peak_error`                | continuous            | `peak_candidate - peak_target`           |
| `peak_time_error_hours`     | continuous            | offset between candidate and target peak times |
| `melt_out_time_error_hours` | snow_depth / swe      | offset between candidate and target return-to-near-zero times |

`status` is `"ok"` when at least one finite pair exists, `"not_applicable"` when
no finite pairs are available (e.g. the candidate omits a variable, or all
observations are missing for the window). `status` is the right column for
filtering before computing summaries.

For the Colbeck multi-source case, `compareSolutions` produces a long-format
table with these same metrics plus `axis_role` (`"formal"` or `"diagnostic"`)
and `target_source` / `candidate_source` columns identifying which pair the
row evaluates. Per-variable RMSE tolerances drive the formal PASS/FAIL summary
(default storage 5 mm, outflow 5e-7 m/s).

`comparecase` also reports two snow-season timing diagnostics on
`snow_depth_m` and `swe_kg_m2` series: `snow_onset_time_error_hours`
(first-rise above the variable's threshold) and
`melt_out_time_error_hours` (post-peak first-return below the same
threshold). Peak SWE timing and magnitude are already captured by the
`peak_*` columns above. Together these cover the snow-season diagnostics
listed in the original full-tier plan.

## Open Architectural Assessments

These are written assessments of long-term design questions exposed by
the verification work. They are referenced from beads but live here so
they can evolve alongside the suite.

### startdate / enddate as a first-class IceModel option

`startdate` / `enddate` are now first-class run-window specifiers
across the icemodel runtime. The decision matrix in
`icemodel.configureRun` (canonicalizeWindow):

| `simyears`     | `startdate`/`enddate` | Behaviour                                                                         |
| -------------- | --------------------- | --------------------------------------------------------------------------------- |
| provided       | not provided          | as before; year-granularity contract                                              |
| not provided   | provided              | derive simyears from years touched by the window                                  |
| both, compat.  | both, compat.         | both preserved; window is authoritative for met subset, simyears for met-naming   |
| both, incompat | both, incompat        | error (`icemodel:configureRun:windowOutsideSimyears`)                             |
| neither        | neither               | error (`icemodel:configureRun:noWindow`)                                          |

Notes on the integrations:

- **Validation**: configureRun owns the canonicalization step.
  Half-windows error (`halfWindow`); inverted windows error
  (`invalidWindow`). Compatibility checks fail loudly rather than
  silently producing zero-row met.
- **Restart compatibility**: `use_restart=true` requires `startdate`
  to fall on a calendar-year boundary
  (`icemodel:configureRun:restartWindowNotAligned`). Non-aligned
  restart logic is deferred until there is a concrete need.
- **Output-year filtering**: `opts.output_years` continues to derive
  from `simyears(n_spinup_years+1:end)`. The strict full-year contract
  has been relaxed: spinup years can include partial-year leading
  windows without erroring. Users who want output_years to reflect
  the window can pass it explicitly via `resetopts`.
- **Spinup interaction**: no longer enforces year-boundary alignment.
  `n_spinup_years` simply strips leading retained simyears from
  `output_years`; whether those simyears are partial is the caller's
  decision.
- **Demo runs**: `setopts(..., simyears=[], 'startdate', "2016-06-01",
  'enddate', "2016-08-31")` derives `simyears = 2016` automatically.
  Output paths still use sitename/smbmodel/testname; the window
  itself is not encoded in the path. Callers who want to distinguish
  windowed and full-year runs of the same simyears must use distinct
  testnames.

Test coverage is in `test_run_contracts` under `test_window_*` and
`test_use_restart_with_misaligned_window_errors`.

### Staging window vs on-the-fly subsetting

`importEsmSnowmip(case_ids=..., startdate=..., enddate=...)` stages
per-site forcing/observation MAT artifacts for the requested window.
By default, each site is staged at its `default_smoke_window` — one
snow water year. Wider windows are staged only when explicitly
requested.

The retired tier vocabulary (smoke / full / custom) coupled staging
to fixed window choices; the explicit `case_ids` + optional
`startdate` / `enddate` contract removes that coupling. Agents and
developers that want a narrow window inside a wider staged range
subset on the fly via `opts.startdate` / `opts.enddate` (already
implemented in `icemodel.loadmet`).

Open work (tracked under `icemodel-fx1`):

- `comparecase` / `plotcase` should subset the staged target by
  `opts.startdate` / `opts.enddate` (or a `comparison_window` kwarg)
  before computing metrics.
- `runIcemodelSnowCandidate` / `run_snow_verification_suite` already
  forward `startdate` / `enddate` to `opts`, so the icemodel side is
  already on-the-fly.

### Hourly forcing interpolated to 15 min at runtime

IceModel runs much better at 15-min timestep than 1 hr. ESM-SnowMIP
provides forcing at 1 hr; historically the project has cached 15-min
forcings produced by linear interpolation from hourly. Caching the
interpolated copy duplicates the source and forces a re-cache when
the hourly file changes.

Recommendation: stage the hourly forcing (one source of truth) and
interpolate to the requested `opts.dt` at model run time, inside
`icemodel.loadmet` or a dedicated `icemodel.surface.resampleForcing`
helper. Linear interpolation is appropriate for tair, swd, lwd, wspd,
rh, psfc; precipitation requires conservative resampling (rate
preservation, not point interpolation) so total accumulation is
preserved across the resample. Albedo and snow_depth interpolate
linearly without conservation concerns.

This change touches `icemodel.loadmet` and the existing 15-min cache
generators in the sibling `runoff` project. It belongs to the
forcing-builder framework epic (`icemodel-fkg`) since the same
resample logic should serve both ESM-SnowMIP staging and the
project's broader forcing builder.
