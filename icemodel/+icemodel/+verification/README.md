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

- `importEsmSnowmip` stages curated ESM-SnowMIP site data for select time windows.
- `importLaughTests` stages selected Laugh-Tests synthetic process cases.
- `prepareCaseRoot`, `writeManifest`, `makeFamilyManifest`, `makeCaseManifestEntry`,
   and `metadataStruct` are setup helpers used by the importers.

The setup functions are intentionally separate because they create or overwrite
MAT artifacts and manifests. Use `overwrite=true` only when deliberately
refreshing staged setup data.

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
