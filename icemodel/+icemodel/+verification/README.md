# Verification Namespace

`icemodel.verification` contains the snow-verification workflow used to inspect
staged validation data and compare future model outputs against those targets.
It is separate from the formal regression and performance suites.

## Normal Workflow

Use the top-level functions for ordinary verification runs:

- `icemodel.verification.listcases` lists staged cases from committed manifests.
- `icemodel.verification.loadCaseManifest` resolves one case manifest and its artifact paths.
- `icemodel.verification.comparecase` compares staged targets with a candidate or smoke reference.
- `icemodel.verification.plotcase` plots staged targets, references, or candidate comparisons.

These functions read staged data under `demo/data/eval/snow` by default. They do
not import raw upstream data or mutate the staged dataset tree.

## Candidate Contract

`comparecase` accepts an in-memory `candidate` bundle or a `candidate_file`.
`run_snow_verification_suite` can also receive a `candidate_provider` function
handle. The provider receives one resolved case manifest row from `listcases`
and must return a candidate bundle with the same format as that case's staged
reference:

- ESM-SnowMIP site cases use `format="timeseries"` with a `data` timetable.
- Laugh-Tests process cases use `format="experiment_bundle"` with named
  experiment timetables.

This lets an agent run a synthetic or real snow-model adapter against the
committed targets without changing the staged data.

## Plotting Contract

Plotting uses three independent controls:

- `make_plots` creates comparison figures.
- `save_plots` exports comparison PNGs under the run artifact folder.
- `plot_visible` controls MATLAB figure visibility.

The runner defaults are agent-safe: plots are generated, saved, hidden, and then
closed. For live visual review, run:

```matlab
results = run_snow_verification_suite(plot_visible="on");
```

## Setup Workflow

Setup and refresh tooling lives under `icemodel.verification.setup`:

- `importEsmSnowmip` stages curated ESM-SnowMIP site windows.
- `importLaughTests` stages selected Laugh-Tests synthetic process cases.
- `prepareCaseRoot`, `writeManifest`, `makeFamilyManifest`, `makeCaseManifestEntry`, and `metadataStruct` are setup helpers used by the importers.

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
