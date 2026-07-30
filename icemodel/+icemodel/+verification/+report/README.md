# Snow and firn visual-QA report

This package turns saved verification artifacts into reviewable PDF and HTML
reports.
It does not run the model, reconstruct forcing, or change baselines.

## Getting started

Install Quarto and its PDF engine once:

```sh
quarto install tinytex
```

### Render a PROMICE gap-fill report

After generating fresh `promice_filled` products, render a focused report:

```matlab
addpath("icemodel")
report = icemodel.verification.report.buildGapFillReport( ...
   sites=["kanl", "kanm", "kanu"]);
```

Render the complete on-disk cohort with:

```matlab
report = icemodel.verification.report.buildGapFillReport();
```

The default is `sites="all"` and `render=true`. Use `render=false` to rebuild
and inspect the figures, ledgers, CSVs, and generated QMD without invoking
Quarto:

```matlab
report = icemodel.verification.report.buildGapFillReport( ...
   sites=["kanl", "kanm"], render=false);
```

Default outputs are:

| Output | Location |
|---|---|
| PDF | `data/preview/report/gapfill-report.pdf` |
| HTML | `data/preview/report/gapfill-report.html` |
| Generated QMD | `data/preview/report/gapfill-report.qmd` |
| Method-detail PDF | `data/preview/report/gapfill-detail-report.pdf` |
| Method-detail HTML | `data/preview/report/gapfill-detail-report.html` |
| Method-detail QMD | `data/preview/report/gapfill-detail-report.qmd` |
| Overview figures | `data/preview/figures/gapfill/overview/` |
| Method-detail figures | `data/preview/figures/gapfill/detail/` |
| Scientific-interpretation figures | `data/preview/figures/gapfill/interpretation/` |
| Summary CSVs, figure ledger, and diagnostics | `data/preview/qa/gapfill/` |

The main report contains the scientific Results and station overview
appendix. The companion method-detail report contains the full set of
station/channel/method gap panels and their per-station diagnostic tables;
separating it keeps the main HTML responsive without discarding evidence.

To rerender either generated QMD without rebuilding figures or ledgers:

```sh
quarto render data/preview/report/gapfill-report.qmd
quarto render data/preview/report/gapfill-detail-report.qmd
```

Each command produces both PDF and HTML because both formats are declared in
the QMD. Use `--to pdf` or `--to html` only when intentionally rebuilding one
format.

The builder reads saved artifacts only. If the reconstruction code or native
manifest fingerprints changed, regenerate the products first; see the
[verification README](../README.md#produce-a-fresh-promice-gap-fill-report).

### Other reports

- `buildTestSuiteReport` presents saved regression or performance results.
- `generateFinalSnowPreview` and `generateFinalFirnPreview` reproduce the
  accepted seasonal-snow and firn preview figures.
- `snow-artifact-qa.qmd` is the combined accepted snow/firn artifact-QA index.

## Report-layer ownership

The test suite and verification suite share presentation code, not execution
semantics. `test/` owns unit discovery, numerical/performance runners, managed
baselines, and comparison artifacts. `icemodel.verification` owns scientific
candidate/reference comparisons and artifact QA. This namespace owns the
Quarto-facing presentation for both: the snow/firn scientific report below and
`buildTestSuiteReport`, which turns saved regression/performance results into
plots, compact CSV, generated QMD, and a rendered report without rerunning a
model or changing a baseline.

`run_regression_suite` and `run_perf_suite` render their reports beneath the
same `test/artifacts/<run>/` directory as their MAT files. Pass
`build_report=false` only for a deliberately artifact-only run.

## PROMICE gap-fill report

`buildGapFillReport` renders the reconstruction report from saved native,
filled, plan, audit, and readiness artifacts.

### Evidence and tables

The builder writes an exact figure ledger and
`gapfill_method_diagnostics.csv`. The diagnostics expose every candidate
season × duration admission or denial, maximum validated duration, selection
baseline and baseline RMSE, uncertainty status, one- and two-sigma coverage,
and selection- and evaluation-year metric components.

The Results section is also artifact-derived:

- cohort readiness verdicts, including the relational-diagnostic policy;
- fill volume by provenance family in `gapfill_fill_by_family.csv`;
- cohort admission and held-out skill aggregates;
- staged native stations without a filled product, including
  `windowRecordDisjoint` refusals; and
- a fixed scientific-interpretation catalog in
  `gapfill_interpretation_catalog.csv`, with explicit selectors, status,
  mechanism, policy basis, and absence rows; and
- curated example figures with their selection rule stated in the report.

Every reconstruction option and runtime donor-pool control (`donor_sites`,
`use_ktransect`, and `use_gcnet`) appears in a complete, producer-pinned JSON
block under Methods. The block must be identical across all selected plans.
The appendix contains the complete readiness ledger, not a display-only prefix.

### Figure organization

Station figures follow POLICY D-31 and have matching report sections:

- `overview/` contains one full-period, eight-channel figure per station.
- `detail/` contains at most six windowed figures per station. Candidates are
  interleaved by duration rank within channel × method-family groups; each
  selected figure remains single-method and contains up to six representative
  gaps for that method.
- `interpretation/` contains the bounded cohort-level Results evidence:
  KANL SWD fill geometry, the longest accepted gap, climatology variability,
  boundary behavior, hourly disaggregation, and residual missingness. A
  category with no material case remains in the catalog without a figure.

Both appendix sections have clickable station links, and the figure ledger's
`file` field records each subfolder-relative path. The table of contents and
channel summaries retain sites and channels with no plotted segment, so figure
availability cannot hide complete-native or unresolved data.

Native curves are provenance-masked so legacy staged fills cannot masquerade as
observations. Observed, raw-fallback, and source-backed negative-clamped
shortwave codes remain native context; darkness and reconstructed codes do not.
Each detail panel accents only its named method. Other fills in the context
window are muted grey, derived by `methodFillLayers` from per-sample provenance
and plan-audit spans. `gapfillFigureStyle` owns every report color.

Each ledger row states whether true observed context exists before and after its
exact filled segment. Gap spans and `gap_end` are end-exclusive, so plotted
widths equal audited durations even for one-posting fills. Scientific
interpretation figures use reproducible month-to-year windows and are ledgered
separately from the per-station detail budget.

### Input verification

Before creating output, the builder:

1. resolves producer-manifest paths relative to the selected data root;
2. contains each artifact within its role-specific data or QA root;
3. verifies the recorded size and SHA-256 for every native, filled, plan/audit,
   and readiness artifact; and
4. accepts proxy filenames as the product window only when each timetable and
   their union cover every 15-minute posting continuously between the encoded
   UTC bounds.

The QMD and HTML link each machine-readable CSV at its actual QA location.

### Publication transaction

Figures, CSVs, QMD, and optional HTML are built in one isolated layout. All
outputs must succeed before a rollback-safe transaction replaces the current
publication. For selected sites, the transaction reconciles `overview/`,
`detail/`, and `interpretation/`, removes obsolete selected-site and flat-layout
PNGs, and leaves unrelated subset figures untouched.

`snow-artifact-qa.qmd` renders one local HTML index for the accepted PROMICE,
ESM-SnowMIP, and Laugh/Colbeck preview figures plus canonical RetMIP, IMAU,
SUMup, and research-site firn figures. The Python generator validates both
first-class artifact-QA results and their figure inventories before writing
relative, lazy-loading image links. It does not rerun MATLAB, modify staged
artifacts or canonical figures, copy PNGs, or embed them in the HTML.
Generation fails unless schema-1.0 QA counts reconcile, the audit passed with
zero errors or blockers, and its exact artifact-derived case set matches
current, nonempty post-QA PNGs beneath `data/preview/figures/`. Documented
warnings remain visible in the generated report but do not invalidate a passed
audit. Before consuming readiness or quantitative coverage, the report verifies
every audited artifact's `artifact_sha256` and `artifact_size_bytes` against its
current bytes on disk. PROMICE readiness is then copied from that byte-verified
staging manifest's `forcing_ready`, `forcing_ready_reason`, and
`forcing_complete_windows` fields. Per-case, per-source, per-variable tables
report PROMICE native/MAR/MERRA-2/RACMO legs and
ESM-SnowMIP staged met directly from QA channel records without inferring or
overriding readiness. Case sections link to these tables in the final coverage
appendix so they do not interrupt the figure sequence. They report
`missing_count` as individual **missing
samples**, `missing_run_count` as **contiguous missing runs**, and
`longest_missing_run_samples` as the longest run length; the three quantities
are never conflated. The generated status links both the concise
`artifact_qa.md` and complete `artifact_qa.json`. Figure captions distinguish
native/readiness views from completeness-gated daily diagnostics. Daily albedo
captions distinguish three source contracts. Solar-screened radiometer ratios
require six represented hours and use an SWD-weighted energy ratio. Native
MAR/MERRA-2 modeled surface-state albedo uses an exact-grid arithmetic
finite-sample mean and retains valid polar-night coverage. RACMO/model-ratio
albedo uses positive SWD as the energy weight without the radiometer six-hour
gate; an all-zero-SWD polar-night day remains undefined. Native daily products
without collocated SWD retain finite values, and raw staged
shortwave/model-albedo components remain unchanged.

MAR userdata intentionally stores native radiation components rather than
duplicating runtime-derived balances. For visual comparison, an energy-balance
caption identifies the in-memory-only derivation of missing
`swn = swd * (1 - albedo)`,
`lwn = icemodel.surface.net_longwave_radiation(tsfc, lwd)`, and
`netr = swn + lwn`. The longwave helper uses the same surface emissivity and
Kelvin temperature convention as `icemodel.processmet`; existing native terms
are preserved, display provenance is recorded, and the staged artifact is not
modified.
Optional QA and output overrides must remain beneath the selected preview
root's `qa/` and `report/` directories.

Only `data/preview/figures/` is the canonical figure namespace. Historical
direct family directories were never report inputs and were removed after the
combined report was accepted. Do not recreate them; plots from current QA and
plotting code belong beneath `data/preview/figures/`.

The concise PROMICE model-development handoff is
`data/preview/qa/promice_snow_model_ready_site_years.csv`; the corresponding
coverage ledger and self-check sit beside it. Firn-family quantitative QA,
the exact figure ledger, and figures live beneath `data/preview/qa/firn/`,
`data/preview/qa/firn_figure_index.csv`, and `data/preview/figures/firn/` after
their promotion and canonical convergence gates pass; the same HTML report
links them without copying or embedding the PNGs. The final checker requires
the exact accepted scope: 62 seasonal cases across the three seasonal families
and 56 firn cases across the four firn families. It also requires the firn
figure ledger to equal every PNG beneath `data/preview/figures/firn/` exactly.

The promotion and bounded RACMO precipitation/sublimation repairs are complete
and must not be rerun. The user does not need to run the final-render workflow.
For maintainers who need to reproduce the already-completed canonical figures
and report, run from the repository root:

```sh
python3 test/unit/test_generate_final_snow_preview.py
matlab -batch "addpath('icemodel'); icemodel.verification.report.generateFinalSnowPreview()"
python3 test/unit/test_generate_final_firn_preview.py
matlab -batch "addpath('icemodel'); icemodel.verification.report.generateFinalFirnPreview()"
python3 test/unit/test_build_snow_artifact_qa.py
python3 test/unit/test_check_snow_artifact_qa.py
python3 icemodel/+icemodel/+verification/+report/build_snow_artifact_qa.py
quarto render icemodel/+icemodel/+verification/+report/snow-artifact-qa.qmd \
  --output-dir ../../../../data/preview/report
python3 icemodel/+icemodel/+verification/+report/check_snow_artifact_qa.py
```

The promotion utility remains the durable recovery path for a future promotion.
If one is interrupted and leaves `promotion_transaction.json` or
`transaction_backup` in its evidence directory, do not rerun `--apply` and do
not delete either path manually. Resolve that sealed transaction first:

```sh
python3 icemodel/+icemodel/+verification/+setup/promote_snow_verification_artifacts.py \
  --candidate-root sandbox/verification/icemodel-hfc.2.38/candidate \
  --repo-root . \
  --evidence-dir sandbox/verification/icemodel-hfc.2.38/evidence \
  --recover
```

If recovery reports `state=complete`, the accepted post-state is unchanged;
continue with the independent final dry-run and render steps. If it reports
`state=recovered`, the exact pre-state was restored. Preserve that recovery
evidence and repeat the dry-run/apply/final-dry-run sequence with one fresh
evidence directory, using that same new directory for all three commands. The
promoter verifies the sealed pre/post hashes and refuses stale plans or an
unresolved transaction before changing canonical data.

The promoter inventories every candidate file and requires a converged,
all-identical second plan. Candidate figures, reports, backups, scratch files,
and logs are classified as reproducible or disposable rather than copied into
canonical paths. Remove the local candidate only after the inventory is
complete and the canonical report has been visually accepted.

The Python report-index generator and checker use only the standard library, so
they do not require a virtual environment or additional Python packages. No
Python execution engine is used by this report. Quarto must be on `PATH`, and
TinyTeX must be installed for PDF output.

The ignored outputs are `data/preview/report/figures.generated.md`,
`data/preview/report/snow-artifact-qa.html`, and Quarto support files. Open the
HTML from that directory so its `../figures/promice`,
`../figures/esm_snowmip`, `../figures/laugh_tests`, `../figures/firn`, and
`../qa` links resolve against the current preview tree.
