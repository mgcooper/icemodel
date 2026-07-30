# Snow and firn visual-QA report

## Report-layer ownership

The test suite and verification suite share presentation code, not execution
semantics. `test/` owns unit discovery, numerical/performance runners, managed
baselines, and comparison artifacts. `icemodel.verification` owns scientific
candidate/reference comparisons and artifact QA. This namespace owns the
Quarto-facing presentation for both: the snow/firn scientific report below and
`buildTestSuiteReport`, which turns saved regression/performance results into
plots, compact CSV, generated QMD, and self-contained HTML without rerunning a
model or changing a baseline.

`run_regression_suite` and `run_perf_suite` render their reports beneath the
same `test/artifacts/<run>/` directory as their MAT files. Pass
`build_report=false` only for a deliberately artifact-only run.

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
they do not require a virtual environment or additional Python packages.
`quarto` only needs to be on `PATH`; no Python execution engine is used by this
report.

The ignored outputs are `data/preview/report/figures.generated.md`,
`data/preview/report/snow-artifact-qa.html`, and Quarto support files. Open the
HTML from that directory so its `../figures/promice`,
`../figures/esm_snowmip`, `../figures/laugh_tests`, `../figures/firn`, and
`../qa` links resolve against the current preview tree.
