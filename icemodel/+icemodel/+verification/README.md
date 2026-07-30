# Verification Namespace

`icemodel.verification` contains the snow-verification workflow used to inspect
staged validation data and compare future model outputs against those targets.
It is separate from the formal regression and performance suites.

## Normal Workflow

Use the top-level functions for ordinary verification runs:

- `icemodel.verification.listcases` lists staged cases from family manifests.
- `icemodel.verification.loadmanifest` resolves one case manifest and its artifact paths.
- `icemodel.verification.auditArtifacts` runs read-only QA/QC over exact
  manifest-referenced files plus the canonical runtime-resolved ESM-SnowMIP met
  artifact required by that family's intentional atomic manifest schema.
- `icemodel.verification.matchObservations` returns reusable matched SUMup
  interval-SMB and dated density/temperature profile rows.
- `icemodel.verification.comparecase` compares staged targets with a candidate or smoke reference.
- `icemodel.verification.plotcase` plots staged targets, references, or candidate comparisons.

These functions read staged data under the top-level
`data/eval/<dataset_family>` tree by default (family-flat taxonomy). They do
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
- SUMup cases use `format="subsurface_profile_bundle"` with accumulated SMB
  interval tables and dated tall density/temperature profile tables.

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

Selection is by `case_ids` and optional `startdate` / `enddate`.
The same vocabulary applies at both staging and runtime:

- **Staging** (`importEsmSnowmip`): with no explicit window, each requested
  site stages its full forcing/observation source record. Pass `startdate` /
  `enddate` kwargs together to stage a shorter shared window. A metadata-only
  `dry_run` uses the site's one-year `default_smoke_window` without reading the
  source cache or writing artifacts.
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
`opts.enddate` — no re-staging required.

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
- `importEsmSnowmip` stages all 10 ESM-SnowMIP site cases via the builders.
  With no explicit window, each site uses its full source record; pass paired
  `startdate` / `enddate` values to stage a shorter window. Metadata-only dry
  runs keep the short per-site `default_smoke_window` preview. Default staged
  met and the manifest runtime cadence are 15 minutes; explicit `dt_out=""`
  keeps both the met artifact and its manifest cadence hourly.
- `importLaughTests` stages selected Laugh-Tests synthetic process cases.
- `prepareCaseRoot`, `writeManifest`, `makeFamilyManifest`, `makeCaseManifestEntry`,
   and `metadataStruct` are setup helpers used by the importers.

### Canonical family API

Every public family importer uses `case_ids` for its case selector. The
non-atomic firn importers also use the same runtime-source surface:
`forcing_sources`, `build_observations`, and `build_forcing`. Their
`build_forcing` default is `false`; a caller must opt in before any runtime
met/userdata artifacts are written. `overwrite=false`,
`overwrite_family=false`, paired `startdate`/`enddate`, `skip_missing`,
`dry_run`, and the staging-root options have the same meaning across families.
Strict cache validation belongs to each `fetch*` function; importers select
fail-fast versus per-case skips with `skip_missing`.

The remaining public differences are source contracts, not aliases:

| Family | Case/source-specific inputs | Why it differs |
|---|---|---|
| ESM-SnowMIP | `case_ids`, `dt_out` | Forcing and observations share one source window and stage atomically; there is no independently attachable RCM leg. |
| Laugh Tests | `case_ids` | Each case is an atomic evaluation/reference experiment bundle with no runtime forcing artifact. |
| PROMICE | `forcing_sources` may include native `"promice"` plus RCM ids | PROMICE owns both observations and an optional native AWS runtime leg. |
| RetMIP | protocol/native source directories | Protocol, Vandecrux, Samimi, PROMICE, and IMAU products require explicit family parsing before optional RCM attachment. |
| IMAU | hourly AWS plus daily-QA cache | The daily product is provenance/QA, not another staged case inventory. |
| SUMup | `points`, `anchors`, `years`, spatial radii | Cases are selected from geospatial observations rather than a fixed site catalog. |
| research_site | fixed `family="research_site"`, `observation_source="sumup"` | The current Humphrey adapter derives observations from SUMup and owns no native station forcing. |

Fetchers use `products` for independently cached products, `stations` for a
physical-station subset, and SUMup's `variables`/`region` for its grouped
geospatial release. Builders likewise keep source-faithful selectors:
`station`/`site` for native AWS readers, `location` for gridded RCM extraction,
and `point` for SUMup observations. Those names do not represent the common
manifest case selector and must not be replaced with compatibility aliases.

The setup functions are intentionally separate because they create MAT
artifacts and manifests. Ordinary `overwrite=false` calls are additive: they
reuse current requested artifacts, add missing artifacts/sites/sources, and
preserve unrelated cases and source legs. Use `overwrite=true` only when
deliberately replacing requested staged setup data.
An enclosing staged window satisfies a narrower repeat. A wider same-case
request makes a fixed-name observation bundle stale, rebuilds it with a visible
replacement warning, and writes widened window-named met/userdata artifacts.
An equal/covered fixed-name `observations.mat` is reusable only when its saved
target metadata has no concrete scalar identity conflict with the requested
producer. Equal facts preserve bytes, concrete producer/product/schema/method
changes select a visible rewrite, and missing legacy facts remain
unknown-compatible. This check reads only the already-staged MAT target during
non-dry observation staging; dry runs do not gain artifact or raw-source reads.

### v1.1 release-data tooling

`test/assets/icemodel-v1.1-data-manifest.json` is the single source of truth for
the required `formal-core` and `verification-showcase` capabilities and the
optional `forcing-integration` capability. `fixtureFileList` selects manifest
rows; `packFixtures` writes one archive per selected capability; and
`fetchFixtures` verifies or transactionally installs selected archives.

Calling `fetchFixtures("v1.1")` explicitly provisions the mandatory capabilities
and downloads missing release archives. Callers that only verify canonical data
must pass `download=false`; local manifest/archive overrides preserve offline
operation. Missing mandatory data in network-free verification reports the exact
MATLAB provisioning command. See `README_FIXTURES.md` for the trust-boundary
contract.

### Source-cache layout

Generated / staged artifacts (top-level data tree):

```sh
# Forcing follows the standard icemodel input layout so configureRun +
# createMetFileNames + loadmet resolve it without verification-only branches.
# writemet/writeuserdata stage into a per-source subfolder (NOT bundled with the
# eval target - the eval is forcing-AGNOSTIC: any forcing usable at runtime).
# Multi-year staged windows produce a single window-stamped file per
# source/product token:
data/input/met/<source>/met_<site>_<source>_<YYYYMMDD>_<YYYYMMDD>_15m.mat
data/input/userdata/<source>/<site>_<source>_<YYYYMMDD>_<YYYYMMDD>.mat

# RCM selectors use the shared public forcing API
# (forcing_sources=["mar","merra","racmo"]), while staged
# RCM products use versioned source/product tokens in folders, filenames, and
# source lists so future versions can coexist:
data/input/met/mar3.11/met_<site>_mar3.11_<YYYYMMDD>_<YYYYMMDD>_15m.mat
data/input/met/merra2/met_<site>_merra2_<YYYYMMDD>_<YYYYMMDD>_15m.mat
data/input/userdata/racmo2.3p3/<site>_racmo2.3p3_<YYYYMMDD>_<YYYYMMDD>.mat

# Eval targets stay in the eval tree. The eval target is a data-only
# observations.mat bundle per case (forcing-agnostic; ESM-SnowMIP, SUMup, and
# freshly staged PROMICE). The analytical families add a reference.mat (the
# frozen SUMMA / analytical solution). The taxonomy is dataset-family-flat:
# families live directly under eval/ with no intermediate snow/ or firn/ level
# (the physical regime is recorded per case in the surface_zone manifest field).
data/eval/<dataset_family>/<case_id>/observations.mat   # eval target
data/eval/<dataset_family>/<case_id>/reference.mat      # analytical only
data/eval/<dataset_family>/manifest.json
```

Local raw source cache (gitignored under `data/verification/**`, same
family-flat layout):

```sh
data/verification/esm_snowmip/    # 10-site PANGAEA NetCDFs
data/verification/laugh_tests/    # Laugh-Tests checkout
data/verification/promice/        # PROMICE/GEUS AWS + co-located RCM sources
data/verification/retmip/         # RetMIP protocol, outputs, scripts, Samimi
data/verification/gcnet/          # Vandecrux GC-Net surface/firn products
data/verification/imau/           # IMAU hourly S21/S22/S23 + daily QA
data/verification/sumup/          # SUMup firn profiles (access-gated)
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
GC-Net/Vandecrux station files may vary only in case and separator spelling:
each registry rule must resolve to exactly one normalized basename. Partial
filenames and duplicate package copies keep the product incomplete rather than
letting later inventory code select an arbitrary file. XML remains optional
provenance and is reported separately from native-data completeness.
Tolerant `gcnetInventory` calls consume the fetch status row's `resolved_files`
station/suffix/filename bindings directly: they may expose valid classes from an
incomplete product, but never reclassify a full path, revive a partial name, or
choose one member of an ambiguous class. They also reuse the status row's optional
XML path, so data and provenance have one cache scan and identity authority.

### Snow setup and restage

Snow verification artifacts are separate from the firn-development staging
workflow. They share the same top-level `data/` root but do not use the RCM
forcing/Data workflow.

```matlab
addpath("icemodel")
icemodel.dependencies()

data_root = string(icemodel.internal.fullpath("data"));

esm_dir = fullfile(data_root, "verification", "esm_snowmip");
laugh_dir = fullfile(data_root, "verification", "laugh_tests");

icemodel.verification.setup.importEsmSnowmip( ...
   esm_dir, ...
   output_root=data_root, ...
   overwrite=true, ...
   overwrite_family=true);

icemodel.verification.setup.importLaughTests( ...
   laugh_dir, ...
   output_root=data_root, ...
   overwrite=true, ...
   overwrite_family=true);
```

### Firn development preview and full staging

Use the short preview when an importer, source cache, or plotting change needs a
fast production-path preflight before a full restage. It stages representative
one-year-or-shorter `build_forcing=true` cases through the canonical importers
into `data/preview/firn_staging`, writes figures beside that preview tree, and
leaves the final `data/eval` and `data/input` trees untouched. It is not required
after a canonical tree has passed QA, and it neither repairs nor promotes that
tree; rerun one family when isolating a future regression.

```matlab
addpath("icemodel")
icemodel.dependencies()

data_root = string(icemodel.internal.fullpath("data"));
preview_root = fullfile(data_root, "preview", "firn_staging");

promice_dir = fullfile(data_root, "verification", "promice");
retmip_dir = fullfile(data_root, "verification", "retmip");
gcnet_dir = fullfile(data_root, "verification", "gcnet");
samimi_dir = fullfile(retmip_dir, "samimi");
imau_dir = fullfile(data_root, "verification", "imau");
sumup_dir = fullfile(data_root, "verification", "sumup");

mar_dir = "/path/to/MARv3.11/RUH2";
merra_dir = "/path/to/MERRA2/1hrly/ncfiles";
racmo_dir = "/path/to/RACMO2.3p3/subsurface";
modis_dir = "/path/to/GEUS/MODIS/albedo/gris";

rcm = { ...
   "mar_dir", mar_dir, ...
   "merra_dir", merra_dir, ...
   "racmo_dir", racmo_dir, ...
   "modis_dir", modis_dir};
```

Every lasting artifact-contract correction must be implemented in the canonical
builder or importer, but accepted artifacts do not always require a full raw
restage. The durable partial-repair tools have deliberately separate scopes:

- refreshManifestSourceLists(manifest_file) rewrites only derived
  forcing_sources and eval_sources from authoritative colocation state. It does
  not read MAT files or rebuild observations.
- repairRcmArtifactMetadata(...) inventories every current-product RCM cache file
  when unscoped and reports files absent from current manifests as unmapped.
  dataset_family restricts the run to exact family-manifest references, while
  source_id filters the product inventory. The function defaults to a dry run,
  synchronizes canonical metadata, validates preservation boundaries, replaces
  files atomically, and reports pass-two identity. Its optional repair_function
  callback is the extension seam for a
  future bounded field/property migration after the same change is canonical in
  the production builder. The caller must declare every variable and UserData
  field or CustomProperty the callback may change; the coordinator rejects
  time-axis, schema-order, and undeclared payload/property changes. New columns
  must be appended after retained columns. Existing source or sampling identity
  that conflicts with the current manifest is never repaired by restamping.
- repairMetTimeSupport(files, ...) remains the explicit recovery path for
  legacy 15-minute met files written with the old
  linear_adjacent_finite_only policy. It proves native rows from saved
  provenance, delegates the corrected hold to
  icemodel.forcing.helpers.resampleMetTimestep, and preserves unrelated MAT
  variables. It is not currently an arbitrary cadence converter.
- stageRcmForcing manifest mode remains the additive path for cache
  reattachment, deriving missing MAR/MERRA met from compatible Data, and staging
  true misses.

Use repair only when the operation is manifest-selected or explicitly
file-selected, dry-runnable, idempotent, atomic, and able to prove preservation
of unrelated data. Observation changes, point/source-grid selection changes,
unsupported cadence changes, and transformations without sufficient provenance
still require canonical staging. Temporary incident wrappers may call the
durable coordinator, but are removed after the migration; their scientific
correction remains in the canonical builder.

For a metadata-only RCM refresh:

~~~matlab
dry = icemodel.verification.setup.repairRcmArtifactMetadata( ...
   dataset_family="promice", source_id="mar3.11");
written = icemodel.verification.setup.repairRcmArtifactMetadata( ...
   dataset_family="promice", source_id="mar3.11", dry_run=false);
second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
   dataset_family="promice", source_id="mar3.11");
assert(second.summary.unchanged == second.summary.total)
~~~

A payload migration supplies
[artifact, actions] = repair_function(artifact, context) plus
allowed_variable_changes, allowed_metadata_changes,
allowed_custom_property_changes, and/or allowed_table_property_changes. The
callback contains only the temporary transform; manifest discovery, mutation
bounds, hashing, atomic replacement, and convergence reporting stay reusable.

For the proven legacy met-time repair:

~~~matlab
files = [
   "/absolute/path/to/met_site_source_20000101_20001231_15m.mat"
   "/absolute/path/to/met_other_source_20000101_20001231_15m.mat"
];
assert(all(isfile(files)), "replace files with exact existing artifact paths")
dry = icemodel.verification.setup.repairMetTimeSupport(files);
written = icemodel.verification.setup.repairMetTimeSupport( ...
   files, dry_run=false);
second = icemodel.verification.setup.repairMetTimeSupport(files);
assert(second.summary.unchanged_count == second.summary.file_count)
~~~

The current userdata writer supports explicit native cadence or canonical
hourly output, not 15-minute userdata. A future userdata-cadence change must
first be defined in writeuserdata; only then should the partial-repair
orchestrator be extended to rewrite artifact filenames and manifest references
under that canonical policy.

Run all preview families and generate one-year figures. The preview artifacts
are already clipped to the lesser of one year or the case period, but the plot
call still uses the same `startdate` / `enddate` options used for full staged
artifacts:

```matlab
preview = icemodel.verification.setup.previewFirnStaging("all", ...
   "output_root", preview_root, ...
   "promice_dir", promice_dir, ...
   "retmip_dir", retmip_dir, ...
   "gcnet_dir", gcnet_dir, ...
   "samimi_dir", samimi_dir, ...
   "imau_dir", imau_dir, ...
   "sumup_dir", sumup_dir, ...
   rcm{:});

preview_summary = icemodel.verification.plotVerificationArtifacts( ...
   "output_root", preview_root, ...
   "dataset_family", "all");
```

Run one family at a time when isolating source problems:

```matlab
icemodel.verification.setup.previewFirnStaging("promice", ...
   "output_root", preview_root, "promice_dir", promice_dir, rcm{:});

icemodel.verification.setup.previewFirnStaging("retmip", ...
   "output_root", preview_root, "retmip_dir", retmip_dir, ...
   "gcnet_dir", gcnet_dir, "samimi_dir", samimi_dir, ...
   "imau_dir", imau_dir, "promice_dir", promice_dir, rcm{:});

icemodel.verification.setup.previewFirnStaging("imau", ...
   "output_root", preview_root, "imau_dir", imau_dir, rcm{:});

icemodel.verification.setup.previewFirnStaging("research_site", ...
   "output_root", preview_root, "sumup_dir", sumup_dir, rcm{:});

icemodel.verification.setup.previewFirnStaging("sumup", ...
   "output_root", preview_root, "sumup_dir", sumup_dir, rcm{:});

plot_start = datetime(2016, 5, 1, 'TimeZone', 'UTC');
plot_end = datetime(2017, 4, 30, 23, 59, 59, 'TimeZone', 'UTC');
icemodel.verification.plotVerificationArtifacts( ...
   "output_root", preview_root, ...
   "dataset_family", "retmip", ...
   "startdate", plot_start, ...
   "enddate", plot_end);
```

After visually checking `preview_root/figures`, rebuild the full final artifact
tree under top-level `data/`:

```matlab
common = { ...
   "output_root", data_root, ...
   "build_forcing", true, ...
   "overwrite", true, ...
   "overwrite_family", true, ...
   rcm{:}};

promice = icemodel.verification.setup.importPromiceSites( ...
   promice_dir, common{:});

imau = icemodel.verification.setup.importImau( ...
   imau_dir, common{:}, "dry_run", false);

retmip = icemodel.verification.setup.importRetmip( ...
   retmip_dir, common{:}, "dry_run", false, ...
   "promice_dir", promice_dir, ...
   "gcnet_dir", gcnet_dir, ...
   "samimi_dir", samimi_dir, ...
   "imau_dir", imau_dir);

research_site = icemodel.verification.setup.importResearchSites( ...
   sumup_dir, common{:}, ...
   "family", "research_site", ...
   "observation_source", "sumup", ...
   "case_ids", "humphrey");

anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
   "output_root", data_root);

sumup = icemodel.verification.setup.importSumup( ...
   sumup_dir, common{:}, ...
   "points", zeros(0, 2), ...
   "anchors", anchors);

% First write a full-period figure set so every staged artifact is represented
% in the production QA output.
full_summary = icemodel.verification.plotVerificationArtifacts( ...
   "output_root", data_root, ...
   "dataset_family", "all", ...
   "figure_root", fullfile(data_root, "figures", "firn_staging_full"));

% For human review, also inspect readable one-year slices. The years below are
% an example subset; before signoff, repeat/extend the loop to cover every year
% present in the staged manifests you intend to use.
plot_years = 2012:2018;
for yy = plot_years
   year_start = datetime(yy, 1, 1, 'TimeZone', 'UTC');
   year_end = datetime(yy, 12, 31, 23, 59, 59, 'TimeZone', 'UTC');
   final_summary = icemodel.verification.plotVerificationArtifacts( ...
      "output_root", data_root, ...
      "dataset_family", "all", ...
      "figure_root", fullfile(data_root, "figures", ...
      "firn_staging_full", string(yy)), ...
      "startdate", year_start, ...
      "enddate", year_end);
end
```

The top-level firn importers default to `data/verification/<family>`,
`data/eval`, and `data/input` using `icemodel.internal.fullpath("data")`. They do
not require calling `icemodel.config()` first, and their defaults do not depend
on the active `ICEMODEL_DATA_PATH`. Pass `output_root`,
`evaluation_data_root`, `input_data_root`, or a nonblank
`icemodel_config_casename` only when deliberately targeting another tree.
`plotVerificationArtifacts` follows the same root convention. With no
`startdate` / `enddate`, it plots the full available period of each staged
artifact. Pass `overwrite=true` to clear prior PNGs only from each selected
case's figure directory before regenerating renamed or regrouped outputs.
All shared time-series panels preserve explicit missing values and break lines
across safely inferred omitted-time gaps. Colors are keyed by source rather than
series order; MAR, RACMO, and MERRA-2 therefore retain their canonical runoff
palette in every panel and source subset.
All staged families use one source-aware daily albedo path. Collocated
radiometer albedo is weighted by positive downwelling shortwave, which is
equivalent to the ratio of daily reflected and incoming shortwave energy and
therefore prevents nighttime or low-sun ratios from dominating the visual QA.
Radiometer source adapters first require solar elevation greater than 20
degrees, and dense radiometer days require at least six hours of valid weighted
support. Native RCM albedo states (MAR ALH and MERRA-2 SNICEALB) use an
exact-grid finite-sample arithmetic daily mean, so valid modeled polar-night
coverage is preserved without the observation support gate; a source day with
no native state remains missing. RACMO 2.3p3 does not provide an albedo state in
the staged source files, so its albedo is the ratio derived from net and
downwelling shortwave. That modeled ratio uses the daily shortwave-energy ratio
whenever at least one positive-shortwave native sample exists, without the
radiometer six-hour gate, and remains undefined through zero-shortwave polar
night. Already-daily albedo products without collocated shortwave, including
GEUS MODIS, retain their finite-sample values. Daily wind direction is the
direction of the mean wind vector rather than an arithmetic degree mean;
collocated wind speed supplies the vector weight. Other uniformly sampled
scalar state and flux channels retain arithmetic daily means because their
regular support gives every sample equal duration. Precipitation and
mass-balance rate channels retain cadence-aware daily totals. Sparse interval
and profile observations remain on their native support rather than being
retimed. These rules are centralized in the shared plotting path, so PROMICE,
ESM-SnowMIP, Laugh Tests, RetMIP, IMAU, SUMup, and research-site figures inherit
the same reduction without family-local copies.
Paired met-forcing scalars use a solid line on the left axis and a dotted line
on the right; the downwelling shortwave/longwave pair uses the same
solid/dotted order. Incremental surface-mass-balance panels display complete
daily totals in `mmWE`. `smb` interval observations remain in `mWE` and are
drawn as disjoint horizontal spans, so unrelated campaigns are never joined.
RetMIP `snowf_subl` is labelled as net snow accumulation: positive values are
snowfall/deposition and negative values are sublimation.

The surface-state figure keeps `surface_height` and `snow_depth` distinct.
`surface_height` is signed relative surface-elevation change (positive upward),
not a snow-depth alias; `snow_depth` remains nonnegative physical snow height.
The plotted `tice10m` is the canonical staged 10 m ice/firn temperature channel,
not a value synthesized by the renderer. Profile legends use family-facing
observation labels and split a small, explicit `name` identity set into separate
colors; profile collections without such an identity remain one honest series.
`plotFirnArtifacts` forwards the same option as a compatibility wrapper, but new
workflow examples should use the family-neutral `plotVerificationArtifacts`.

### RCM met/userdata workflow

RCM staging is owned by `stageRcmForcing` and is called by the firn importers
when `build_forcing=true`. The public selector is `forcing_sources`; values use
the short source-family names:

```matlab
forcing_sources = ["mar", "merra", "racmo"];
```

The staged artifacts use product ids:

- MAR 3.11 writes full Data to `data/input/userdata/mar3.11/` and icemodel met
  to `data/input/met/mar3.11/`.
- MERRA-2 writes full Data to `data/input/userdata/merra2/` and icemodel met to
  `data/input/met/merra2/`.
- RACMO 2.3p3 writes Data only to `data/input/userdata/racmo2.3p3/`; the
  available RACMO files do not include the near-surface state variables needed
  for a valid icemodel met file.

All model-met write paths default to `dt_out="15m"`, including native station
met, ESM-SnowMIP met, cached-met recovery, and MAR/MERRA RCM staging. Pass
`dt_out=""` through the public builder/importer/preview to retain native
model-met cadence. The known RetMIP Samimi Dye-2 override therefore writes and
resolves a `_30m.mat` model-met artifact; the default remains `_15m.mat`.
Userdata write paths default to `dt_out="1h"`: native hourly
PROMICE and RCM products remain unchanged, finer products are aggregated to
clock-hour means, coarser products are interpolated to hourly support, and wind
direction follows circular rather than scalar arithmetic. Pass `dt_out=""`
directly to `writeuserdata` only for an explicitly native-cadence artifact.

Manifests keep stable colocation keys (`colocation.mar`, `colocation.merra`,
`colocation.racmo`) so call sites can still request model families. Each staged
RCM leg also carries `source` (the short family selector) and `source_id` (the
product id). Manifest `forcing_sources` and `eval_sources` use product ids such
as `mar3.11`, `merra2`, and `racmo2.3p3`, because those are the runnable or
comparable staged products.

Partial rebuilds are source-local. PROMICE imports use
`forcing_sources=["promice","mar","merra","racmo"]` to select which runtime
artifacts are written; sibling importers and shared RCM helpers use
`forcing_sources=["mar","merra","racmo"]`. If a case already has MAR, MERRA,
and RACMO legs, a later
source-limited restage refreshes only the requested source and preserves every
omitted case/source leg verbatim, even if the patched point or period changes.
When that requested source needs a wider window, an overlapping but shorter
cache is not sufficient: staging rebuilds/writes the full requested window. If
raw rebuilding is unavailable, the clipped fallback remains only with a warning
and manifest note. Adding an unrelated site does not auto-widen earlier sites.
Directory discovery attaches cached RCM files only when their saved point and
sampling-method metadata match the request. A prior manifest may preserve its
explicitly referenced unstamped legacy files after a failed refresh when every
reference exists, its bounded window overlaps the case, and its declared method
matches. Missing referenced files and concrete saved point/method conflicts are
validated requested-leg replacements: the staging layer emits a transient
replacement signal that the additive manifest merge consumes, so invalid prior
references cannot be restored and the signal is never written to the manifest.
Compatible requested failures and every omitted leg remain additive. Staging and
fallback resolve the same default met/userdata roots, so omitting output paths
does not change this behavior. Manifest staging keeps two window concepts
separate: the full case period is the discovery/ranking window for compatible
cached files, while the source-clipped resolved leg is the required cache
coverage, rebuild/write window, and staged manifest `leg.window`.

PROMICE evaluation and runtime availability are separate manifest facts. The
native `colocation.promice.staged` flag is true only when that leg references at
least one met or userdata artifact; `eval_staged` records availability of the
case-level `observations.mat` bundle. Thus observation-only and fresh RCM-only
calls retain `promice_obs` in `eval_sources` without advertising a nonexistent
PROMICE forcing leg. Source-list refreshes prefer `eval_staged` when present and
fall back to `staged` for legacy manifests. An additive source-selective refresh
preserves an existing valid PROMICE leg verbatim.
`overwrite_family=true` is the explicit destructive boundary: the family is
rebuilt from the current request, prior omitted state is not carried forward,
and replacement emits a warning when it actually removes prior cases, sources,
family extension fields, period coverage, evaluation/reference artifacts, staged
leg state, or met/userdata references. Pure additions remain warning-free.
Requested RCM sources are stable-deduplicated and checkpointed one source at a
time in both dataset and standalone manifest staging, so a later unexpected
source failure leaves every earlier completed source linked in the manifest.

For RetMIP, IMAU, SUMup, and research-site cases,
`build_observations=false` is a guarded forcing-only fast path. Every requested
case must already exist in the target family manifest; its observation/native
entry is retained while the current call attaches `forcing_sources`. This path
does not reread the family observation cache. SUMup additionally requires
explicit `case_ids`, because point discovery is part of observation staging.
SUMup reuse and fresh point staging converge on the same final manifest path.
With `build_forcing=false`, importers do not create input directories, discover
or derive cached met/userdata, write input artifacts, or attach new forcing
legs; the call is an observation/manifest patch only.
For an IMAU observation refresh, that patch merges only compatible prior native
met/data references, readiness, and artifact provenance into the freshly built
IMAU leg. Current source paths, observation window, evaluation reference, daily
QA, and observation metadata win. RetMIP uses the analogous selective split:
fresh protocol-owned fields overlay retained native runtime state. Compatibility
is producer-specific: IMAU compares source family, station, and any finite saved
point; RetMIP compares native family, source id, relationship, and the shared
scalar identity fields (including the documented method alias). An identity
conflict leaves old bytes untouched but removes their manifest references and
marks the fresh native leg as requiring `build_forcing=true`; stale provenance
cannot be relabelled as the new source.
ESM-SnowMIP stays atomic because its source forcing and observations share one
case/window conversion. Laugh Tests stay atomic because they stage only an
evaluation/reference bundle and have no independently attachable runtime
forcing.
Optional `startdate`/`enddate` inputs are one paired UTC contract: both or neither,
chronological, and validated before source/cache reads or staging-root creation.
RCM `method` accepts only `nearest` or `natural` at the same non-mutating boundary.
RACMO point extraction uses only cells with the native fractional IceMask rounded
to ice (fraction at least 0.5; `Promicemask > 0` fallback) and accepts the nearest
valid cell only within one native grid diagonal. A site with no local valid cell
is explicitly unavailable; staging does not snap it to a distant inland cell or
preserve an older off-mask artifact. Conservative polygon remapping uses the
same mask.

### Top-level staging roots and demo data

Normal verification staging targets the top-level `data/` tree. The tracked
`demo/data` tree serves only the three public examples and is not a
verification discovery, import, forcing-builder, comparison, or plotting root.
Formal KAN fixtures belong to the provisioned `test/data` capability:

```sh
test/data/input/met/met_kanl_kanl_2015_15m.mat
test/data/input/met/met_kanm_kanm_2015_15m.mat
test/data/input/met/met_kanm_kanm_2016_15m.mat
```

Standard unit tests must use temporary/synthetic fixtures or provisioned
`test/data` explicitly. Demo scripts target `demo/data` only through
`config("demo")`.

### PROMICE co-located firn staging (`importPromiceSites`)

`importPromiceSites` stages PROMICE + MAR + MERRA-2 met + RACMO Data anchored on
each PROMICE AWS site. The eval target IS bundled as a data-only
`observations.mat` (the PROMICE-obs timeseries; same contract as
ESM-SnowMIP/SUMup), referenced from the case `evaluation_file`. The manifest is
**forcing-AGNOSTIC**: the FORCING is NOT bundled with the eval target and is NOT
stipulated - any forcing file may be used at the site at runtime without
rewriting `observations.mat`. Forcing/Data live in individual
`<source>/met_<site>_<source>_<window>.mat` / window-stamped userdata files via
the standard naming convention (`writemet`/`writeuserdata` stage into
per-source/product subfolders such as `met/promice/`, `met/mar3.11/`, and
`userdata/racmo2.3p3/`) and are recorded in the manifest `colocation` record by
source id. `forcing_sources` declares the runnable met source ids used by runtime
helpers; it must match the staged product ids when forcing is available. The
per-leg windows are
**decoupled**:
PROMICE met/eval defaults to **all available years for each station** (read live
from the L3 record, including partial years) and is **never gated by RCM
coverage**; the MAR/MERRA met legs use the PROMICE window intersected with each
source's on-disk availability, while the RACMO Data leg uses its own coverage
(FGRN11 surface ~2012-2015, subsurface ~2012-2018) independently. A leg with no
on-disk overlap is skipped-with-reason (recorded at `colocation.<model>`), and a
requested-vs-actual coverage table is printed at the start of every run. Each
leg's actual staged window is recorded at `colocation.<model>.window`.

RACMO is staged as eval/reference Data, NOT a met source: the available RACMO
2.3p3 source files carry radiation (swd, lwd, derived albedo), turbulent fluxes
(shf, lhf), precip and SMB components, but LACK the near-surface meteorological
STATE variables (tair, wspd, rh, psfc). (RACMO in general carries these; only
the available source files omit them - obtain the full set from the RACMO
developers, or borrow tair/wspd/rh/psfc from MAR/MERRA/PROMICE at the point, for
a RACMO-forced run.) The same applies to the RACMO leg of the SUMup co-located
forcing (`buildSumupForcing` / `importSumup`).

PROMICE `tice10m` is the primary subsurface-temperature comparison channel. It
is the GEUS standardized temperature 10 m below the evolving surface, not a
fixed-depth thermistor from installation, so model comparisons sample the
model's surface-following temperature column at 10 m. See
`+forcing/README.md` for the full thermistor and `tice10m` protocol.

The `output_root` kwarg is the explicit staging-root switch: eval manifests go
to `<output_root>/eval`, forcing/Data to `<output_root>/input`.
`output_root=<repo>/data` (or the unset default, which resolves to the repo-root
data tree) is the normal verification-data target. Verification staging must
not target `demo/data`; that tracked tree is limited to the public demo inputs.

**Incremental staging (MERGE by default).** Staging one site **adds or updates
only that site's case entry** in the family `manifest.json` and **preserves
every other site's committed case + staged files byte for byte** (shared helper
`icemodel.verification.setup.writeFamilyManifestMerge`, used by the dataset
family importers). The existing manifest at the target
root is read, the requested sites' cases replace/append (matched by `case_id`),
and untouched cases re-encode identically (raw decode, no field reordering);
hand-added family fields like `schema: "metadata_only"` survive. Re-staging the
same site updates exactly its entry (idempotent), and a stale `skipped[]` entry
for a now-staged site clears while other sites' skips are preserved. So adding
DY2/EGP into the family root that already holds the KAN fixtures never churns or
drops them. A shorter same-identity refresh preserves the enclosing case/leg
window and its prior artifact references; an enclosing/equal rebuild replaces
those references, while a partial extension retains both file sets supporting
the overlapping union. A disjoint refresh is warned and leaves the prior scalar
window/artifact state intact rather than inventing continuous coverage across an
unproven gap. Blank all-available coverage is never narrowed by an ordinary
bounded call, and source/product, schema, sampling-method, or point conflicts
replace rather than union unlike artifacts. Production `method` and manifest
`sample_method` are one documented identity alias. Native producer
`source_family`, `station`, `doi`, and `bundle_doi` fields compare by those exact
spellings when both records provide them; missing legacy values remain unknown,
and no additional cross-field aliases are implied. Known cadence is compared per
artifact class, so incompatible met or Data refs replace rather than union while
a valid 15-minute-met/hourly-Data leg remains compatible. Durable filename cadence
is used only when every ref in that class is unambiguous; unknown legacy names
remain compatible. Source lists are derived again from the final colocation graph.
Pass `overwrite_family=true` only to deliberately
rebuild the entire family root from the requested sites alone.

```matlab
% Observation-only full PROMICE refresh: ALL stations found under
% data/verification/promice/hour, all available years, into the repo-root data/
% tree (data/eval + data/input). With no `case_ids` the driver defaults to the
% full station list. Use the firn/ablation production staging command above
% for the build_forcing=true PROMICE + MAR/MERRA/RACMO workflow.
icemodel.verification.setup.importPromiceSites( ...
   output_root=icemodel.internal.fullpath("data"), ...
   build_forcing=false, ...
   overwrite=true);
```

Each case carries THREE orthogonal descriptors, all single-sourced from
`promiceSiteCatalog(site)`:

- `surface_zone` is the glaciological zone ONLY (`ablation`/`percolation`/
  `wet_snow`/`dry_snow`/`accumulation`/`land`/`tundra`/`unknown`).
- `eval_target` names the model capability a case exercises
  (`seasonal_snow`/`bare_ice`/`firn`/`ablation`).
- `permafrost_zone` is the permafrost EXTENT class of the GROUND, ORTHOGONAL to
  `surface_zone` (`continuous`/`discontinuous`/`sporadic`/`isolated`/`none`/
  `unknown`): off-ice land/tundra sites carry the Obu et al. (2019) extent,
  ice-sheet/glacier sites carry `none`.

The KAN anchors are AUTHORITATIVE (KAN_L/M=`ablation`, KAN_U=`percolation`). The
non-KAN PROMICE stations and the ESM-SnowMIP sites carry **authoritative**
`surface_zone` + `permafrost_zone` values HARD-CODED from a data-driven analysis
(`test/interactive/site_classification/classify_site_facies.m`):

- `surface_zone` (ablation vs accumulation): MODIS end-of-summer **Bare Ice
  Extent** 2000-2018 frequency (ablation when bare ice in a majority of years;
  reproduces the KAN anchors KAN_L/M f_bare=1.00). Within the accumulation area
  the firn facies is refined: `dry_snow` for cold high interior (elev >= 2500 m,
  no bare ice), `percolation` where a **SUMup_2025** GrIS density profile lies
  within 15 km (firn observed), else `accumulation` (facies unresolved).
- `permafrost_zone`: **Obu et al. (2019)** permafrost-zone map, point-in-polygon
  on lon/lat (off-ice sites; ice-sheet/glacier -> `none`). Replaces the v1 Brown
  et al. (1997) source.

> The Obu zone polygons are read through `activelayer.readobuzones`, which lives
> in the external [`activelayer`](https://github.com/mgcooper/activelayer) dev
> repo and depends on shared helpers from
> [`matfunclib`](https://github.com/mgcooper/matfunclib). Both are placed on the
> path by `icemodel.dependencies` (called automatically by the test bootstrap);
> set `ICEMODEL_ACTIVELAYER` / `ICEMODEL_MATFUNCLIB` (or `ICEMODEL_PROJECTS_ROOT`)
> to use non-default locations. When these repos are absent the Obu-dependent
> classification skips cleanly. See the repo README "Optional external
> dependencies" section.

All cataloged sites now resolve. ZAC_A/L/U (A.P. Olsen / Zackenberg GlacioBasis
transect, NE Greenland) carry no installation coordinates in the CSV; their
lon/lat are sourced from the per-station L3 NetCDF global `latitude`/`longitude`
attributes (the same coordinates `readPromiceAws` reads), classifying them as
local-glacier AWS -> `ablation`, `permafrost_zone=none` (on glacier ice).
Uncataloged station ids fall back to the legacy elevation-band heuristic flagged
`classification="first_pass"`.

### SUMup firn cases (`importSumup`)

`importSumup` stages SUMup_2025 firn observation profiles (density `rho(z)`,
subsurface temperature `T(z,t)`, accumulation/SMB) as first-class
`firn_observational` cases under the `sumup` dataset family, enumerated by
`listcases`/`loadmanifest` alongside `promice` and the snow families. The eval
target IS bundled as a data-only `observations.mat` profile bundle (referenced
via `colocation.sumup.obs_file`, resolved to the case `evaluation_path`); the
manifest is **forcing-AGNOSTIC** - the co-located MAR met + RACMO Data forcing is
NOT bundled with the eval target and is recorded by source id only. A SUMup case
is a manifest (`case_id`, `case_type="firn_observational"`, `site`,
`surface_zone`, `eval_target`, `permafrost_zone`, `period`,
`comparison_variables`, the observation-bundle reference, and the colocation
record). KAN-co-located SUMup cases inherit the anchor classification from
`promiceSiteCatalog` (KAN_L/M=`ablation`, KAN_U=`percolation`).

SUMup cases are staged under `data/eval/sumup/` by default. The case ids and
on-disk case dirs match the PROMICE convention for co-located sites
(`kanl`/`kanm`/`kanu`, no redundant `sumup` prefix - the family folder is
already `sumup`), and the observation bundle is named `observations.mat`
consistently in the file tree and the manifest references. The SUMup family
shares those compact ids with PROMICE, so disambiguate with
`loadmanifest("kanl", dataset_family="sumup")`. The candidate adapter
(`candidateFromIcemodelOutput`) maps the SUMup profile-bundle comparison
variables (`density`, `subsurface_temperature`, `smb`) to a
`subsurface_profile_bundle` candidate. Density and temperature preserve every
model time-by-depth column as dated tall rows; behavior stays soft (diagnostic,
no hard gate).

SUMup has no demo-local subset path. Full interactive verification uses the
top-level `data/` archive; automated tests use isolated temporary roots and
provisioned `test/data` capabilities.

Use `matchObservations` when model development needs the actual paired values
rather than only the report overview or aggregate `comparecase` metrics:

```matlab
manifest = icemodel.verification.loadmanifest( ...
   "kanu", dataset_family="sumup");
loaded = icemodel.verification.helpers.loadArtifact( ...
   manifest.evaluation_path, "targets");
candidate = icemodel.verification.candidateFromIcemodelOutput( ...
   ice1, ice2, opts, manifest);
matches = icemodel.verification.matchObservations( ...
   loaded, candidate, variables=string(manifest.comparison_variables));
```

`matches.intervals` contains one row per completely covered half-open SMB
interval `[start_date,end_date)`. Model rate samples are integrated only when
their finite timestamps continuously cover the interval. `matches.profiles`
contains one row per common observed/modelled depth after pairing physical
profile identity and UTC calendar date; candidate values interpolate only
inside common depth support. Multiple candidate identities on one date require
an exact identity match and are never pooled. When a full model history has
multiple timestamps on that date, the matcher selects the exact observation
timestamp or the unique nearest timestamp, never repeated depth rows from
separate model states. Both tables retain source ids, variables, units,
observed values, and modeled values. Set `make_plot=true`
with `plot_kind="profile"` or `"interval"` for one focused plot built from those
same matched rows. `comparecase` delegates profile-bundle matching to this API,
so direct model-development use and suite metrics share one numeric contract.

Selected SUMup 2025 rows are deduplicated by exact variable-specific scientific
identity before datetime shaping. Location, depth/support, time/period,
reference, method, value, and uncertainty participate; `name_key`, resolved
`name`, elevation, and source measurement ids remain first-row provenance rather
than identity. Each artifact records flat `*_raw_rows`, `*_unique_rows`, and
`*_duplicate_rows_removed` metadata plus the same counts in its provenance note.
Repeated missing uncertainty compares equal on pre-R2026a MATLAB without
colliding with a genuine zero.

The fully implicit SUMup 2025 default catalog applies the source-audited canonical
case map: MIT owns SER-B's subset, TAS-L owns identical TAS-U, THU-U owns THU-L
and THU-L2 subsets, ZAC-U owns ZAC-A/ZAC-L subsets, and the externally staged
`research_site/humphrey` target owns Humphrey. A mapped loser is removed only
when its keeper is present; DY2 and Dye-2-long both remain because each adds
distinct identities. Explicit `points`, `case_ids`, or `anchors` bypass this
fixed release map.

### Additive compatibility discovery

`comparisonCompatibility` derives possible model-data, model-model, and
data-data comparison pairs from staged artifacts. It reads family manifests,
uses `whos -file` for MAT-file headers, and uses `ncinfo` for NetCDF headers so
large arrays are not loaded just to discover comparison possibilities.

This helper is additive. Manifest fields such as `comparison_variables`,
`observation_variables`, `forcing_sources`, `eval_sources`, and `colocation`
remain complete provenance declarations written by the importers. Compatibility
discovery can help future consumers choose possible pairings from staged
artifact contents, but it is not a reason to under-populate manifests.

### Artifact QA/QC gate

Run `auditArtifacts` before accepting regenerated verification previews. Pass
the paired roots and selected families explicitly so the audit cannot fall back
to a different staged tree:

```matlab
report = icemodel.verification.auditArtifacts( ...
   evaluation_data_root=fullfile(pwd, "data", "eval"), ...
   input_data_root=fullfile(pwd, "data", "input"), ...
   families=["promice", "esm_snowmip", "laugh_tests"]);
```

The default call is non-mutating and returns stable `artifacts`, `channels`,
`findings`, and `summary` structs. ESM-SnowMIP intentionally omits
`forcing_sources`, `colocation`, and `met_files`; audit and plotting therefore
reuse the standard runtime resolver to select its exact nested-or-flat met
artifact without changing the manifest contract. Findings distinguish `error`,
`blocker`, `warning`, and `placeholder`; intentional all-NaN placeholders are
recorded but not counted as observations. Seasonal 15-minute model met must also
carry the writer's source-gap provenance; the audit reports an error when a
derived channel contains fewer missing values than its native support requires.
For MAR userdata, the same source-light pass verifies the optional `SUH`/`SU`/
`RZ`/`ME` provenance, daily `/24` support, signed sublimation/deposition contract,
conservative unit/range bound, and the saved daily-`ME` versus hourly-`MEH`
ledger. The explicitly combined `subl_evap` and `refreeze_deposition` channels
plot in their own group; only hourly `SUH -> subl` shares a pure-sublimation
comparison panel with RACMO. Optional absence is recorded and does not make a
model-met artifact unready.
To persist results for Quarto or another report,
pass a generated output folder such as `report_dir=fullfile(pwd, "data",
"preview", "qa")`. That explicit side effect writes `artifact_qa.json` and
`artifact_qa.md`; it never rewrites manifests, met files, or userdata.

Generate the annual PROMICE snow-development handoff from the same explicit
roots after staging or promotion:

```matlab
icemodel.verification.setup.writePromiceSnowModelReadyYears( ...
   evaluation_data_root=fullfile(pwd, "data", "eval"), ...
   input_data_root=fullfile(pwd, "data", "input"), ...
   output_dir=fullfile(pwd, "data", "preview", "qa"));
```

The utility audits only complete years inside each mandatory MAR 3.11
precipitation-leg window. Practical rows require complete MAR model forcing and
an exact hourly PROMICE `snow_depth` grid with at least 95% finite samples and
one finite sample every day; strict rows require 100% finite snow depth. The two
CSVs, grouped Markdown, and self-hashed JSON use relative provenance and are
byte-deterministic for unchanged staged inputs.

### PROMICE forcing readiness and completion

PROMICE staging preserves corrected source values and source gaps. The canonical
model-forcing channels are `tair`, `swd`, `lwd`, `albedo`, `wspd`, `rh`, `psfc`,
and `ppt`. A native met artifact has `forcing_ready=true` only when all eight are
finite at every posting and its time coordinate is contiguous at the inferred
native cadence. `forcing_sources` still records that the artifact exists and is
selectable; it does not imply readiness. In particular, native PROMICE met keeps
`ppt` as an intentional all-NaN placeholder and therefore remains unready until
an explicit derived layer supplies precipitation or a different forcing leg is
selected.

`icemodel.verification.setup.metForcingReady` returns an inclusive table of
complete contiguous windows as its third output. Each reported window has 100%
coverage for every canonical forcing channel. No minimum duration is hidden in
the helper: the caller chooses a window long enough for the intended experiment.
An omitted timestamp splits a window just like an explicit missing value.
PROMICE, RetMIP, and IMAU importers diagnose the exact scalar window path
returned by `writemet`, including an exact or broader enclosing no-overwrite
reuse, rather than the pre-write request timetable. Their colocation leg stores
the resulting portable UTC records in `forcing_complete_windows`; empty,
singleton, and multiple records remain JSON arrays across additive rewrites.

Completion is separate from the source artifact:

- A complete-window run permits zero missing postings in any channel. Require all
  eight channels, a regular time coordinate, and a suitable window length;
  retain native source metadata and source uncertainty.
- A source-authoritative correction is not a gap fill. Use only a published
  corrected channel, preserve source flags and physical bounds, fall back to raw
  only where corrected data are missing, and record corrected/raw counts and
  source version.
- Empirical longwave is opt-in for `lwd` only when the selected source window has
  no finite observed longwave. It may span that whole window after validation
  against overlapping stations and physical units/ranges. Treat its model-form
  uncertainty as high and retain `lwd_estimated` plus its policy metadata.
- RCM substitution supplies `ppt` by default. Another wholly unavailable
  required channel needs source-specific validation. The maximum gap is the
  selected window within RCM coverage. Preserve rate integrals/signs, cadence
  support, bounds, overlap diagnostics, model/product id, remap method, and a
  per-channel substitution mask in a new derived artifact.
- Interpolation and statistical fills are not approved for canonical
  verification forcing. They need a separately accepted maximum-gap rule and
  held-out physical/conservation validation. They never rewrite the native
  artifact; ARIMA/seasonality alone is not evidence.

The current source audit separates coverage artifacts from source-reproduced
extremes. KAN_M's largest daily SWD mean (2012-06-03) uses 18/24 hourly samples;
KAN_U's (2019-06-05) uses 16/24. Coverage-aware daily plots require all 24
expected postings and therefore mask both partial-day artifacts. EGP 2024-06-12
instead has
24/24 corrected SWD and SWU samples: daily means are 562.213 and 472.461 W m-2,
and the corrected 15:00 UTC SWD is 1008.684 W m-2 versus 989.066 W m-2 raw.
Approximate solar geometry at the source latitude/longitude gives 799.925 W m-2
horizontal top-of-atmosphere flux at that posting. It is therefore a complete-
day, source-reproduced, physically review-worthy extreme, not a plotting gap
artifact. Preserve and flag it; do not clip it without an upstream quality rule.
The file labels raw `dsr`/`usr` as physical measurements and corrected
`dsr_cor`/`usr_cor` as model results, but supplies no per-sample acceptance flag;
"source-reproduced" therefore does not mean independently validated.
The reproducible inputs are `time`, `dsr`, `dsr_cor`, `usr`, and `usr_cor` in
`data/verification/promice/hour/{EGP,KAN_M,KAN_U}_hour.nc`, using the files'
`time:units`, latitude, longitude, pypromice version, and data-issues hash.

Coverage-aware plotting classifies a subdaily series as dense only when its full
unfiltered axis safely establishes a repeated cadence and at least one expected
day of postings. Each dense daily bin then requires the exact unique, monotonic,
on-grid native timestamps and finite values; row counts cannot conceal a
duplicate-plus-omitted or shifted posting. Partial boundary days remain missing.
Ordinary sparse amounts/means and explicit interval observations retain their
native timestamps and values, while a sparse rate with no defensible duration is
not converted into an invented daily total.

## Family Adapter Architecture

Family entry points are thin source adapters around shared control flow, not
independent staging implementations. Importers validate an optional paired
window, resolve their family catalog and staging roots, reuse or stage requested
cases, record provenance once, return through `runDatasetFamilyDryRun`, or
persist through `runDatasetFamilyImport`. `stageDatasetFamilyCases` owns the
common skip/error loop and `stageDatasetRcmForcing` owns delegated RCM work.
Research-site dry runs intentionally return before root/source resolution so a
metadata-only preview works on a clean machine; fixed Laugh-Test cases have no
caller-selected window. Those are source-contract differences, not alternate
manifest pipelines.

Fetch adapters retain family-specific file discovery and provenance, while the
shared registry helpers own selector extraction, validation, ordered status
construction, and retrieval banners where the upstream package model matches.
GC-Net remains station-aware and therefore does not use the simpler IMAU/RetMIP
product-registry adapter. Data-backed met builders all convert through
`icemodel.forcing.helpers.data2metCollection`; source readers/builders remain
separate because their variables, grids, and validation evidence differ.

`fillwithmissing` is a direct builder validation option. Verification importers
always request `fillwithmissing=true` so a staged native met artifact has the
canonical channel schema with unavailable channels represented explicitly as
NaN. Direct builder callers can pass `fillwithmissing=false` to require the
source itself to satisfy the complete met contract. It is deliberately not an
importer option because changing it would make the persisted schema depend on
which family wrapper was called.

## Source Catalogs and Staged Schemas

Family source catalogs are inventories, not competing artifact schemas. Their
names make that distinction explicit:

- `promiceSiteCatalog` and `esmSnowmipSiteCatalog` describe source sites and
  source-specific coverage/classification fields;
- `imauSiteCatalog` and `researchSiteCatalog` describe the source anchors those
  importers can stage;
- `retmipCaseCatalog` remains a case catalog because RetMIP protocol aliases,
  windows, and source associations describe protocol cases rather than only
  physical sites.

The catalogs deliberately retain source-specific fields instead of padding a
union struct with fields that do not apply. Importers normalize staged snow
cases through `makeCaseManifestEntry` and staged firn observational cases
through `makeFirnCaseManifestEntry`; those factories and their field-name
helpers are the canonical persisted schemas. Catalog selectors only choose
source rows and must not define a new on-disk contract.

## Support Namespaces

- `helpers` contains normal workflow helpers for path discovery (`evaluationDataRoot`,
  `inputDataRoot`, `esmRuntimeMetFiles`), manifest reads, artifact loading, candidate
  resolution, metric schema definition, the per-run markdown report writer
  (`writeRunReport`), the per-site default window (`default_smoke_window`), and
  per-site default window (`default_smoke_window`). The standard-contract
  opts builder used by `runIcemodelSnowCandidate` is
  `icemodel.test.helpers.setModelOptsForCase`, which accepts both formal-case
  rows and verification manifests via input dispatch.
- `setup` contains the consistently named family source catalogs listed above,
  their shared strict site-id selector (`selectSiteCatalogEntries`), and the
  canonical staged-case factories. RetMIP keeps alias-aware case selection in
  `retmipCaseCatalog`; PROMICE retains its documented first-pass fallback and
  ESM-SnowMIP retains its scalar site lookup because those semantics differ.
- `namelists` contains canonical selector lists for dataset families, case ids,
  case types, surface zones (`surfacezone`, the per-case physical-regime
  vocabulary stamped onto case manifests), the ESM-SnowMIP site-name namelist
  (`snowmipsite`), and the Laugh-Tests case-id namelist (`laughtests`). `caseid` dispatches uniformly
  across families using these per-family namelists. The richer per-site
  ESM-SnowMIP catalog query helper lives at
  `icemodel.verification.setup.esmSnowmipSiteCatalog`.
- `validators` contains argument-block validators that consume the namelists,
  including `mustBeSnowmipSite` for per-site builders.

## Data Contract

Each dataset family has one `manifest.json` under:

`data/eval/<dataset_family>/manifest.json`

The per-case folder layout is split by `case_type`:

- **Analytical families** (`laugh_tests`; `synthetic_process`) bundle a computed
  reference: each case folder stores `evaluation.mat` (the staged targets) and
  `reference.mat` (the analytical / frozen-SUMMA solution the case is gated
  against). This is the reference, not a smoke copy, so it is KEPT.
- **Observational families** (`esm_snowmip`/`esm_site`, `promice`/`sumup`/
  `firn_observational`) are FORCING-AGNOSTIC: the case folder stores one
  data-only `observations.mat` bundle (the eval target). The manifest is
  forcing-agnostic - it records which forcing/eval sources are available (by id,
  informational only), but the forcing is NOT bundled and NOT stipulated, so any
  forcing usable at runtime without rewriting `observations.mat`. No bundled
  `reference.mat` smoke copy is written — the default candidate, with no model
  output supplied, falls through to the soft diagnostic lane. Forcing always
  lives separately under per-source subfolders `data/input/met/<source>/` and
  `data/input/userdata/<source>/` (standard icemodel naming via
  `writemet`/`writeuserdata`), never in the eval folder. (PROMICE demo fixtures
  staged before this contract carry no
  `observations.mat`; the workflow functions fall back to reconstituting the
  PROMICE-obs target from the per-year userdata files those manifests declare.)

Manifests keep case paths relative to the dataset-family folder. Normal workflow
functions resolve those paths to absolute paths at read time. For esm_snowmip and
freshly staged promice the `observations.mat` bundle is referenced from
`evaluation_file` (and `observation_variables.obs_file` for esm_snowmip); SUMup
references it via `colocation.sumup.obs_file`. `reference_file` is empty for all
observational families.

### Target schema variants

Two staged-target shapes are supported (`evaluation.mat` for the analytical
families, `observations.mat` for the observational families):

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
`peak_*` columns above.
