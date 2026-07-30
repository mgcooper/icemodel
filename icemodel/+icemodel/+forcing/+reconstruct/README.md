# icemodel.forcing.reconstruct — provenance-preserving gap filling

Build a separately named, provenance-preserving forcing product from incomplete
native station data. The main entry point is `fillPromiceStation`.

## Getting started

Run MATLAB from the repository root and put the package on the path:

```matlab
addpath("icemodel")
result = icemodel.forcing.reconstruct.fillPromiceStation("kanm");
```

The default call reads the staged native PROMICE, donor, proxy, and MODIS
artifacts from `data/`, then writes:

| Output | Default location |
|---|---|
| Filled 15-minute met | `data/input/met/promice_filled/` |
| Station plan, audit, seam QA, and native-sensor findings | `data/preview/qa/gapfill/plans/` |
| Per-year readiness ledger | `data/preview/qa/gapfill/ledger/` |
| Selection/evaluation split | `data/preview/qa/gapfill/splits/` |

The native PROMICE met file is an input and is never modified. Before the first
run against older staging, refresh its manifest fingerprints:

```matlab
icemodel.verification.setup.refreshPromiceMetIdentities()
```

That command hashes the staged met files already on disk; it does not fetch or
rebuild source data. See the
[verification README](../../+verification/README.md#refresh-promice-met-fingerprints)
for custom roots and the full-cohort workflow. To build a report after the
products exist, see the
[report README](../../+verification/+report/README.md#promice-gap-fill-report).

Common overrides are explicit:

```matlab
opts = icemodel.forcing.reconstruct.setopts();
result = icemodel.forcing.reconstruct.fillPromiceStation("kanm", ...
   opts=opts, donor_sites=["kanl", "kanu"], ...
   use_ktransect=true, use_gcnet=true);
```

Use `met_dir`, `modis_dir`, `out_dir`, and `qa_dir` together when operating on
an alternate data tree. The selected native met directory determines the
enclosing data root; reconstruction never falls back to another tree.

## Contract at a glance

This family-generic engine produces `promice_filled`, the canonical runnable
PROMICE forcing. Native PROMICE is incomplete for most station-years and is
retained unmodified as the provenance source.

- The approved contract is [`POLICY.md`](POLICY.md) and the DesignSpec
  `2026-07-23-promice-gap-filling-and-ktransect`.
- Every filled sample carries a `uint8` provenance code and every contiguous
  segment has an audit row.
- `setopts` is the single source for scalar knobs, channel lists, and proxy
  source mappings. Dedicated functions own per-channel bounds, admission caps,
  precipitation names, and bucket edges.
- Every required non-precipitation output is planned. Precipitation is excluded
  from the statistical method set at option validation and at the direct
  planner boundary.
- No product ships without the provenance path required by policy.

## Role contract

Engine functions take target, donor, and proxy SERIES as role-typed inputs
(timetable + site metadata) — never family names. The same code fills a
PROMICE station from K-transect donors or a K-transect station from PROMICE
donors (the bidirectionally demonstrated role-reversal acceptance).

## How the engine is organized

### Validation and admission

Methods admit only through this validation harness:

- `gapCensus` — per-channel missing-run census (joint-core in-record bounding,
  daylight-only shortwave option, `bucketEdges` strata assigned by the
  right-closed `gapDurationBucket` convention).
- `validationSplit` — persisted whole-year selection/evaluation split
  (a schema-valid manifest wins only while it remains a disjoint, complete
  partition of the current record years).
- `syntheticMissingness` — blocked synthetic gaps drawn from the real
  run-length distribution, inserted only into observed spans; time-sorted.
- `validationMetrics` — bias/RMSE/correlation/within-gap observed spread and
  anomaly variability ratio/coverage (admission requires measurable spread
  evidence for every included gap), boundary jumps vs the shared `stepScale`, scalar and relational bounds
  via `physicalValidity`, sigma calibration, typical stratum magnitude,
  and filled-sample provenance accounting. Planner calls mask every held-out
  draw before deriving the native seam scale.
- `admissionGate` — the approved per-channel thresholds (bias caps, ≥10%
  common-support RMSE improvement over persistence for ≤6 h or station
  day-of-year climatology otherwise, zero bound violations, coverage
  usefulness floor, and bounded variability ratio — see POLICY.md).

### Method tiers

- `fillShortGaps` — tier 1 bounded interior interpolation (≤6 h except
  the observed-only-supported ≤9 h SWD and RH caps; SWD is CSI first and
  flux-linear only on the post-final pass; observed-only-supported linear
  albedo ≤30 h; never independent swu).
- `toaIrradiance` / `physicalValidity` — shared solar geometry and physical
  checks, including `swd <= max(1.05*TOA, 5 W/m²)` plus the approved
  civil-twilight allowance classified from the complete posting interval,
  matching staging and darkness, and
  `swu <= swd` (tier-1 CSI and the swd darkness pre-pass must agree on
  what "meaningful sun" means).
- `fitDonorTransfer` / `applyDonorTransfer` — tier 2 per-season transfer
  with the monotone piecewise knot hyperparameter, lag search, elevation
  adjustment (`elevationAdjust`), and refusal beyond the allowed
  extrapolation fraction of the fitted donor range (knobs in `setopts`).
  Overlap is recorded in elapsed hours by integrating aligned target
  intervals, so support-held coarse postings are not multiplied by row count.
  SWD lag search, fitting, elevation brackets, and application operate in
  station-specific clear-sky-index space and return flux on the target TOA
  geometry.
- `fitProxyCalibration` / `applyProxyCalibration` — tier 4 overlap bias
  calibration (additive for state channels; multiplicative for shortwave and
  wind speed, preserving their shape and wind's nonnegative support).
  Shortwave overlap is screened by target-station TOA rather than proxy
  magnitude; fitted proxy
  corrections persist in the station plan for calibrated last-resort use,
  whether or not the candidate wins an admitted stratum, but zero-overlap
  corrections are ineligible for both competition and last resort.
  `lwdEstimator` is
  the calibrated empirical lwd candidate.
- `climatologyFill` — day-of-year × exact posting-time lookup climatology
  (long-gap admission baseline and last-resort method), preserving hourly,
  half-hourly, or quarter-hourly target cadence.
- `deriveUpwardShortwave` — the ordered algebraic `swu = albedo * swd`
  fill, applied only after both operands finish every reconstruction tier
  and stamped with its dedicated provenance code.
- `partitionPrecipitation` — canonical total-to-rain/snow split at the
  configured air-temperature transition; phase channels are never
  independently reconstructed.

### Orchestration and production

- `reconstructSeries` — zero-fills missing no-sun swd samples first (the
  policy's darkness rule: a known value with its own provenance code,
  decomposing multi-day outages into daylight fragments), then composes
  ADMITTED methods per season-contained missing segment (tier-1 refuses a
  cross-season bridge except SWD within its evidenced nine-hour cap, which
  borrows an adjacent anchor across, or itself straddles, the calendar edge;
  RH has the same calendar-label exception only within its
  evidence-backed nine-hour cap;
  later tiers split at the boundary), blends
  excess anchored-seam mismatch across each contiguous residual time segment
  (policy rule
  7), enforces physical bounds, stamps
  `provenanceCodes`, and emits one audit row per contiguous segment through
  `auditSegments`, including per-run refusal reasons and a stable context id
  joining planned fills to their exact persisted fit/validation record. A
  method cannot be applied beyond the longest held-out duration supporting
  its admitted bucket. Fitting/admission happen upstream; this only composes.
- `blendSeams` — the shared linear taper used by both tier-1 interpolation
  and later admitted methods; physical-bound failures still refuse a fill.
- `fillTwilightClimatology` — final SWD-only repair for one residual
  civil-twilight posting beside known all-interval darkness; reuses the
  station day-of-year/posting median only with native support >=30 and the
  50 W m-2 twilight ceiling, with dedicated twilight-climatology
  provenance and exact audit; the later seam pass cannot replace or relabel
  these separately validated samples.
- `fillShortGaps` also performs the opt-in post-final D-45 edge repair:
  2--9 still-sunlit postings beside one policy-known darkness zero use the
  opposite daylight anchor's CSI over target TOA. The configured darkness
  threshold and shared SWD physical-validity ceiling still apply; deep
  darkness and one-posting D-44 cases remain separate.
- `fillPromiceStation` reruns that same post-final closer after guarded
  hourly or half-hourly postings return to the delivered 15-minute axis
  (D-47). It then derives newly available `swu` from final `albedo * swd`
  and reconciles the audit, so a coarse low-sun refusal cannot leave four
  fillable quarter-hours or a stale dependent channel in the product.
- `smoothShortwaveSeams` — post-final SWD boundary QA against untouched
  native steps from the same season and fine solar-elevation bin. It repairs
  at most one reconstructed posting beside an empirical outlier, preserves
  native/raw/clamped/darkness/twilight-climatology values, accepts smooth
  same-direction ramps within one coarse solar regime and season, and
  excludes the one-posting neighborhood of physical solar-regime edges.
  Zero native-reference steps do not themselves cause `review`: a record
  without a method boundary needs no boundary reference, and a boundary
  independently proven continuous by the following rules remains supported.
  A nonzero boundary without matching reference support returns `review`
  unless its step is independently supported by an adjacent native-only
  slope, the unavoidable native-anchor bridge slope, or an exact linear
  `bounded_interp` run between finite anchors. Exact-bridge support is
  boundary-local: the boundary must independently pass the season, fine-bin,
  coarse-regime, and physical-edge-neighborhood guards even if the enclosing
  bridge begins in another regime. Fine-bin native evidence is preferred; a
  sparse bin may use only same-season, same-coarse-regime native evidence
  meeting the same sample minimum. A nonlinear or arbitrary synthetic
  neighbor never supplies support; an exactly continuous provenance change
  is not a seam.
- `blendFallbackSeams` — applies that same taper to proxy, precipitation,
  and constant fallback segments using the frozen native step scale.
- `stationMethodPlan` — the per-station selection experiment: the
  policy geometry gate (`setopts`), deterministic selection draws created
  only from years with jointly finite core-channel support before fitting
   and withheld from every candidate, selection-only
   climatology, union-masked persistence anchors, persistence-relative ≤6 h
   and climatology-relative longer-gap
  admission on an adequate common candidate/baseline sample set with recorded
  denial reasons, evaluation of every geometry-eligible
  channel donor before the distinct-donor cap is applied by held-out RMSE
  independently within each season × duration stratum,
  the additional requirement that wind donors beat calibrated MAR wind in
  the same station/season stratum,
  overlap-fitted elevation corrections (documented fallback when
  underdetermined), and evaluation-year regrade of admitted winners. A
  one-year record with no disjoint evaluation year admits no production
  method. Support-held donor values stop after the donor's final cadence
  interval. SWD census and synthetic draws use the same central meaningful-sun
  TOA threshold, so polar-night zeros cannot enter validation ranking.
- `fillPromiceStation` — the production driver. Its responsibilities are:

  - **Inputs and donor pool.** It considers inventory-selected PROMICE
    neighbors, hourly-aggregated K-transect cases, and GC-Net
    origin-observed samples. Direct planner calls and the production driver
    exclude the target station and its canonical aliases. Approved
    K-transect `promice_neighbor` relations remain distinct donor hypotheses.
    K-transect geometry is cross-checked between the byte-pinned observation
    artifact and its manifest before donor gating.
  - **Proxy discovery.** It loads only the canonical per-site RCM met selected
    by the acceptance window: `data/input/met/mar3.11` and
    `data/input/met/merra2`, using versioned storage names rather than short
    tokens. `selectedDataRoot` normalizes canonical per-source, canonical flat,
    compact, and custom met paths once. Proxy and family-donor discovery then
    use that same root and the shared per-source-first/flat-fallback order.
  - **Residual outages and precipitation.** One source-consistent proxy covers
    each whole residual outage. Total precipitation adopts MAR first, then
    MERRA-2, with seam control and preserved PROMICE corrected-liquid
    observations. Phases come from native components, exact complements, or
    the proxy's rescaled split; reconstruction never partitions by
    temperature. The runtime `precip_phase_source` option resolves
    phase-unknown samples (POLICY A10/D-18).
  - **Native honesty.** Legacy-filled staged albedo is masked to raw-source
    finite support before planning. Raw PROMICE albedo support uses the native
    closed `[0,1]` source contract, not the narrower post-fill bounds.
    `flatRunScreen` also identifies only sustained, multi-channel
    burial/rime signatures; implicated channels are masked in the working copy
    before target training or donor use and recorded in the plan sidecar while
    the staged artifact remains unchanged.
    `boom_height` gaps remain missing in the product: no donor or proxy supplies
    geometry. The runtime A3 chain in `loadmet` resolves measured →
    interpolated → nominal geometry.
  - **Cadence.** Guarded PROMICE 15-minute staging is collapsed to hourly
    source postings for planning and reconstruction. Observed values and their
    provenance remain exact held copies over the original four-sample support;
    filled values use the policy-approved mean-preserving disaggregation. The
    canonical runtime artifact is published only at 15 minutes, and runtime
    discovery rejects every other `opts.dt`. When several staged windows
    exist, saved timetable coverage—not MAT-file size—selects the widest.
  - **Readiness.** Provisional `unfilled` audit rows are reconciled before two
    per-station-year verdicts are written: `verdict_icemodel` grades the A5
    seven-channel forcing set, and `verdict_snowmodel` additionally grades
    total precipitation or snowfall availability. Boom geometry grades
    neither verdict (POLICY A3/A5).
  - **Publication and verification.** Artifacts ship through `stampMetadata`
    and `writemet`. Every native target and PROMICE donor must match its
    upstream staging-manifest byte identity before training. Each written
    station records size and SHA-256 identities for the native, filled,
    plan/audit, readiness, and proxy-window artifacts consumed by the report.
    Paths are relative to the selected data root, so moving the complete staged
    tree does not invalidate provenance.
  - **Runtime binding.** The filled artifact records the engine version, policy
    SHA-256, donor list, and exact planned channel set. Runtime loading verifies
    the readiness ledger and every configured filled met file against that
    transaction's producer manifest before trusting the site/product identity,
    canonical provenance registry, or typed provenance arrays.

## Hard rules

- Never overwrite a finite native sample.
- Reconstructed donor samples are never donors (per-channel GC-Net origin
  masks; a channel without an origin flag is ineligible).
- In-sample overlap statistics are never admission evidence.
- Proxy coverage is read from the selected staged artifacts, never assumed from
  a global product date. The common MERRA-2 window is 2012–2018 and MAR commonly
  ends in 2019, while some sites have later staged products. Outside each
  station's actual proxy window, `verdict_snowmodel` honestly reports missing
  precipitation or snowfall input; precipitation never gates
  `verdict_icemodel` (POLICY A5).
- PROMICE runtime observation heights come from the time-varying
  `boom_height` where measured; unresolved samples resolve through the
  runtime fallback chain (measured → interpolated → nominal 2.6 m,
  POLICY A3) and never block readiness or loading.
- Tests: `test/unit/test_reconstruct_{harness,engine,production}.m`.
