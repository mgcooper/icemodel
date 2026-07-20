# MAR density profile contract

MAR density profiles are an optional verification sidecar. They are not model
meteorology, scalar userdata, observations, or a forcing-readiness requirement.

The public product is source variable `RO1`, sampled only on explicit requested
UTC calendar dates. Each selected snapshot preserves the native 18-level
`OUTLAY` coordinate and `OUTLAY_bnds`, identifies rows by `profile_id` and source
`datetime`, uses positive-down depth in metres and density in kg m-3, and labels
the values as MAR model output. A reader never substitutes a nearest date or
extrapolates outside the yearly file. Reduced sources, missing years, non-ice
cells, and invalid profiles produce explicit omission statuses.

The reader validates source metadata before reading or stamping public values.
`TIME` must have only its `TIME` axis, `OUTLAY` only its `OUTLAY` axis, and
`LON`, `LAT`, and `SRF` must share exactly one native X axis and one native Y
axis (full RUH2 names include their index extents). `RO1` must contain exactly
those grid axes plus `TIME` and `OUTLAY`. Native `OUTLAY` units must be `m`;
native density accepts only MAR's equivalent `kg/m3` or `kg m-3` spellings and
is stamped canonically as `kg m-3`. Missing, duplicated, unknown, or extra
public axes and missing or wrong units return `invalid_source`.
Only a source that genuinely omits the RO1/OUTLAY profile product is optional;
once RO1/OUTLAY exist, missing or unusable `TIME`, missing bounds, or malformed
bounds dimensions are authoritative public-schema failures and also return
`invalid_source`.

Native daily time is decoded from `TIME` plus its declared elapsed-time units
and, in full RUH2 files, checked exactly against `YYYY`, `MM`, `DD`, `HH`, and
`MIN`. All 40 full years place the daily profile at 12:00 UTC. The float32
packed `DATE(YYYYMMDDHH)` field cannot represent ten decimal digits: it
quantizes and sometimes duplicates adjacent dates. It is never decoded or used
for profile selection. Exact components validate but do not replace the native
coordinate; a missing `TIME`, partial component set, or disagreement with
`TIME` fails closed.

`DZSN1` thickness and `ROSN1` density remain internal source QA. Their positive
activity masks must agree. Active source indices are reversed from MAR's
numeric bottom-to-surface storage order to construct diagnostic surface-down
midpoints, and active `DZSN1` thickness must match permanent-ice `SHSN3` sector
1 within 2e-5 m. Interpolation of that diagnostic reconstruction to `OUTLAY`
reports bias/RMSE/max mismatch only; it never changes `RO1`. Dynamic profiles
must not become a public artifact until independent source documentation
confirms the vertical-order interpretation.

Dynamic reads are equally fail-closed but remain nested: `DZSN1` and `ROSN1`
must contain exactly the native X/Y grid, `SNOLAY`, and `TIME` axes; `SHSN3`
must contain exactly the native X/Y grid, `SECTOR`, and `TIME` axes, with at
least sector index 1. `DZSN1` and `SHSN3` require native `m`; `ROSN1` requires
one of the two accepted native density spellings above. The reader explicitly
selects permanent-ice `SECTOR=1`. A dynamic schema, unit, or read failure is
reported as `source_read_error` and cannot invalidate an otherwise valid
public `RO1` snapshot.

The low-level reader currently accepts nearest-cell sampling only: its native
`[row column]` index is an exact grid-cell identity, not proof of a natural-
neighbour or polygon collapse. A later staging caller can reuse
`buildMarData`'s existing `grid_start`, `grid_count`, and `method` metadata when
`grid_count=[1 1]`; non-nearest legs must remain explicitly profile-unavailable
until the shared colocation collapse is applied to the profile arrays.

The low-level contracts are:

- `icemodel.forcing.helpers.readMarDensitySnapshots`: dimension-name-aware,
  exact-date RO1/OUTLAY snapshot selection with profile and provenance output.
- `icemodel.forcing.helpers.marDynamicProfileQa`: dynamic-layer activity,
  thickness, ordering, and diagnostic-only RO1 reconstruction checks.

SUMup staging attaches the optional sidecar through the existing MAR colocation
leg after primary RCM staging. The leg uses `model_output_files`,
`model_output_format`, and `model_output_variables`; it never enters
`data_files` or changes forcing readiness. The stable per-case sidecar filename
supports additive surgical refreshes: valid requested `profile_id` groups are
replaced, unrequested groups are preserved, and empty/failed refreshes do not
delete an earlier valid product. Comparison and plotting group by both
`profile_id` and date, pair exact UTC calendar dates, and compare MAR with SUMup
only over common midpoint depth support.

Replay the source and grid audit with
`icemodel.verification.setup.auditMarSemanticsAndGrid()`.
Its durable variable, representative-site, and all-PROMICE mapping ledgers live
under `data/preview/qa`.

Replay the bounded three-case staging, comparison, plot, and artifact-audit proof
with `icemodel.verification.setup.runMarProfileBundleProof()`. It writes only the
compact profile proof table and three review figures under `data/preview/qa`; its isolated
temporary stage is removed after the checks finish.

Run the complete focused gate with
`matlab -batch "addpath('icemodel'); icemodel.verification.setup.runMarProfileValidation()"`.
The gate runs
the dedicated unit files, selected shared-path regressions, and both source-backed
replays in dependency order.
