# icemodel.forcing

Forcing and evaluation-data builders for icemodel. Each gridded/AWS source has
a `read*` (raw -> canonical channels), a `build*Met` (met-contract timetable),
and a `build*Data` (evaluation/userdata timetable) function.

The conservative (area-weighted) polygon remap in
`icemodel.forcing.helpers.remapPolygon` needs the external
[`exactremap`](https://github.com/mgcooper/exactremap) dev repo on the path. It
is wired by `icemodel.dependencies` (called by the test bootstrap); set
`ICEMODEL_EXACTREMAP` (or `ICEMODEL_PROJECTS_ROOT`) to use a non-default
location. When `exactremap` is absent the conservative-remap tests skip cleanly.
See the repo README "Optional external dependencies".

## Model-met cadence and source gaps

`icemodel.forcing.helpers.resampleMetTimestep` owns model-met resampling.
Data-backed met builders use the shared
`icemodel.forcing.helpers.data2metCollection` adapter; it preserves scalar or
cell input shape, delegates the canonical channel conversion to `data2met`, and
then applies the same resampling helper used by `writemet`. The writer emits
15-minute files by default; `dt_out=""` is the explicit native-cadence override.
The proven Samimi Dye-2 override is named `_30m.mat` and is discoverable by the
standard runtime resolver; no other new native cadence is implied.
Canonical source rows are interval-start means, so upsampling uses zero-order
hold over `[t,t+source_dt)` rather than treating them as point samples. This
preserves every interval integral, includes the complete final source interval,
and leaves explicit NaNs and omitted-time gaps missing. Gap steps must span whole
native intervals. Two-row sources support a finite 15-minute-multiple cadence up
to one hour; longer lone steps fail closed because they cannot distinguish native
cadence from an outage. Saved metadata records the
exclusive support end and source-derived missing counts for
`icemodel.verification.auditArtifacts`. For yearly output, the full native source
is validated/resampled before calendar slicing so Jan 1 gaps remain explicit.
Compact per-year summaries record overlapping source rows, intersecting source
gaps, missing counts, and support bounds; each saved file replaces the parent
aggregates with its matching summary. Guarded 15-minute input must carry valid
summaries and constant cadence blocks or yearly writing fails closed. A native
`dt_out=""` window write remains an exact no-op, while guarded yearly slices keep
their rows and localize only the aggregate provenance. Source-row and gap counts
are artifact-local overlap counts; a source interval crossing Jan 1 can belong to
both adjacent summaries, so yearly counts are not necessarily additive.

`writeuserdata` owns the analogous public userdata boundary. It defaults to
hourly output while leaving already aligned hourly products unchanged. Finer
native averages are aggregated into clock-hour means; coarser products are
linearly interpolated to hourly support; wind direction is handled on the unit
circle in both directions. Builders remain source-native, and `dt_out=""` is
the explicit writer override for a deliberately native-cadence artifact. Saved
metadata records the source cadence, output cadence, and resampling policy.
Hourly userdata retains the legacy suffix-free name; native variants use an
explicit cadence suffix such as `_30m`, including yearly files. A verification
manifest's selected `data_files` are resolved into `opts.userdatafname`, which is
authoritative in `loadmet`; this prevents legacy name discovery from choosing an
hourly sibling when a native Data artifact was selected. A stale recorded path
therefore fails at that exact path instead of silently falling back to a sibling.
Manifest-selected `met_files` are also authoritative in verification runtime
setup. A recorded file is selected only when its saved `met.Time` encloses the
requested run. Otherwise, the existing recorded files that meet the run are
sorted by saved support and must form one non-overlapping, gap-free list at one
uniform saved cadence. Missing, overlapping, gapped, mixed-cadence, or
noncovering explicit lists fail with their recorded paths before `configureRun`
can invoke name-based discovery. Runtime `dt` comes from the selected saved Time
cadence, not an earlier list member, filename suffix, or conflicting caller
override. An explicit caller `dt` may repeat that authoritative saved cadence,
but a different value fails with the recorded met paths; legacy cases without an
explicit staged met list keep caller-override behavior. A 30-minute run that swaps
channels from another source
prefers 30-minute staged met, then falls back to compatible 15-minute and hourly
met in that order.

Window-named `writemet` and `writeuserdata` outputs are additive across repeated
calls. A narrower request reuses an enclosing current artifact. When a newly
written wider window strictly contains older windows in the same exact
site/source/cadence/provenance class, the writers warn and remove only those
stale shorter files. Concrete source/product, schema, sampling-method, and point
metadata must also agree for reuse or pruning; an exact-name conflict requires
`overwrite=true`. Different sources, sites, cadences, identities, overlapping
windows, and enclosing windows remain unchanged. Production `method` and repaired
`sample_method` metadata are one documented sampling-identity alias. Both writers
stamp the actual uniform payload cadence as `artifact_cadence_seconds`; correctly
sampled legacy files without that top-level field are checked from their saved
timetable, so a forged `_15m` suffix cannot authorize reuse or pruning.
Both writers validate an existing exact target before considering a compatible
broader window. Model-met writes return a compatible exact target because runtime
met resolution also gives that name precedence; userdata may still return the
widest compatible enclosing file after the exact target has passed validation.
The pure scalar identity comparison is shared with manifest merging so family,
source, source id, native `source_family`/`station`, product, DOI/bundle DOI,
schema, relationship, and the documented method alias cannot drift between reuse
and merge boundaries. Native producer keys compare by their exact production
spellings; missing legacy values remain unknown-compatible, and no cross-key
alias is invented beyond `method`/`sample_method`.
`writemet` completes resampling, cadence derivation, channel validation, and
guarded-year provenance checks before creating its per-source directory. Rejected
input therefore leaves no empty artifact directory. `writeuserdata` already uses
the same prevalidation-before-mutation ordering.

## GC-Net / Vandecrux RetMIP surface forcing

`buildGcnetVandecruxData` and `buildGcnetVandecruxMet` read the Vandecrux
GC-Net gap-filled surface/SEB files used for RetMIP Dye-2-long and Summit
(`DYE_2_surface.nc`, `Summit_surface.nc`). They map `SRin/SRout/LRin/LRout`,
`Ta_2m`, `RH_2m`, `WS_10m`, `Pres`, and the source mass-balance terms onto
canonical icemodel channels, attach userdata location metadata, and preserve the
source variable/provenance map in `Data.Properties.UserData`.

The surface files are hourly by row; their fractional-day time coordinate drifts
inside leap years and jumps at several year boundaries, so the builder
regularizes the axis to a continuous hourly sequence from the first source
sample.

`buildGcnetVandecruxFirnTemperature` reads the companion
`*_T_firn_obs.nc` firn-temperature observation products for Dye-2 and Summit.
It preserves the source `T_firn`/`Depth` names and DOI provenance while exposing
canonical `subsurface_temperature` [K] and `depth` [m] matrices shaped
time-by-level. The reader windows NetCDF payload reads after reading the small
time coordinate, so tests and staging checks do not need to load the full source
array unless the caller requests the full product.

Source policy is explicit: `LRin` is source-filled regional-climate-model
longwave, not an observed GC-Net longwave sensor. The builder records both the
settled RetMIP/Vandecrux HIRHAM5 context and the local surface-package RACMO2.3p2
attribute text. The files carry snowfall but no rain channel; snowfall is
converted from source timestep amount to `snowf` [m s-1], `rainf` is all-NaN,
and `ppt` remains all-NaN in native met products so downstream substitution can
provide total precipitation without confusing absence with zero.

## IMAU hourly AWS forcing

`buildImauHourlyData` and `buildImauHourlyMet` read the Van Tiggelen et al.
PANGAEA hourly S21/S22/S23 tables. Corrected air temperature, relative
humidity, wind speed, pressure, radiation, albedo, surface temperature, boom
height, and surface-height change are mapped into canonical Data/userdata
channels. The daily 19-station SEB product is QA/provenance for importers, not
a first-pass case inventory.

An IMAU observation refresh with `build_forcing=false` retains only a compatible
prior native artifact contract (`met_files`, `data_files`, readiness, and artifact
metadata). Compatibility requires the same producer family/station and any finite
recorded point. Fresh hourly source paths, observation window, evaluation
reference, and daily-QA provenance remain owned by the current observation build.
When producer identity changes, the refreshed manifest unlinks the prior native
references without deleting their files and marks the leg for an explicit forcing
rebuild.

The hourly files do not provide a precipitation split. The builder records this
policy and writes `rainf`/`snowf` as all-NaN placeholders; native met products
therefore carry all-NaN `ppt` through `data2met(fillwithmissing=true)`.

## MAR native daily runoff and SMB

`buildMarData` and `buildMarMet` constrain hourly `RUH`/`SMBH` with MAR's
native daily `RU`/`SMB`. The MAR archive carries two
surface sectors: sector 1 is permanent ice and sector 2 is tundra. Point builds
select the sector from `SRF` (`4` means permanent ice); ice-masked polygon
builds use sector 1. Native daily `RU` is the delayed runoff product, while
`RU2` is explicitly labelled "without delay" in the NetCDF metadata.

Raw-source comparison across current PROMICE pixels shows that summed `RUH`
matches delayed `RU` on every selected-pixel May--September day but often
follows no-delay `RU2` outside that window. For each complete UTC day and
channel, the builder preserves the native hourly structure when its sum agrees
with the daily reference (absolute tolerance `1e-5 mWE/day`, relative tolerance
`1e-4`). Missing, incomplete, or inconsistent hourly days use daily/24; a
missing daily reference preserves the hourly source as explicitly unverified.
Partial artifact-boundary days receive the known daily rate on available rows.

MAR hourly fields do not declare a `999` fill value. `readMar3p11` therefore
masks only signed fill magnitudes of at least `1e30`; it does not apply the
legacy positive `>=999` cut, which removed valid no-delay runoff and one side of
paired SMB pulses. Builder and repaired artifacts record the
`daily_constrained_hourly` method, sector, tolerances, and a compact per-day
status/reference ledger (`preserved`, `replaced`, or `unverified`) in
`Properties.UserData` / `artifact_metadata`. An artifact created by the former
constant-daily policy requires full restaging to recover discarded RUH/SMBH;
metadata-only repair never relabels it as hourly-preserving.

Intentionally reduced legacy/test sources that omit either `RU` or `SMB`
retain their available hourly `RUH`/`SMBH` values. They are not presented as
corrected: metadata records `mar_qc_status=not_applicable` and the explicit
`hourly_RUH_SMBH_retained_native_daily_unavailable` fallback token.

## MAR optional mass-balance diagnostics

The 1980--2019 RUH2 archive declares `SUH`, `SU`, `RZ`, `ME`, and `MEH` in all
40 yearly files. Their NetCDF definitions are not interchangeable:

- `SUH` is hourly sublimation (`mmWE/h`) and maps to canonical `subl` after the
  reader's `mmWE`-to-`mWE` conversion. Positive values are surface mass loss;
  negative values are deposition.
- `SU` is daily "Sublimation and evaporation" (`mmWE/day`) with the archive's
  two-sector axis. It selects the same point/polygon sector as `RU`/`SMB`, is
  divided by 24 and previous-held over its UTC day, and maps to the deliberately
  combined `subl_evap` channel. It is not relabelled as pure `subl` or `evap`:
  sampled PROMICE pixels show daily `SU` differs from summed `SUH` because the
  former also includes evaporation.
- `RZ` is daily "Meltwater refreezing and deposition" (`mmWE/day`) on the
  singleton permanent-ice sector. It maps to combined `refreeze_deposition`, not
  to pure RACMO `refreeze`, using the same daily/24 previous-hold support.
- `MEH` remains the public hourly `melt`. Daily `ME` is read from the singleton
  permanent-ice sector only as an independent daily-sum validation and is not
  emitted as a duplicate melt channel. Sampled PROMICE pixels agree to numerical
  precision.

Artifacts record exact source variables, sectors, signs, daily support, the
`ME`/`MEH` validation ledger, and a conservative symmetric `0.05 mWE/h`
absolute sanity bound that catches unconverted `mmWE` values. `RZ` is preserved
as a native signed combined term rather than clipped or relabelled as pure
refreezing. Across the 43 selected artifacts, 44 strict-negative UTC days occur
at 15 sites; 36 are tiny source-exact values between about `1e-11` and
`1.6e-9 mWE/h`. The source investigation's reporting threshold of `1e-8 mWE/h`
separately identifies eight material negative days at three sites in both
deposition-sign and active-melt regimes. Metadata records the canonical signed
policy, strict-negative day count/minimum, and reporting-only material threshold
plus material day count/minimum. That threshold never clips values and is not an
acceptance tolerance. `auditArtifacts` recomputes both statistic sets while
retaining the broad unit and support gates. A reduced or mixed-year source records
`mar_diagnostic_status=not_available` or `partial`; this does not change the
essential model-forcing readiness contract.

MAR also carries fixed-depth density `RO1(OUTLAY)` and dynamic-layer
`ROSN1(DZSN1)`. Both are scientifically useful, but neither is a scalar time
series and the committed reduced fixture omits them. They are therefore deferred
to a separate profile-bundle staging task rather than flattened into userdata or
manufactured as long time-series panels.

## MAR bounded daily diagnostics and snow depth

Daily MAR state diagnostics (`SHSN2`, `CC`, `ST`, and `SP`) use linear
interpolation only inside their native support. Before the first and after the
last daily sample of a separately processed year they hold the nearest endpoint,
rather than extrapolating the last slope. Cloud fraction additionally enforces
the source and output range `[0,1]`, preventing the former December 31 overshoot.

MAR `snowd` remains `SHSN2`, whose source definition is "Snow Pack Height above
Ice". It is never replaced by `SHSN3`, which is total multilayer snow/firn
thickness. An archive-calibrated annual-boundary screen identifies an interior
pixel-year only when severe jumps bracket it and its annual median is separated
from both neighbouring years. That isolated year is masked in the artifact;
one-sided or ambiguous boundary years remain finite but are explicitly marked
unverified. The calibration and affected years are saved in artifact metadata
and checked by `icemodel.verification.auditArtifacts`.

## RACMO point selection and ice mask

`buildRacmoData` applies the companion FGRN11 topography mask to point sampling
and conservative polygon extraction. A point uses the nearest cell whose fractional `IceMask`
rounds to ice (fraction at least 0.5; `Promicemask > 0` when `IceMask` is
unavailable). The valid cell must also lie within one native grid diagonal of
the requested point. Otherwise the point is explicitly unavailable: verification
staging does not substitute a distant inland cell or preserve an older off-mask
artifact. Natural-neighbour point interpolation excludes masked cells from its
local neighbourhood; conservative polygon remapping uses the same mask. Point
sampling fails closed when the companion topography/mask file is absent. A
legacy point artifact without `racmo_ice_mask_applied` and
`racmo_point_max_distance_m` provenance is not cache-compatible; one canonical
restage is required because metadata repair cannot change its sampled payload.

## RACMO precipitation numerical undershoot

RACMO 2.3p3 labels `precip` as a 3-hour-mean precipitation flux in
`kg m-2 s-1`. The native model field contains small finite negative numerical
undershoots at otherwise valid ice cells. Because precipitation is physically
nonnegative, `buildRacmoData` enforces zero as a source invariant after point or
polygon sampling/remapping and optional hourly interpolation, once `ppt` is in
the canonical `m s-1` unit. There is no magnitude threshold: every finite
negative becomes exactly zero, while NaN, zero, legitimate positive values,
time, and unrelated channels are unchanged.

Artifacts record the source variable, normalization stage, input minimum, and
replacement count in `racmo_ppt_qc_*` metadata. The exact-reference repair path
uses the same helper, so a second dry run is unchanged and source-light rather
than silently relaxing the generic nonnegative precipitation QA bound.

## GEUS MODIS albedo coverage provenance

MAR, MERRA-2, and RACMO builders can add the optional daily GEUS Greenland
Reflectivity 5 km C6 albedo as `modis`. Fresh artifacts record the same exact
coverage contract used by metadata repair and artifact QA:

- `modis_product = "GEUS Greenland Reflectivity 5km C6"`;
- `modis_status = "source_coverage"` when every artifact year is covered,
  `"partial_source_coverage"` when covered and missing years coexist, or
  `"no_source_coverage"` when none is covered;
- `modis_coverage_years` is the sorted source-backed subset of artifact years.

Missing years are never filled. Partial artifacts retain `modis` with NaNs in
uncovered years. A no-source artifact carries the explicit metadata but omits the
physical channel, matching the idempotent repair representation. A file-matched
covered year must contribute at least one finite real albedo value in `[0,1]` on
the artifact axis; otherwise staging fails instead of claiming false coverage.
MAR/MERRA Data and their derived met share one physical read and provenance
record, and the writers save identical payload `Properties.UserData` and
top-level `artifact_metadata`.

## PROMICE / GC-Net AWS evaluation data

`readPromiceAws` reads a pypromice Level-3 station NetCDF into canonical
icemodel channels; `buildPromiceMet` produces a met-contract forcing file;
`buildPromiceData` produces the observational evaluation/userdata timetable.
Everything below is governed by the product readme
(`data/verification/promice/AWS_data_readme.pdf`) and the variable dictionary
(`AWS_variables.csv`). Read those before changing the builders.

### Shortwave source and placeholder policy

`readPromiceAws` keeps the pypromice source channels distinct: raw pyranometer
`dsr`/`usr` are exposed as `swd`/`swu`, and the product's tilt/bias-corrected
`dsr_cor`/`usr_cor` are exposed as `swd_cor`/`swu_cor`. The PROMICE variable
dictionary permits raw `dsr` down to -10 W m-2, so small negative nighttime
sensor offsets are source-valid measurements but are not physical public
radiative fluxes.

Both `buildPromiceMet` and `buildPromiceData` therefore use the corrected
product where finite, fall back to the raw channel where the correction is
missing, and clamp any remaining finite negative selected value to zero. The
reader remains the source-faithful inspection boundary. Saved metadata names
the raw/corrected channels and records finite corrected use, raw fallback, raw
negative, corrected-negative, raw-fallback-negative, and final clamp counts plus
the raw minimum.

The builders also resolve a narrow source-level ambiguity needed by usable
forcing: a missing hourly `swd` or `swu` sample is set to physical zero only
when the complete source interval `[t,t+1h)` stays at or below -6 degrees solar
elevation. The NOAA fractional-year approximation is evaluated in UTC at the
interval start, quarter-hours, and end using the station coordinates. Missing
twilight/daylight samples remain `NaN`; the builders do not interpolate them or
copy the broad nonzero fills found in older demo artifacts. Saved metadata
records the method, threshold, derived-zero count, and remaining-missing count
for each channel. This rule is applied identically to hourly observational Data
and to the hourly source that is expanded onto the default 15-minute met grid.
Source-backed availability is measured on the complete station product before
window selection, so full and surgical builds return identical values for the
same timestamps; a channel absent or all-missing across that complete product
remains an intentional placeholder.

If shortwave or albedo is absent or entirely missing, the met schema retains an
all-NaN channel for a later forcing swap and records a channel-specific
`swd_policy` or `albedo_policy` containing an explicit `NaN placeholder`
reason. No observation or constant radiation is invented. Observational Data
omits a source-absent channel and carries the same provenance metadata; the
`observations.mat` bundle reuses that exact Data timetable.

### Precipitation policy

PROMICE precipitation availability is station- and period-specific. In the
current staged hourly product, 35 of 51 station files carry `rainfall_cor_u`;
33 contain positive samples, while 16 omit the channel. The variable dictionary
defines it as corrected liquid rainfall within the timestep [mm], derived from
the tipping-bucket gauge; the gauge is not a reliable solid-precipitation
measurement. `buildPromiceMet` therefore converts available hourly amounts to
canonical `rainf` [m s-1], but never treats liquid-only rainfall as total
precipitation. `snowf` and `ppt` remain all-NaN placeholders for runtime
fill/source swapping, and `rainf` is also a placeholder when the station lacks
the source channel. Native PROMICE artifacts are not runnable on precipitation
alone until total precipitation is supplied by that later step.

### Time convention (the one canonical rule)

**Canonical convention (single source of truth): a forcing/eval hourly value
represents the interval `[t, t+dt)` and is LABELLED AT THE INTERVAL-START `t`.**
This matches PROMICE ("the timestamp of the hourly averages indicate the start
of the averaged hour") and the model's own `[t, t+dt)` implicit integration. The
convention is applied **at the builder, where the source stamping is known** —
it is NOT special-cased downstream.

The L3 **hourly timestamp is the START of the averaged hour** (readme,
"Temporal averaging"). The NetCDF encodes it as integer hours since a station
epoch, so `epoch + hours(t)` reproduces the bin-START stamp exactly;
`readPromiceAws` snaps with `dateshift('start','hour')` defensively (idempotent,
not a re-bin) and returns a UTC axis. `buildPromiceMet` / `buildPromiceData`
inherit this axis unchanged. icemodel's met/Data axis is this **same bin-START
hourly grid (UTC)**: the timestepping loop treats a met row's `Time` as the
forcing valid AT that timestamp and integrates forward over `[t, t+dt)`, so a
START-of-hour averaged forcing is the correct mean to drive that step and **no
half-hour recentring is applied or needed**.

#### Comparison protocol (cumulative/flux vs instantaneous state)

Align simulation to observation **on the START-of-hour stamp** — do not recentre
to the hour middle. Two cases, by quantity type:

- **Cumulative / flux quantities** (ablation, accumulation, SMB, energy fluxes,
  any per-interval mean): the model's per-interval output is itself the mean/
  total over `[t, t+dt)`, so it aligns **exactly** to the same interval label as
  the observation. No offset.
- **Instantaneous STATE** (e.g. temperature): the model uses implicit
  backward-Euler, so the solved state is the **end-of-interval** value (valid at
  `t+dt`) while it is stored against the interval-start label `t`. This is a
  sub-`dt` labelling offset: at hourly resolution it is **negligible** (one
  hour against a multi-day-to-seasonal subsurface thermal signal) and is **not
  corrected**. Document, do not special-case.

#### Other forcing sources (native stamping vs this convention)

The gridded met builders carry their native stamping; only PROMICE is snapped:

- **MAR** (`readMar3p11`): yearly files start `Jan 1 00:00 UTC`; the hourly axis
  is `t0 + (0:23)h` per day — **interval-start, matching this convention**.
- **MERRA-2** (`readMerra2` / `buildMerraData`): the reader preserves native
  snapshot or interval-center coordinates exactly. The application builder
  consumes only time-averaged collections: `slv`/`rad`/`flx` are tavg1 centered
  at `:30`, while `glc` is tavg3 centered at `01:30`, `04:30`, and so on. The
  builder decodes and validates every daily native coordinate against its
  filename day, explicitly subtracts half of each declared support interval,
  then holds each mean over its one- or three-hour `[t,t+dt)` support. It never
  shifts an instantaneous collection or linearly bridges missing source
  intervals.
  `hasCanonicalMerraTimeSupport` checks the declared time-support provenance,
  `hasProvenMerraTavg3SourceGrid` records the exact native glc timestamp
  inventory (including any omitted 3-hour stamps),
  while `hasConstantMerraTavg3Support` independently requires the staged
  `runoff`, `albedo`, `snowd`, and `swe` values to remain constant within every
  available three-hour support block. All three checks are required: clock
  alignment and policy metadata alone cannot prove that a legacy regularizer did
  not invent a value at an omitted native source stamp.

### Ablation vs accumulation channel semantics (the core rule)

`buildPromiceData` branches on the **presence of `z_ice_surf`** in the L3 file
(the operational ablation-site signal; recorded in
`metadata.site_surface_type`, and agreeing with the readme Table 1 "Site type"
surfaced via `promiceSiteCatalog`):

- **Ablation sites** (z_ice_surf present): two surface channels —
  - `ablation` [m, +down] = `-(z_ice_surf - z_ice_surf(window start))`,
    cumulative ice-surface lowering relative to installation.
  - `snow_depth` [m] = `snow_height` clamped `>= 0` (the readme calls
    snow_height "strictly positive"; negatives are counted in
    `metadata.snow_depth_negatives` then clamped, never deleted).
- **Accumulation / percolation / bedrock sites** (no z_ice_surf): one channel —
  - `surface_height` [m, **+up**] = `z_surf_combined - z_surf_combined(start)`,
    a NET surface-height change. This is **not** an ablation channel and is not
    relabeled as one. No `snow_depth` is fabricated (the L3 `snow_height` at
    these sites is not a true snow-over-ice depth).

The branch keeps each channel physically meaningful: accumulation sites emit
only the net `surface_height` they actually measure, and no ablation/snow_depth
channel is synthesized where the underlying `z_ice_surf`/snow-over-ice depth does
not exist (e.g. EGP, which carries only `surface_height`).

### Surface heights are relative to installation

**All L3 surface heights are relative to (zero at) the surface height at the
initial station installation** (readme). Before comparing to a simulation,
**subtract the window-start value** of the surface-height channel (the builder
already zeroes `ablation`/`surface_height` at the first finite window sample).

### What is GEUS-provided vs ours (we flag, never silently fix)

The authoritative GEUS L3 product provides **all QC, manual flagging/fixing, and
gap-filling** (slope-bridging the surface height across sensor outages), the
multi-sensor `z_surf_combined`, the `z_ice_surf` re-derivation, `snow_height`,
the depth-tagged thermistor string, and the standardized 10 m `t_i_10m`. **We
provide only**: read the L3 channels; derive per `site_type`
(`ablation`/`surface_height`); clamp `snow_depth >= 0`; unit conversions; and
per-sample FLAGS. **No GEUS data value is modified** — we flag, never silently
fix. The one correction (de-stepping) is an opt-in transform applied at analysis
time, never baked into the staged data (see below).

### Gap flag (sensor-derived; trend usable, exclude only RATE diagnostics)

Surface-height change across data gaps is **not a direct observation**: the
readme bridges gaps by a manual slope when **all** surface-ranging sensors fail.
But the readme also states the **trend over the whole period is unaffected by the
gaps** — so the cumulative series stays usable; only per-timestep RATES through a
gap are unreliable.

`buildPromiceData` derives `surface_height_flag` from FIRST PRINCIPLES via
`icemodel.forcing.helpers.surfaceFlags`: a sample is gap-bridged (`1`) when the
surface series is finite **but every underlying L3 surface-ranging sensor**
(`transducer_depth`, `boom_height`, `stake_height` as available) **is NaN** —
i.e. the value is slope-interpolated, not measured. Leading/trailing samples
(before first / after last finite surface value) are flagged too. The old
heuristic flagged only samples where the surface value itself was NaN, which
**missed the slope-bridged segments** (the surface is finite there, manufactured
by interpolation) — e.g. MIT: ~6.8k slope-bridged samples the old flag missed.

**Comparison guidance (revised):**

- **CUMULATIVE and visual comparison use the FULL series** (the trend is
  preserved through gaps; do NOT drop gap-bridged samples from a cumulative or
  visual comparison).
- **Only RATE-based diagnostics** (per-timestep rate) flag/exclude
  `surface_height_flag == 1` segments.

Daily retime aggregates the flag by `max`. Data are FLAGGED, never deleted.
`metadata.gap_flagged_samples` records the count. (The daily `max` can spread a
flag to a whole day touched by one bridged hour; the sensor-derived flag now
fires on far fewer, genuinely-bridged samples, so the figure markers no longer
overlay observed samples.)

### Station-transition flag and de-stepping (opt-in correction)

A PROMICE **site** can merge several **stations** (AWS); at a handover the
surface or subsurface series can carry an expected discrete offset — a **step**,
not a NaN, so the gap flag never sees it. `buildPromiceData` stages two further
flag families:

- `station_transition_flag` marks samples within a station-handover window. The
  staged product carries only a SITE-level install date (not per-station
  handover dates), so this flag is currently inert; the FACT that a site merges
  multiple AWS is recorded in `metadata.is_multistation` /
  `metadata.composing_stations` (from `promiceSiteCatalog.stations`).
- `step_detected_flag`, `step_correctable_flag`, `step_magnitude` stage the
  de-stepping DETECTION run by `icemodel.forcing.destepSurface` in `detect`
  mode. Detection scores each single-timestep jump (between adjacent finite
  samples; a gap interior is interpolation, not a step) on multiple evidence
  lines — rate-implausible magnitude (gate), gross single-step implausibility,
  station-transition coincidence, and melt-season inconsistency — and classifies
  it **UNAMBIGUOUS** (magnitude + a corroborating line) or **AMBIGUOUS**
  (magnitude alone). Ambiguous steps are flagged, never corrected.

**The staged `.mat` is faithful** (raw surface values + the flags). The
de-stepping CORRECTION is applied **opt-in at analysis time** by calling
`icemodel.forcing.destepSurface(t, surf)`; its **default mode corrects
UNAMBIGUOUS steps only** (so unambiguous installation jumps — e.g. MIT's two
~5.9 m jumps in Aug 2009 — are leveled by default at consumption, while the
staged data stays unaltered). See `test/interactive/promice_qaqc/figures/`
`promice_step_screening.md` for the all-site screening table.

### Thermistors (subsurface temperature) and the tice10m comparison protocol

- `tice10m` [K] is the **PRIMARY** subsurface-temperature evaluation channel.
  It is GEUS's **standardized 10 m-BELOW-the-EVOLVING-SURFACE** temperature: GEUS
  builds it by tracking each thermistor's **time-dependent depth below the
  current surface** (`d_t_i_*`, which changes as the surface ablates or
  accumulates), discarding surfaced thermistors, and depth-interpolating the
  surviving subsurface string to 10 m below the **current** surface at each time
  step (`t_i_10m`). It is therefore a **moving (Lagrangian) 10 m depth**, not a
  fixed 10 m from installation.

  **Model sampling protocol:** the model must be sampled at **10 m below its OWN
  current surface** — a moving Lagrangian depth at ablation sites, where the
  surface lowers over time — **NOT** at a fixed 10 m from the start of the run.
  Sampling a fixed installation-referenced 10 m at an ablation site compares
  different physical material as the surface drops and is wrong. Track the
  model's current surface and sample 10 m beneath it to match `tice10m`.

- `tice1..ticeN` [K] is the depth-tagged raw string (`t_i_1..N`), each with a
  companion `dtice1..dticeN` [m] depth (`d_t_i_*`, positive = below surface).
  `readPromiceAws` clamps each to the dictionary range `[-80, 1] C` and
  **discards surfaced thermistors** (`d_t_i_N <= 0`, sensor at/above the
  surface as ablation-site firn melts out) per the readme. Use the depth tag to
  place each sensor vertically; a NaN tice sample is either out-of-range,
  surfaced, or simply missing. The raw string is **secondary / diagnostic**: to
  compare it the model must be sampled at each sensor's time-dependent depth
  `d_t_i_N`, so prefer `tice10m` as the standardized primary channel.

### Specific humidity (shum)

L3 `qh_u` is specific humidity reported in **g/kg** despite the product's
"kg/kg" label (both the dictionary and the NetCDF units attribute are wrong;
KAN_M median 1.68). `readPromiceAws` divides by 1000 so `shum` is **kg/kg**,
matching the MAR/MERRA convention the vapor kernel
(`icemodel.vapor.relative_humidity_from_specific_humidity`) consumes. (PROMICE
forcing uses its own measured `rh` channel; the conversion is for cross-source
consistency.)

### Station vs site (expected step-shifts)

A **site** can merge several **stations** (AWS) over time. `promiceSiteCatalog`
carries the composing `stations` per site and the `site_type`. A step-shift in
the surface or subsurface series that coincides with a station transition is
**expected**, not a defect (e.g. CEN ~2021 CEN2/CEN1/GITS, QAS_U QAS_Uv3/QAS_U).
The `station_transition_flag` records the handover window when per-station dates
are available, and the de-stepping detector uses transition coincidence as one
evidence line. With only a SITE-level install date staged, the actual handover
TIME is recovered from the data (a coincident step), and a step that cannot be
corroborated stays AMBIGUOUS (flagged, not corrected).

### Diagnostics

`test/interactive/promice_qaqc/plot_promice_eval.m` regenerates per-site figures and the QA
summary into the gitignored `test/interactive/promice_qaqc/figures/` (site-type surface
channel; surface-height with gap-bridged samples in red, station-transition
windows in green, detected steps marked unambiguous/ambiguous; `tice10m`
heavy-black primary over the thermistor string). `promice_site_cause_classification.md`
there is the manual cause classification, and `promice_step_screening.md` is the
all-site step-screening table.
