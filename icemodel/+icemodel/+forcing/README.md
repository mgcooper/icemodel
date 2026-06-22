# icemodel.forcing

Forcing and evaluation-data builders for icemodel. Each gridded/AWS source has
a `read*` (raw -> canonical channels), a `build*Met` (met-contract timetable),
and a `build*Data` (evaluation/userdata timetable) function.

## PROMICE / GC-Net AWS evaluation data

`readPromiceAws` reads a pypromice Level-3 station NetCDF into canonical
icemodel channels; `buildPromiceMet` produces a met-contract forcing file;
`buildPromiceData` produces the observational evaluation/userdata timetable.
Everything below is governed by the product readme
(`data/verification/promice/AWS_data_readme.pdf`) and the variable dictionary
(`AWS_variables.csv`). Read those before changing the builders.

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
- **MERRA-2** (`readMerra2`): `tavg1` hourly bins are stamped at the **bin
  CENTER** (`:30`), kept native by `merraTime`. This is a half-hour offset from
  the interval-start convention; it is small at hourly resolution and currently
  left native (noted, not silently shifted).

### Ablation vs accumulation channel semantics (the core rule)

`buildPromiceData` branches on the **presence of `z_ice_surf`** in the L3 file
(the operational ablation-site signal; recorded in
`metadata.site_surface_type`, and agreeing with the readme Table 1 "Site type"
surfaced via `promicesiteinfo`):

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
  `metadata.composing_stations` (from `promicesiteinfo.stations`).
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
staged data stays unaltered). See `test/interactive/figures/`
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

A **site** can merge several **stations** (AWS) over time. `promicesiteinfo`
carries the composing `stations` per site and the `site_type`. A step-shift in
the surface or subsurface series that coincides with a station transition is
**expected**, not a defect (e.g. CEN ~2021 CEN2/CEN1/GITS, QAS_U QAS_Uv3/QAS_U).
The `station_transition_flag` records the handover window when per-station dates
are available, and the de-stepping detector uses transition coincidence as one
evidence line. With only a SITE-level install date staged, the actual handover
TIME is recovered from the data (a coincident step), and a step that cannot be
corroborated stays AMBIGUOUS (flagged, not corrected).

### Diagnostics

`test/interactive/plot_promice_eval.m` regenerates per-site figures and the QA
summary into the gitignored `test/interactive/figures/` (site-type surface
channel; surface-height with gap-bridged samples in red, station-transition
windows in green, detected steps marked unambiguous/ambiguous; `tice10m`
heavy-black primary over the thermistor string). `promice_site_cause_classification.md`
there is the manual cause classification, and `promice_step_screening.md` is the
all-site step-screening table.
