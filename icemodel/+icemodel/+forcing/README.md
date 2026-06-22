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

### Time convention (model-vs-observation alignment)

The L3 **hourly timestamp is the START of the averaged hour** (readme,
"Temporal averaging"). The NetCDF encodes it as integer hours since a station
epoch, so `epoch + hours(t)` reproduces the bin-START stamp exactly;
`readPromiceAws` snaps with `dateshift('start','hour')` defensively (idempotent,
not a re-bin) and returns a UTC axis.

icemodel's met/Data axis is this **same bin-START hourly grid (UTC)**. The
timestepping loop treats a met row's `Time` as the forcing valid AT that
timestamp and integrates forward over `[t, t+dt)`; a START-of-hour averaged
forcing is the correct mean to drive that step, so **no half-hour recentring is
applied or needed**. When comparing a simulation to these observations, align
on the START-of-hour stamp — do not recentre to the hour middle.

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

This is a correctness fix: the previous builder fabricated an `ablation` and
`snow_depth` for every site (falling back to z_surf_combined for accumulation),
which produced mirror-image spurious "ablation" spikes (e.g. EGP).

### Surface heights are relative to installation

**All L3 surface heights are relative to (zero at) the surface height at the
initial station installation** (readme). Before comparing to a simulation,
**subtract the window-start value** of the surface-height channel (the builder
already zeroes `ablation`/`surface_height` at the first finite window sample).

### Gap flag (exclude gap-bridged segments)

Surface-height change across data gaps is **not a direct observation** (the
readme bridges gaps by a manual slope). `buildPromiceData` attaches a per-sample
`surface_height_flag` to whichever surface channel the site carries:

- `0` = direct observation (underlying L3 sensor finite here),
- `1` = gap-bridged / interpolated (sensor missing here).

Daily retime aggregates the flag by `max` (a day touched by any gap-bridged hour
stays flagged). **Comparison code must exclude `surface_height_flag == 1`
samples** before scoring the surface-height channel. Data are FLAGGED, never
deleted. `metadata.gap_flagged_samples` records the count.

### Thermistors (subsurface temperature)

- `tice10m` [K] is the **PRIMARY** subsurface-temperature evaluation channel:
  the standardized 10 m firn temperature GEUS depth-interpolates (`t_i_10m`).
  Prefer it for site-to-site model comparison.
- `tice1..ticeN` [K] is the depth-tagged raw string (`t_i_1..N`), each with a
  companion `dtice1..dticeN` [m] depth (`d_t_i_*`, positive = below surface).
  `readPromiceAws` clamps each to the dictionary range `[-80, 1] C` and
  **discards surfaced thermistors** (`d_t_i_N <= 0`, sensor at/above the
  surface as ablation-site firn melts out) per the readme. Use the depth tag to
  place each sensor vertically; a NaN tice sample is either out-of-range,
  surfaced, or simply missing.

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

### Diagnostics

`test/interactive/plot_promice_eval.m` regenerates per-site figures and the QA
summary into the gitignored `test/interactive/figures/` (site-type surface
channel, surface-height with gap-bridged samples in red, tice10m + thermistor
string). `promice_site_cause_classification.md` there is the manual cause
classification for the user-flagged sites.
