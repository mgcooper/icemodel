# icemodel.forcing

General-purpose builders that turn raw meteorological source data
(PROMICE weather stations; MAR, RACMO, and MERRA2 gridded products) into
the standard icemodel forcing and evaluation artifacts.

## Per-source builders

| Builder | Source | Produces |
| ------- | ------ | -------- |
| `buildPromiceMet` / `buildPromiceData` | PROMICE v3 AWS NetCDF (any station, any year; reader `readPromiceAws`) | met file; evaluation Data (SEB channels, cumulative ablation, snow depth, ice-temperature string) |
| `buildMarMet` / `buildMarData` | MAR v3.11 yearly NetCDF (reader `readMar3p11`, grid `marGridInfo`) | met file; forcing/userdata Data (point or polygon) |
| `buildRacmoData` | RACMO2.3p3 per-variable 3-hourly archive | evaluation Data (SMB components; RACMO carries no met channels) |
| `buildMerraData` | MERRA-2 daily collection files (slv/rad/flx/glc) | forcing/userdata Data; convert with `data2met` |

Locations: PROMICE takes a station name; gridded sources take `[lat lon]`
(point — `method="nearest"` default, or `"natural"` natural-neighbour) or
a `polyshape` in EPSG:3413 meters (catchment). For polygons,
`remap="conservative"` (default) does exact overlap-area-weighted remapping
via the `exactremap` toolbox — with the source ice mask inpainted and, for
RACMO, the true FGRN11 cell areas; `remap="equal"` is a dependency-free plain cell mean.
Windows/years are arbitrary within the source archive; calendars derive
from the files present.

Dependency: conservative polygon remap requires the `exactremap` toolbox on
the path (`https://github.com/mgcooper/exactremap`); point builds and
`remap="equal"` do not. The dependency-resolution pass
(`icemodel.internal.installRequiredFiles`) will add it permanently.

## Conservative polygon remap (per source)

For polygon (catchment) locations, `remap="conservative"` computes the exact
overlap-area-weighted mean in each source's **native grid frame** — never by
reprojecting the native grid to a common projection, which would render it
curvilinear. Each source supplies the geometry needed:

- **MERRA-2** — regular geographic lon/lat grid; remapped directly in lon/lat
  (`UseGeoCoords=true`). A static land-ice mask derived from glacier-tile (glc)
  validity restricts the average to ice.
- **RACMO2.3 FGRN11** — regular in rotated-pole `rlon`/`rlat`; remapped in the
  rotated frame using the CF `rotated_pole` parameters, the shipped true cell
  areas (`gridarea`), and the `IceMask`, all from the FGRN11 topography file.
- **MAR v3.11** — regular in native projected km (`X`/`Y`); the catchment
  polygon is mapped from EPSG:3413 into native km via the shipped 2-D `LON`/
  `LAT` (no MAR projection formula is required), then clipped against the
  regular native grid with the `SRF==4` ice mask. Cells are uniform 225 km².

Source no-data is screened before remap (MAR `_FillValue` ~1e36 and the 999
mass-flux sentinel; MERRA ~1e15 fill), and masked/missing cells are inpainted
so partial-coverage catchments are not biased. EPSG:3413 is the default
geometry for the Greenland workflow; true cell areas are used where shipped, so
catchment totals are unaffected by projection area distortion. The
dependency-free `remap="equal"` is a plain in-polygon cell mean — adequate when
many native cells fall inside, weaker for sub-few-cell catchments.

For point locations, `method="nearest"` returns the native cell value (no
blending across surface-type boundaries); `method="natural"` is a deterministic
natural-neighbour interpolation over a fixed symmetric 5×5 cell pad, giving a
continuous, horizontally elevation-aware estimate where a site sits between
cells of differing model elevation. Natural-neighbour is horizontal blending
only; it does not lapse-correct for vertical elevation differences — forcing is
stored source-faithful at the model cell elevation.

## Source observation heights and channels

- **MAR** near-surface hourly fields are at standard heights: air temperature
  and humidity (`TTH`, `QQH`) at 2 m, wind (`UUH`, `VVH`) at 10 m. Humidity is
  derived from `QQH` + surface pressure + `TTH`; MAR's multi-level products
  (sigma, pressure, and height levels) are not used. `icemodel.setopts` sets
  `opts.z_tair=2`, `opts.z_wind=10` for `forcings="mar"`, so the turbulent-flux
  scheme receives the correct measurement heights.
- **MERRA-2** uses 2 m temperature/humidity and 2 m winds (`U2M`, `V2M`);
  `opts.z_tair=opts.z_wind=2` for `forcings="merra"`. `swe` is kept in native
  `kg m-2` (snow mass), distinct from `snowd` (snow depth, m).
- **RACMO** ships no albedo variable, so surface albedo is recovered from the
  shipped shortwave fluxes (`1 - swn/swd`) only where `swd ≥ 10 W m-2`, leaving
  low-sun timesteps to the QA/QC gap fill; where filled, `swd ≈ 0` so the
  downstream `swn = swd*(1-albedo)` reconstruction is insensitive to the value.

## Met-file contract

Required variables (see `helpers.metvariables`): `tair` [K], `swd`
[W m-2], `lwd` [W m-2], `albedo` [-], `wspd` [m s-1], `rh` [%], `psfc`
[Pa], `ppt` [m water equivalent]. Optional pass-through diagnostics:
`rainf`, `snowf`, `snow_depth`, `melt`, `runoff`, `shf`, `lhf`, `tsfc`,
`cfrac`, `snowd`, `wdir`. `helpers.validatemet` enforces the contract at
the write boundary.

## QA/QC (`helpers.metchecks`)

Every builder runs `metchecks` before writing: per-variable NaN/complex
counts (returned for provenance), linear gap fill with nearest-value end
fill, and physical-range clamps:

| variable | clamp |
| -------- | ----- |
| albedo   | [0.05, 0.98] |
| rh       | [5, 99.99] % |
| wspd     | >= 0.1 m s-1 |
| wdir     | wrapped to (0, 360] |
| tsfc     | <= 273.16 K (or <= 0 C) |

Wind direction is a circular variable and is gap-filled through its
unit-vector (sin/cos) components, so fills that cross the 360/0 wrap are
correct.

## Method notes

- PROMICE rh comes from the v3 `rh_wrtwater` channel (relative humidity with
  respect to water, matching the `icemodel.vapor` convention; already percent
  despite its units attribute). The plain `rh` channel is relative humidity
  with respect to ice (legitimately exceeds 100% when subfreezing) and is not
  used.
- PROMICE ablation is derived from the transducer-depth record: curated KAN
  service windows plus generic removal of positive depth jumps (re-drills).
  It is first-order; treat cumulative totals as approximate.
- PROMICE snow depth is a September-median boom-height reference minus boom
  height, clamped at zero (first-order).
- MERRA swd is read directly from `SWGDN` (downwelling shortwave).
- Precipitation-rate units are currently source-conventional (mWE/h for
  MAR/MERRA, m/s in ESM-SnowMIP staged forcing) and not yet harmonized.

## Artifact types

- **met files** — model forcing. A timetable named `met` saved as
  `met_<site>_<forcings>_<YYYYMMDD>_<YYYYMMDD>_<dt>.mat` (window form)
  or `met_<site>_<forcings>_<YYYY>_<dt>.mat` (per-year form). Loaded by
  `icemodel.loadmet` via `icemodel.createMetFileNames`.
- **Data files (userdata / met-swap)** — a timetable named `Data` saved
  as `<site>_<source>_<YYYY>.mat` under `<input>/userdata/`. When
  `opts.userdata` / `opts.uservars` are set, `icemodel.loadmet` swaps
  the named met-file variables with the Data file's columns. Data
  timetables carry location metadata as table CustomProperties:
  X, Y, Lat, Lon, Elev, Slope, ScalarUnits.
- **eval files** — observational evaluation data (e.g. PROMICE ablation)
  under the eval root resolved by `icemodel.getpath('eval')`.

## Conventions

- Humidity is computed with the `icemodel.vapor.*` kernels
  (`relative_humidity_from_specific_humidity`).
- Outputs are double precision; channels are not rounded.
- Projection: WGS 84 / NSIDC Polar Stereographic North (EPSG:3413) via
  `helpers.psnProjection`.
- Paths resolve through `icemodel.getpath` / `icemodel.config`; builders
  accept explicit `source_dir`-style arguments for raw inputs. No
  environment variables are introduced by this namespace.
- Locations: station/site names for PROMICE; projected points or a
  catchment polygon for gridded sources (polygon averaging produces
  catchment-mean Data files).

## Raw source data

Raw inputs are not committed. Each builder takes a `source_dir` argument
pointing at the corresponding raw archive (MAR v3.11 NetCDF; PROMICE v3 AWS
NetCDF; GEUS AWS v3; RACMO2.3; MERRA-2; GEUS MODIS albedo). An optional
gitignored cache under `data/forcing/<source>/` may mirror subsets of these,
following the manual-fetch pattern of
`icemodel.verification.setup.fetchEsmSnowmip`.
