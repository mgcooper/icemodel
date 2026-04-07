# icemodel.vapor

Purpose: vapor thermodynamics and atmosphere-surface vapor conversions shared across domains.

Public entry points:
- `icemodel.vapor.vappress`
- `icemodel.vapor.vappress2rh`
- `icemodel.vapor.vapordensity`
- `icemodel.vapor.vapordiffusivity`
- `icemodel.vapor.vapork`
- `icemodel.vapor.tdewpoint`
- `icemodel.vapor.twetbulb`
- `icemodel.vapor.vaporinit`
- `icemodel.vapor.moist_air_density`
- `icemodel.vapor.specific_humidity_from_vapor_pressure`

Contents:
- saturation vapor pressure relations
- RH/dew-point/wet-bulb conversions
- vapor density, diffusivity, and conductivity helpers

Rules:
- keep general thermodynamic transforms here
- keep column microstructure transport, such as `VAPORTRANSFER`, in `icemodel.column`

Migration status: active vapor thermodynamics and atmosphere-surface conversion helpers now live here.
