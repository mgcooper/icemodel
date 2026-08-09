# icemodel.vapor

Purpose: vapor thermodynamics and atmosphere-surface vapor conversions shared across domains.

Public entry points:
- `icemodel.vapor.saturation_vapor_pressure`
- `icemodel.vapor.relative_humidity_from_vapor_pressure`
- `icemodel.vapor.saturation_vapor_density`
- `icemodel.vapor.vapor_diffusivity`
- `icemodel.vapor.vapor_thermal_diffusion_coefficient`
- `icemodel.vapor.dew_point_temperature`
- `icemodel.vapor.wet_bulb_temperature`
- `icemodel.vapor.initialize_vapor_model`
- `icemodel.vapor.moist_air_density`
- `icemodel.vapor.specific_humidity_from_vapor_pressure`

Contents:
- saturation vapor pressure relations
- RH/dew-point/wet-bulb conversions
- vapor density, diffusivity, and conductivity helpers

Rules:
- keep general thermodynamic transforms here
- keep column microstructure transport, such as `vapor_mass_transfer`, in `icemodel.column`

Migration status: active vapor thermodynamics and atmosphere-surface conversion helpers now live here.
