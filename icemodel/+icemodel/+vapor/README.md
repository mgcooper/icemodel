# icemodel.vapor

Purpose: vapor thermodynamics and atmosphere-surface vapor conversions that
every domain shares.

Public entry points:
- `icemodel.vapor.saturation_vapor_pressure`
- `icemodel.vapor.relative_humidity_from_vapor_pressure`
- `icemodel.vapor.saturation_vapor_density`
- `icemodel.vapor.vapor_diffusivity`
- `icemodel.vapor.vapor_thermal_conductivity`
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
- keep column transport and microstructure updates, such as
  `vapor_transport_faces` and `update_grain_radius`, in
  `icemodel.column`
