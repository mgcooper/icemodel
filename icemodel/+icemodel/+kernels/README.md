# icemodel.kernels

Purpose: Self-contained physical-property formulas. Some are production
kernels. The others are reference formulas that plots and docstrings cite.

Contents:

- Production kernels:
  - `air_kinematic_viscosity` approximates the kinematic viscosity of air.
    `icemodel.surface.initialize_surface_state` and the bulk Richardson
    diagnostics call it.
  - `potential_surface_vapor_demand` converts latent heat to a surface vapor
    demand. `icemodel.surface.potential_surface_vapor_demand` calls it.
- Reference vapor formulas:
  - `buckVaporModel` holds the Buck (1981) vapor model.
    `icemodel.plot.vaporModel` calls it, and
    `icemodel.vapor.saturation_vapor_pressure` cites it. It computes:
    - saturation vapor pressure
    - vapor density
    - vapor thermal conductivity
    - dew point
  - `saturationVaporPressure` holds the Romps and Ambaum saturation vapor
    pressure. `icemodel.plot.vaporModel` calls it, and
    `icemodel.vapor.initialize_vapor_model` cites it.
  - `latentEnthalpyWater` holds the Romps and Ambaum latent enthalpy.
    `icemodel.vapor.initialize_vapor_model` cites it.
- Conductivity formulas, each plotted by the `icemodel.plot` function with
  the same name:
  - `thermal_conductivity_air` returns the conductivity and its derivative.
  - `thermal_conductivity_water` returns the conductivity and its derivative.
  - `thermal_conductivity_ice`
  - `thermal_conductivity_firn`
  - `thermal_conductivity_snow` is an archived multi-option snow
    conductivity helper.
