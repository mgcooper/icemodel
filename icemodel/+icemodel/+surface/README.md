# icemodel.surface

Purpose: public surface-energy, surface-state, and surface mass-diagnostic contracts.

Public entrypoints:
- `diagnose_turbulent_heat_fluxes`
- `diagnose_surface_energy_balance`
- `solve_surface_energy_balance`
- `surface_energy_balance_residual`
- `surface_energy_balance_terms`
- `surface_flux_linearization`
- `physical_surface_temperature`
- `apply_surface_vapor_mass_change`
- `diagnose_surface_ablation`
- `diagnose_surface_runoff`
- `incoming_longwave_radiation`
- `empirical_incoming_longwave_radiation`
- `incoming_shortwave_radiation`
- `terrain_adjusted_shortwave_radiation`
- `atmospheric_pressure_from_elevation`

Allowed dependencies:
- `icemodel.surface.turbulence.*`
- shared physics kernels and model constants
- `icemodel.numerics.*` for generic math only

Rules:
- keep public call sites here short and stable
- push scheme-specific turbulence details down into `+turbulence`
- do not reintroduce root `SEB*` or `fSEB` compatibility wrappers

Migration status: active and mostly complete for the THF/SEB stack and the
surface-owned mass-diagnostic helpers on `feat/bulk-mo-thf-scheme`.
