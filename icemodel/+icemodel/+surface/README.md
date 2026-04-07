# icemodel.surface

Purpose: public surface-energy and surface-state contracts.

Public entrypoints:
- `diagnose_turbulent_heat_fluxes`
- `diagnose_surface_energy_fluxes`
- `solve_surface_energy_balance`
- `surface_energy_balance_residual`
- `surface_energy_balance_terms`
- `surface_flux_linearization`
- `physical_surface_temperature`

Allowed dependencies:
- `icemodel.surface.turbulence.*`
- shared physics kernels and model constants
- `icemodel.numerics.*` for generic math only

Rules:
- keep public call sites here short and stable
- push scheme-specific turbulence details down into `+turbulence`
- do not reintroduce root `SEB*` or `fSEB` compatibility wrappers

Migration status: active and mostly complete for the THF/SEB stack on `feat/bulk-mo-thf-scheme`.
