# icemodel.surface

Purpose: public surface-energy, surface-state, and surface mass-diagnostic
contracts.

Public entry points:
- `diagnose_turbulent_heat_fluxes`
- `diagnose_surface_energy_balance`
- `solve_surface_energy_balance`
- `surface_energy_balance_residual`
- `surface_energy_balance_terms`
- `numerical_surface_flux_linearization`
- `surface_flux_linearization`
- `physical_surface_temperature`
- `diagnose_surface_ablation`
- `diagnose_surface_runoff`
- `potential_surface_vapor_demand`
  - diagnoses turbulent latent heat at the physical surface temperature and
    converts it to the liquid-water-equivalent energy demand `d_pevp`
- `potential_surface_vapor_exchange`
  - partitions one `d_pevp` demand into signed liquid- and ice-phase
    liquid-water-equivalent increments before storage limits. Wet evaporation
    exhausts mobile liquid before converting the remaining energy demand to
    ice at `Lv/Ls`. Wet condensation targets liquid, and dry exchange targets
    ice. The partition closes the mixed-latent energy identity exactly
- `surface_vapor_mass_flux`
  - converts a liquid-water volume fraction to the surface face flux
    [kg m-2 s-1]. It carries no latent heat: any phase correction happens
    upstream. Keeping conversion and correction apart is what lets callers
    hold fractions and never handle a mass
- `apply_surface_vapor_exchange`
  - owns surface boundary policy. It partitions the energy demand, invokes
    `icemodel.column.apply_vapor_transfer` for the shared state limits, routes
    rejected positive liquid condensation to runoff, and translates realized
    and unapplied increments to the surface budget terms
- `empirical_incoming_longwave_radiation`
- `outgoing_longwave_radiation`
- `net_longwave_radiation`
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
