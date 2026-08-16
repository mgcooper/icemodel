# icemodel.surface

Purpose: public surface-energy, surface-state, and surface mass-diagnostic
contracts.

Public entrypoints:
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
- `potential_surface_vapor_exchange`
  - phase-corrects one substep's surface vapor demand. `d_pevp` is an energy
    demand on a fixed `Lv` basis, and the same energy sublimates less mass
    than it evaporates, so this rescales it to the liquid-water volume
    fraction of the mass actually moved. It owns the wet/dry read through
    `icemodel.column.vapor_exchange_is_wet` and must run before the exchange,
    because a merge can replace the top cell
- `surface_vapor_mass_flux`
  - converts that accumulated fraction to the Neumann face flux
    [kg m-2 s-1]. It carries no latent heat: the phase correction already
    happened. Keeping the two apart is what lets the driver accumulate a
    fraction and never handle a mass
- `apply_surface_vapor_exchange`
  - applies the surface vapor exchange to the top cell under three limits:
    condensation capped by pore capacity, deposition capped by available air
    space, and sublimation floored at `f_ice_min`. The surface energy balance
    fixes that exchange as an ENERGY, so this may spend part of one demand on
    liquid and the rest on ice. Condensation the top cell cannot hold leaves
    as runoff. Its interior counterpart is
    `icemodel.column.apply_vapor_transport`, which applies a MASS instead,
    because Fick's law fixes the mass and lets the energy follow
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
