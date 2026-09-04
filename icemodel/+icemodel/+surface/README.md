# icemodel.surface

Purpose: surface-energy, surface-state, and surface mass-diagnostic functions.

Contents:

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
  - Diagnoses turbulent latent heat at the physical surface temperature and
    converts it to the liquid-water-equivalent energy demand `d_pevp`.
- `potential_surface_vapor_exchange`
  - Partitions `d_pevp` into liquid and ice increments that satisfy
    `Lv * d_pevp = Lv * d_vap_liq + Ls * d_vap_ice_lwe` before storage limits
    are applied.
- `surface_vapor_mass_flux`
  - Converts a liquid-water volume fraction to the surface face mass flux
    [kg m-2 s-1].
- `apply_surface_vapor_exchange`
  - Partitions the surface vapor-energy demand and calls
    `icemodel.column.apply_vapor_transfer` to enforce liquid and ice storage
    limits. Evaporation exhausts mobile liquid before converting the remaining
    energy demand to ice at `Lv/Ls`. Rejected liquid condensation becomes
    runoff. Top-cell capacity limits ice deposition.
    Removal that the top cell cannot supply passes to the cells below as a
    latent-energy-equivalent demand, converted between cell thicknesses.
    Applied liquid and ice increments are returned for each cell. `d_vap`
    contains the realized surface exchange as a liquid-water-equivalent
    fraction of the top cell.
- `empirical_incoming_longwave_radiation`
- `outgoing_longwave_radiation`
- `net_longwave_radiation`
- `incoming_shortwave_radiation`
- `terrain_adjusted_shortwave_radiation`
- `atmospheric_pressure_from_elevation`
