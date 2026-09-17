# icemodel.column

Purpose: Vertical-column physics.

- column-state initialization
- enthalpy/temperature solves
- shortwave source-term assembly
- remeshing/layer-merging
- liquid transport / infiltration
- vapor transport
- runoff diagnostics
- mass-balance bookkeeping

Contents:

- `icemodel.column.shortwave_source_term`
- `icemodel.column.initialize_column_state`
- `icemodel.column.solve_column_enthalpy`
- `icemodel.column.solve_column_temperature`
- `icemodel.column.assemble_enthalpy_system`
- `icemodel.column.meltzone_transform`
- `icemodel.column.liquid_fraction_derivative`
- `icemodel.column.liquid_fraction_function`
- `icemodel.column.bulk_enthalpy`
- `icemodel.column.bulk_thermal_conductivity`
- `icemodel.column.firn_thermal_conductivity`
- `icemodel.column.integrate_column_budget`
  - Integrates column solid and liquid mass on a MWE basis and optional
    enthalpy in J m-2 from `T_ice`, `f_ice`, `f_liq`, and `dz`.
- `icemodel.column.initialize_budget_state`
  - Returns the initialized 21-channel budget for one forcing step and records
    the storage start endpoints from the entry state. On each substep, applier
    functions pass per-substep phase change increments to `accumulate_*`
    functions, which the budget records over the full forcing step.
- `icemodel.column.finalize_budget_state`
  - Records the storage end endpoints after the substep loop.
- `icemodel.column.accumulate_phase_budget`
  - Records the solid and liquid storage changes produced by the accepted
    enthalpy solve. The model driver calls it.
- `icemodel.column.accumulate_vapor_exchange`
  - Records the surface vapor demand, applied solid and liquid changes, and
    condensation overflow. It stores `d_rof * dz(1)` as the overflow because
    `d_rof` contains the forcing-step total through the current substep.
    `budget_surface_mass_balance` calls it.
- `icemodel.column.accumulate_vapor_transport`
  - Records the applied solid and liquid storage changes from interior vapor
    transport. `couple_vapor_step` calls it.
- `icemodel.column.accumulate_remesh_budget`
  - Records the numerical storage changes in the remesh event ledger.
    `merge_thin_layers` calls it.
- `icemodel.column.initialize_remesh_ledger`
  - Returns the zeroed event ledger that `merge_thin_layers` fills during one
    remesh pass and sends to `accumulate_remesh_budget`.
- `icemodel.column.budget_surface_mass_balance`
- `icemodel.column.apply_vapor_transfer`
  - Applies signed liquid- and ice-phase vapor increments for surface exchange
    and interior transport. Control-volume water capacity limits liquid
    addition. `f_res` limits liquid removal. Control-volume water capacity
    limits ice deposition, and `f_ice_min` limits sublimation. The function
    returns each rejected phase change increment on a liquid-water-equivalent
    basis.
- `icemodel.column.vapor_exchange_is_wet`
  - The residual-mobility wet/dry decision for vapor exchange.
    `potential_surface_vapor_exchange` uses it to partition the surface
    demand, and `vapor_transport_terms` uses it to pick each face's donor
    latent heat.
- `icemodel.column.vapor_transport_terms`
  - Constructs vapor transport terms at control volume faces for the enthalpy
    solve. It applies `fn`-weighted harmonic interpolation to the vapor-free
    node conductivity and effective vapor diffusivity. It returns `k_eff_faces`,
    `k_vap_faces`, `q_vap_deferred_faces`, `U_vap_faces`, and the face donor
    latent heat `L_vap_faces`. Both vapor boundary faces are closed. The matrix
    and deferred terms preserve `Q_vap = L_face * U_vap_faces`. The solver
    returns the accepted mass flux and its donor latent heat to the surface
    coupler and driver.
- `icemodel.column.couple_vapor_step`
  - Computes interior vapor transfer once per accepted substep after the
    surface budgets close. It converts the `U_vap` face flux to node increments
    using the `L_vap` value returned with that flux. It calls
    `apply_vapor_transfer` and records the applied transport increments.
- `icemodel.column.max_liquid_fraction_change`
  - Returns the largest `f_liq` increase a control volume accepts:
    `ro_ice/ro_liq * (1 - f_ice) - f_liq`, which is `f_wat_max - f_wat` on the
    `water_fraction` basis. The pore volume `1 - f_ice` is scaled to water
    equivalent as if it were ice. The bound is less than the pore volume by
    `(1 - ro_ice/ro_liq) * (1 - f_ice)`. That difference gives liquid room to
    expand if it refreezes. `apply_vapor_transfer`, `infiltration`, and
    `assert_max_water` use the same bound.
- `icemodel.column.potential_sublimation`
  - Converts a potential vapor demand from a liquid-water volume fraction
    to an ice volume fraction with the same latent-heat demand.
    `merge_thin_layers` calls it to predict whether surface demand will take
    a layer below the minimum allowed ice fraction on the next substep,
    triggering a layer merge.
- `icemodel.column.merge_thin_layers`
  - Combines two equal-thickness cells and stores their mean state in the
    surviving cell. This operation removes half of the pair's mass. The
    function shifts the remaining cells and clones the bottom cell.
  - `df_lyr` records the water-equivalent fraction removed by each merge.
    Scale each cell by `dz` to obtain metres water equivalent.
    `mass_budget_merge_export_solid_mwe` records the solid part.
    `mass_budget_top_export_solid_mwe` and
    `mass_budget_top_export_liquid_mwe` record the surface-removal subset.
    `mass_budget_top_deletion_height_m` records the associated surface
    lowering. Pair averaging and bottom-cell cloning can also change stored
    mass because the remap is not conservative. `icemodel-pla` tracks the
    conservative-remap follow-up.
- `icemodel.column.infiltration`
- `icemodel.column.liquid_flux`
- `icemodel.column.update_grain_radius`
  - Updates thermal grain radius using the vapor transfer face-flux magnitudes,
    the surface vapor exchange, and the current liquid fraction.
- `icemodel.column.updatestate`
- `icemodel.column.merge_layer_indices`
- `icemodel.column.merge_layers`
- `icemodel.column.enforce_control_volume_balance`
- `icemodel.column.available_liquid_water`
- `icemodel.column.control_volume_mesh`
- `icemodel.column.diagnose_column_runoff`
- `icemodel.column.water_fraction`
