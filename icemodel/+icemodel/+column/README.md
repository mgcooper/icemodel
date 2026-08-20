# icemodel.column

Purpose: vertical-column physics, state evolution, phase transforms,
hydrology, vapor transport, mass-balance bookkeeping, and mesh ownership.

Planned public areas:
- column-state initialization
- enthalpy/temperature solves
- shortwave source-term assembly
- mesh/layer-merge helpers
- liquid transport / infiltration
- vapor transport
- runoff diagnostics
- mass-balance bookkeeping

Rules:
- own column physics here rather than in root all-caps files
- keep runtime control in `icemodel.timestepping`, not here

Entry points:
- `icemodel.column.shortwave_source_term`
- `icemodel.column.initialize_column_state`
- `icemodel.column.solve_column_enthalpy`
- `icemodel.column.solve_column_temperature`
- `icemodel.column.assemble_enthalpy_system`
- `icemodel.column.meltzone_transform`
- `icemodel.column.liquid_fraction_derivative`
- `icemodel.column.liquid_fraction_function`
- `icemodel.column.bulk_enthalpy`
- `icemodel.column.updatestate`
- `icemodel.column.bulk_thermal_conductivity`
- `icemodel.column.firn_thermal_conductivity`
- `icemodel.column.integrate_column_budget`
  - integrates column solid mass, liquid mass, and optional enthalpy from
    `T`, `f_ice`, `f_liq`, and `dz`. It uses the solver's physical intrinsic
    density and physical-water MWE basis. `use_ro_glc` changes only the
    initialization fractions
- `icemodel.namelists.budgetoutputs` (the budget channel list, in `+namelists`)
- `icemodel.column.initialize_budget_state`
  - returns the zeroed 22-channel budget a forcing step starts from and
    records the storage start endpoints from the entry state. The transient
    `budget.substep` struct carries within-substep baselines and is not a
    channel; the driver strips it at output emission
- `icemodel.column.finalize_budget_state`
  - records the storage end endpoints after the substep loop
- `icemodel.column.accumulate_phase_budget`
- `icemodel.column.accumulate_vapor_budget`
- `icemodel.column.accumulate_remesh_budget`
- `icemodel.column.initialize_remesh_ledger`
  - returns the zeroed event ledger for one remesh pass. `merge_thin_layers`
    is its only caller and folds the filled ledger into the budget through
    `accumulate_remesh_budget`
  - the accumulators update budget state rather than returning standalone
    event terms, which is why they are named `accumulate_*_budget`. The
    budget accumulates on every profile; only output emission is bound to
    the diagnostic profile. `accumulate_phase_budget` runs in the driver;
    `accumulate_vapor_budget` runs inside `budget_surface_mass_balance`;
    `accumulate_redistribution_budget` runs inside `couple_vapor_step`;
    `accumulate_remesh_budget` runs inside `merge_thin_layers`. One
    exception: `accumulate_vapor_budget`
    ASSIGNS the condensation-overflow channel rather than adding to it, because
    `d_rof` is reset once per forcing step and accumulated across substeps, so
    it already arrives as the step total
- `icemodel.column.budget_surface_mass_balance`
- `icemodel.column.apply_vapor_transfer`
  - the one phase-state mutator for surface exchange and interior transport.
    It applies signed, phase-resolved liquid-water-equivalent increments under
    the shared liquid and ice storage limits. Phase selection and energy-to-mass
    conversion happen upstream. Rejected increments return on their original
    phase basis
- `icemodel.column.vapor_shortfall_ice_equivalent`
  - converts rejected liquid- and ice-phase LWE increments to the common
    ice-fraction energy basis. The surface vapor budget uses this basis for
    its unapplied-energy channel
- `icemodel.column.vapor_exchange_is_wet`
  - the one owner of the residual-mobility wet/dry decision for vapor
    exchange. Surface
    demand partitioning, accepted-state interior transfer, and face latent heat
    use this predicate
- `icemodel.column.vapor_transport_terms`
  - the one face construction for the enthalpy solve. It applies `fn`-weighted
    harmonic interpolation to the vapor-free node conductivity and effective
    vapor diffusivity. It returns `k_eff_faces`, `k_vap_faces`,
    `q_vap_deferred_faces`, `U_vap_faces`, and the face donor latent heat
    `L_vap_faces`. Both vapor boundary faces are closed. The matrix and
    deferred terms preserve `Q_vap = L_face * U_vap_faces`. The solver
    returns the accepted mass flux and its donor latent heat to the surface
    coupler and driver without a second calculation. This vapor-specific
    constructor does not own liquid infiltration or `U_liq`
- `icemodel.column.couple_vapor_step`
  - the production entry point for interior vapor transfer. It runs once per
    accepted substep after the surface budgets close. It converts the accepted
    `U_vap` face flux to node increments and routes each face's mass to the
    phase whose latent heat the accepted sweep carried (`L_vap == Lv` marks a
    wet donor). It calls `apply_vapor_transfer` and records the redistribution
    increments. Interior transport never enters `d_liq` or the surface vapor
    channels
- `icemodel.column.accumulate_redistribution_budget`
  - records the per-phase storage increments from interior transport.
    Cross-phase transport conserves
    mass while moving solid and liquid storage in opposite directions. The
    per-phase storage closures include these increments. The surface vapor
    closure identity excludes them. Transport rejected by per-cell limits is
    visible through `apply_vapor_transfer`'s unapplied outputs
- `icemodel.column.max_liquid_fraction_change`
  - the one owner of the largest `f_liq` increase a control volume accepts,
    `ro_ice/ro_liq * (1 - f_ice) - f_liq`, which is `f_wat_max - f_wat` on the
    `water_fraction` basis. The pore volume `1 - f_ice` is scaled to water
    equivalent as if it were ice. The bound is less than the pore volume by
    `(1 - ro_ice/ro_liq) * (1 - f_ice)`. That difference gives liquid room to
    expand if it refreezes. `apply_vapor_transfer`, `infiltration`, and
    `assert_max_water` use the same bound
- `icemodel.column.potential_sublimation`
  - converts a potential vapor demand from a liquid-water volume fraction
    to the ice volume fraction that carries the same latent-heat demand. The
    surface wrapper and the `merge_thin_layers` look-ahead call it, so every
    ice-fraction diagnostic uses the same conversion
- `icemodel.column.merge_thin_layers`
  - provides three nested views of one remeshing export. `df_lyr` (ice2,
    standard and diagnostic profiles) totals the mass removed by all merges.
    The caller scales this water-equivalent fraction by `dz`.
    `mass_budget_merge_export_solid_mwe` (diagnostic profile) is the solid
    part of that total, which the solid closure identity requires.
    `mass_budget_top_export_solid/liquid_mwe` is the surface-removal subset,
    separated because interior merges move mass without lowering the grid.
  - none of the three is a surface-loss comparator. A merge gives the joined
    cell the
    MEAN of the pair that it replaces. Removing a nearly empty top cell still
    exports about half the pair's mass. The export therefore exceeds the mass
    in the removed cell. The PROMICE ablation evaluation scores melt and runoff
    diagnostics. It does not score or plot merge export. `icemodel-pla` tracks
    the conserving remap that would correct this behavior
- `icemodel.column.infiltration`
- `icemodel.column.liquid_flux`
- `icemodel.column.update_grain_radius`
  - advances thermal grain radius once per accepted substep from that
    substep's face-flux magnitudes, the realized surface exchange, and the
    concurrent liquid fraction. It does not diagnose vapor transport or
    saturation
- `icemodel.column.merge_layer_indices`
- `icemodel.column.merge_layers`
- `icemodel.column.enforce_control_volume_balance`
- `icemodel.column.available_liquid_water`
- `icemodel.column.control_volume_mesh`
- `icemodel.column.diagnose_column_runoff`
- `icemodel.column.water_fraction`
