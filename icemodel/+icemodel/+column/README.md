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

Migration status: phase-1 state / enthalpy / thermodynamics migration is
complete; phase 2 mass-transfer, mesh, and runoff helpers are active.
Current migrated entry points:
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
  - integrates column solid mass, liquid mass, and optionally enthalpy from
    `T`, `f_ice`, `f_liq`, and `dz`, on the solver's physical intrinsic-density
    and physical-water MWE basis; `use_ro_glc` changes initialization fractions
    only
- `icemodel.namelists.budgetoutputs` (the ledger field list, in `+namelists`)
- `icemodel.column.initialize_budget_state`
  - returns the zeroed fixed-schema ledger a forcing step starts from
- `icemodel.column.accumulate_phase_budget`
- `icemodel.column.accumulate_vapor_budget`
- `icemodel.column.accumulate_remesh_budget`
- `icemodel.column.initialize_remesh_ledger`
  - returns the zeroed per-event remesh ledger. `merge_thin_layers` is its
    only caller and returns the filled struct as its opt-in eighth output;
    it lives here so that schema has one definition the tests can assert
    against
  - the three accumulators update ledger state rather than returning standalone
    event terms, which is why they are named `accumulate_*_budget`. They keep
    the diagnostic ledger out of the timestep driver; call them once per
    accepted substep, after the enthalpy solve, after the surface vapor
    exchange, and after remeshing. One exception: `accumulate_vapor_budget`
    ASSIGNS the condensation-overflow channel rather than adding to it, because
    `d_rof` is reset once per forcing step and accumulated across substeps, so
    it already arrives as the step total
- `icemodel.column.budget_surface_mass_balance`
- `icemodel.column.merge_thin_layers`
  - three views of the same remeshing export, which nest rather than
    duplicate. `df_lyr` (ice2, standard and diagnostic profiles) totals
    the mass all merges
    removed, as a water-equivalent fraction the caller scales by `dz`.
    `mass_budget_merge_export_solid/liquid_mwe` (diagnostic profile) is that
    same total split by phase, which the closure identities require.
    `mass_budget_top_export_solid/liquid_mwe` is the surface-removal subset,
    separated because interior merges move mass without lowering the grid.
  - none of the three is a surface-loss comparator. A merge gives the joined
    cell the MEAN of the pair it replaces, so removing a nearly empty top cell
    still exports about half the pair's mass. The export therefore over-counts
    what the removed cell held. The PROMICE ablation evaluation scores melt and
    the runoff diagnostics instead, and it neither scores nor plots merge
    export. `icemodel-pla` tracks the conserving remap that would fix this.
- `icemodel.column.infiltration`
- `icemodel.column.liquid_flux`
- `icemodel.column.vapor_mass_transfer`
- `icemodel.column.merge_layer_indices`
- `icemodel.column.merge_layers`
- `icemodel.column.enforce_control_volume_balance`
- `icemodel.column.available_liquid_water`
- `icemodel.column.control_volume_mesh`
- `icemodel.column.diagnose_column_runoff`
- `icemodel.column.water_fraction`
