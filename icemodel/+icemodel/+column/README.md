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
    only caller and returns the filled struct as its opt-in eighth output
  - the three accumulators update ledger state rather than returning standalone
    event terms, which is why they are named `accumulate_*_budget`. They keep
    the diagnostic ledger out of the timestep driver; call them once per
    accepted substep, after the enthalpy solve, after the surface vapor
    exchange, and after remeshing. One exception: `accumulate_vapor_budget`
    ASSIGNS the condensation-overflow channel rather than adding to it, because
    `d_rof` is reset once per forcing step and accumulated across substeps, so
    it already arrives as the step total
- `icemodel.column.budget_surface_mass_balance`
- `icemodel.column.apply_vapor_transport`
  - applies interior vapor transport, which Fick's law fixes as a MASS, so
    each cell receives exactly the mass that arrived. The same three limits
    apply, and what a cell cannot take is recorded rather than assumed. The
    two appliers are separate because they conserve different quantities:
    routing transport through the energy path lets a limited cell apply a
    different mass than arrived
- `icemodel.column.vapor_exchange_is_wet`
  - the one owner of the wet/dry decision for vapor mass exchange. Both
    appliers ask it, so neither can land a cell in a band where two criteria
    disagree. `couple_vapor_transport` does not ask it: that function moves
    mass and makes no phase decision
- `icemodel.column.vapor_face_quantities`
  - the one face rule for vapor transport: an fn-weighted harmonic mean of
    `De` with no porosity factor, and the secant `ro_vap` difference across
    the face. The mass flux and the vapor energy flux are built from these
    same face quantities, which is what makes `energy = L * mass` hold
    discretely rather than approximately
- `icemodel.column.vapor_face_conductance`
  - the energy side of that same face rule, in two pieces. The matrix part
    is the donor-tangent interface conductivity [W m-1 K-1], positive at
    every face, which the enthalpy assembly adds to its face conductivity
    in coupled mode; the deferred-correction flux [W m-2] joins the source
    vector, so the converged face flux is exactly `L_face * U_vap` while
    the assembled system keeps its diagonal dominance (Patankar 1980,
    section 7.2). The donor phase comes from `vapor_exchange_is_wet`, the
    mass applier's predicate. Both boundary faces are zero in both parts:
    the bottom is closed, and the surface exchange enters through the
    surface energy balance rather than as a diffusive flux
- `icemodel.column.couple_vapor_step`
  - the one entry point for the coupled interior vapor path. Called once per
    accepted substep, after the surface budgets close, so interior
    transport never lands in `d_liq` or in the surface vapor channels. It
    evaluates the node quantities at the state the solve converged on,
    moves vapor across the interior faces, applies the arriving mass,
    threads the gross face accumulation `d_vap_faces` in and out for
    grain growth, and records the per-phase redistribution increments and
    their shortfall. The driver carries only the accumulator and the
    pre-exchange ice and liquid snapshots (`f_ice_solve`, `f_liq_solve`);
    every other vapor intermediate lives here
- `icemodel.column.couple_vapor_transport`
  - moves vapor between the cells by Fick's law, on a mass basis. Both its
    boundaries are closed, so it redistributes and creates nothing. The
    surface exchange is deliberately not part of it: the flux divergence is
    linear in the faces, so closing the top face separates the two exactly
    and each keeps the invariant it actually has
- `icemodel.column.vapor_face_diffusivity`
  - the fn-weighted harmonic mean of the node diffusivities at each face.
    Both the mass flux and the energy conductance call it, and the discrete
    `energy = L * mass` identity holds only while they share it
- `icemodel.column.accumulate_redistribution_budget`
  - records the per-phase storage increments interior transport causes,
    and the shortfall its per-cell limits rejected. Cross-phase transport
    conserves mass while moving solid and liquid storage in opposite
    directions, so the per-phase storage closures consume these
    increments; the surface vapor closure identity never carries them
- `icemodel.column.accepted_vapor_quantities`
  - evaluates the saturation vapor density and the effective diffusivity at
    an accepted substep state, once, for the coupled path to reuse
- `icemodel.column.max_liquid_fraction_change`
  - the one owner of the largest `f_liq` increase a control volume accepts,
    `ro_ice/ro_liq * (1 - f_ice) - f_liq`, which is `f_wat_max - f_wat` on the
    `water_fraction` basis. The pore volume `1 - f_ice` is scaled to water
    equivalent as if it were ice, so the bound falls short of the pore volume
    by `(1 - ro_ice/ro_liq) * (1 - f_ice)`. That shortfall is the room the
    liquid needs to expand if it refreezes. Both vapor appliers,
    `infiltration`, and `assert_max_water` use the same bound
- `icemodel.column.potential_sublimation`
  - converts a potential vapor tendency from a liquid-water volume fraction
    to the ice volume fraction that carries the same latent-heat demand. The
    surface applier, the `merge_thin_layers` look-ahead, and
    `apply_vapor_transport`'s wet branches all call it, so a prediction
    cannot use a different factor than the application
- `icemodel.column.merge_thin_layers`
  - three views of the same remeshing export, which nest rather than
    duplicate. `df_lyr` (ice2, standard and diagnostic profiles) totals the
    mass that all merges removed, as a water-equivalent fraction the caller
    scales by `dz`.
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
  - solves the diffusive vapor flux and grows grains from it. Two calling
    conventions: with nine arguments it evaluates its own saturation
    density and diffusivity and uses a Dirichlet ghost node at `Ts` for
    the top face; with the trailing `d_vap_faces` vector it grows grains
    from the gross face exchange the accepted substeps applied, converted
    to step-mean magnitude fluxes, and evaluates no saturation state and
    no ghost node at all
- `icemodel.column.merge_layer_indices`
- `icemodel.column.merge_layers`
- `icemodel.column.enforce_control_volume_balance`
- `icemodel.column.available_liquid_water`
- `icemodel.column.control_volume_mesh`
- `icemodel.column.diagnose_column_runoff`
- `icemodel.column.water_fraction`
