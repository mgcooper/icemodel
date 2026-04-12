# icemodel-lbk Repo-wide namespace reorganization

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This document must be maintained in accordance with `.agents/PLANS.md` and `.agents/plans/PLANS.local.md`.

## Purpose / Big Picture

After this change, the production physics code will no longer be concentrated in the root `icemodel/` folder behind legacy all-caps entrypoints. Instead, the codebase will expose semantic namespaces for surface physics, surface turbulence, column physics, vapor thermodynamics, radiation, timestepping, numerics, and coupling. A contributor should be able to navigate directly to the owning namespace for a behavior, and the repo should no longer rely on compatibility shims or thin wrappers once the migration is complete. The first user-visible proof is that the active turbulent-flux and surface-energy-balance stack on `feat/bulk-mo-thf-scheme` moves into the new `icemodel.surface` and `icemodel.surface.turbulence` layout while continuing to pass the focused unit suites.

## Progress

- [x] (2026-04-07 03:41Z) Reviewed `bd ready`; confirmed there was no existing ready issue for the namespace reorganization and that the work must be opened explicitly.
- [x] (2026-04-07 03:42Z) Created the umbrella Bead `icemodel-lbk` plus follow-up Beads `icemodel-bwq`, `icemodel-vik`, `icemodel-hrx`, `icemodel-x6z`, `icemodel-npa`, `icemodel-bto`, `icemodel-yvb`, `icemodel-r43`, `icemodel-2ex`, and `icemodel-gdk` so another agent can resume without thread history.
- [x] (2026-04-07 03:44Z) Wrote the authoritative ExecPlan in `.agents/plans/execution/`.
- [x] (2026-04-07) Marked the umbrella Bead `icemodel-lbk` in progress. As of 2026-04-09 it remains `IN_PROGRESS`.
- [x] (2026-04-07) Completed the initial surface/turbulence, coupler, timestepping, numerics, vapor, and radiation namespace extraction passes that moved the active physics stack out of the root all-caps entrypoints and removed the temporary wrappers for those completed slices.
- [x] (2026-04-07) Closed Beads `icemodel-yvb` and `icemodel-bto` after the vapor and radiation extraction slices landed. The umbrella Bead remains open; `icemodel-x6z`, `icemodel-npa`, `icemodel-r43`, `icemodel-2w6`, and `icemodel-gdk` remain open.
- [x] (2026-04-07) Renamed the Monin-Obukhov scheme identifier from `'bulk_mo'` to `'monin_obukhov'` across `setopts`, `configureRun`, production dispatch, and tests.
- [x] (2026-04-07) Renamed the vapor namespace entry points from short fused names to descriptive names:
  - `saturation_vapor_pressure`
  - `relative_humidity_from_vapor_pressure`
  - `saturation_vapor_density`
  - `vapor_diffusivity`
  - `vapor_thermal_diffusion_coefficient`
  - `dew_point_temperature`
  - `wet_bulb_temperature`
  - `initialize_vapor_model`
  while keeping `moist_air_density` and `specific_humidity_from_vapor_pressure` in `icemodel.vapor`.
- [x] (2026-04-07) Renamed the radiation namespace entry points from the initial fused-lowercase extraction names to the current descriptive names:
  - `get_scattering_coefficients`
  - `get_solar_spectrum`
  - `initialize_spectral_model`
  - `update_extinction_coefficients`
  - `spectral_extinction_coefficients`
  - `rescale_spectral_extinction_coefficients`
  - `bulk_extinction_coefficients`
  - `bulk_extinction_coefficients_lookup`
  - `solvetwostream`
  - `smooth_twostream_fluxes`
  and kept the column-owned absorbed-shortwave assembler at `icemodel.column.shortwave_source_term`.
- [x] (2026-04-08) Completed the remaining core surface migration items required by the first revision request: migrated `CONDUCT`, `ENBAL`, `MFENERGY`, and `QADVECT` into `icemodel.surface`; moved `bulk_richardson/solve_surface_temperature` to `icemodel.surface.solve_surface_temperature`; renamed `surface_fluxes` to `evaluate_surface_flux`; removed the passthrough `T_sfc` output from `surface_energy_balance_terms`; updated production callers, couplers, tests, and docs; and created Bead `icemodel-2w6` plus the matching `open-items.md` theme for `scalar_exchange_diagnostics`.
- [x] (2026-04-08) Unified the SEB call-contract argument ordering across the surface stack, couplers, tests, and debug helpers around the canonical atmospheric ordering `ea_atm, br_coefs, roL, liqflag, chi`.
- [x] (2026-04-08) Renamed `METINIT` to `icemodel.surface.initialize_surface_forcings` and `LOADMETDATA` to `icemodel.timestepping.getforcings`, then updated production callers, fixtures, tests, and the surrounding forcing-flow documentation. The separate De / `br_coefs` computation also moved out of the old met-loader path and into the main model flow.
- [x] (2026-04-08) Added `icemodel.surface.turbulence.bulk_richardson.richardson_number` and then extended it on 2026-04-09 to support both canonical and precomputed-coefficient fast-path use.
- [x] (2026-04-09) Completed the second revision-request surface cleanup: `conductive_heat_flux` now returns `dQc/dT_sfc`; `solve_surface_temperature` uses explicit `Qc` and `dQc/dT_sfc` rather than the older legacy branches; `evaluate_surface_energy_balance` now accepts precomputed `Qsn` and `Qln`; `surface_energy_balance_terms` exposes `Qsn` and `Qln`; the Robin/Dirichlet conduction distinction is documented in both bulk-Richardson and Monin-Obukhov linearization helpers; and `icemodel.m` now passes current-step state into the Dirichlet coupler instead of stale `x*` copies.
- [x] (2026-04-09) Added `icemodel.helpers.rmttleapinds` and the KANM/KANL THF scheme-comparison infrastructure (`run_thf_cases`, `plot_thf_cases`, and the thin `validate_thf_cases` wrapper).
- [x] (2026-04-08 to 2026-04-09) Recorded successful native MATLAB validation in git history. The latest explicit validation evidence is that all 16 unit test files (157 tests total) passed after the vapor/radiation rename pass, and later commits continued updating tests to the revised SEB/turbulence contracts.
- [x] (2026-04-09) Re-synced the active continuation state for `icemodel-x6z`, re-audited the remaining root column phase-1 stack (`ICEINIT`, `ICEENBAL`, `SKINSOLVE`, `GECOEFS`, `MZTRANSFORM`, `FREEZECURVE`, `MELTCURVE`, `TOTALHEAT`, `TOTALHEAT2`, `UPDATESTATE`, `BULKTHERMALK`, `THERMALK`), and confirmed the current production/test blast radius before editing.
- [x] (2026-04-09) Moved the column phase-1 files into `icemodel.column` using the selected target names and updated the direct production/test call sites away from the removed root all-caps entry points:
  - `initialize_column_state`
  - `solve_column_enthalpy`
  - `solve_column_temperature`
  - `assemble_enthalpy_system`
  - `meltzone_transform`
  - `liquid_fraction_derivative`
  - `liquid_fraction_function`
  - `total_enthalpy`
  - `total_enthalpy_2`
  - `updatestate`
  - `bulk_thermal_conductivity`
  - `firn_thermal_conductivity`
- [x] (2026-04-09) Verified the column phase-1 move slice with `git diff --check`, `rg` confirmation that no direct `ICEINIT(...)` / `ICEENBAL(...)` / `SKINSOLVE(...)` / related root call sites remain in production or tests, and Octave `which(...)` resolution of the new `icemodel.column` entry points. MATLAB is available via the explicit app binary at `/Applications/MATLAB_R2025b.app/bin/matlab`; rerun the focused native-MATLAB unit suites before closing `icemodel-x6z`.
- [x] (2026-04-09) Completed the agreed descriptive rename pass for the moved column thermodynamics helpers: `bulkthermalk -> bulk_thermal_conductivity`, `thermalk -> firn_thermal_conductivity`, `freezecurve -> liquid_fraction_derivative`, `meltcurve -> liquid_fraction_function`, `totalheat -> total_enthalpy`, and `totalheat2 -> total_enthalpy_2`. `liquid_fraction_function` and the inlined derivative paths in `solve_column_enthalpy`, `updatestate`, and `total_enthalpy_2` now call `liquid_fraction_derivative`.
- [x] (2026-04-09) Reran focused native MATLAB validation for the phase-1 column rename pass with `/Applications/MATLAB_R2025b.app/bin/matlab`. The following suites all passed: `test_phase_and_state_kernels`, `test_flux_and_thermo_kernels`, `test_surface_solver_kernels`, `test_core_thermo_and_timestep`, and `test_reduced_solver_runs` (48/48 tests passed).
- [x] (2026-04-09) Completed the SEB-suite contract cleanup requested after the column phase-1 pass: `surface_energy_balance_residual` now owns the conductive-flux evaluation and accepts `T_ice`, `k_eff`, and `dz` instead of a precomputed `Qc`; `surface_energy_balance_terms` now returns only net fluxes (`Qsn`, `Qln`) plus the direct conductive / advective / turbulent terms; `diagnose_melt_freeze_energy` consumes `Qsn` and `Qln` directly for the melt-cap path; and the repeated inlined shortwave / longwave formulas were centralized into `icemodel.surface.net_shortwave_radiation` plus `icemodel.surface.net_longwave_radiation`.
- [x] (2026-04-09) Removed the legacy root `LONGOUT` helper by moving it to `icemodel.surface.net_longwave_radiation`, added `icemodel.surface.net_shortwave_radiation`, and updated the SEB diagnostics, couplers, debug tooling, Monin-Obukhov linearization helper, and focused tests to the new surface-owned radiation bookkeeping.
- [x] (2026-04-09) Reran focused native MATLAB validation for the SEB cleanup with `/Applications/MATLAB_R2025b.app/bin/matlab`. The following suites all passed: `test_flux_and_thermo_kernels`, `test_surface_solver_kernels`, `test_turbulent_flux_schemes`, and `test_reduced_solver_runs` (43/43 tests passed).
- [x] (2026-04-10) Simplified the top-level SEB diagnostic contract by renaming `diagnose_surface_energy_fluxes` to `diagnose_surface_energy_balance`, removing `diag_turbulent` from `surface_energy_balance_terms`, returning the full term set `[Qe, Qh, Qc, Qsn, Qln, Qa, Qm, Qf, Qbal]`, and collecting `thf_diag` directly from `diagnose_turbulent_heat_fluxes` only in the `output_profile='diagnostic'` branches of `icemodel.m` and `skinmodel.m`.
- [x] (2026-04-10) Reran focused native MATLAB validation for the SEB contract rename / simplification pass with `/Applications/MATLAB_R2025b.app/bin/matlab`. The following suites all passed: `test_flux_and_thermo_kernels`, `test_surface_solver_kernels`, `test_turbulent_flux_schemes`, and `test_reduced_solver_runs` (43/43 tests passed).
- [x] (2026-04-10) Opened and completed Bead `icemodel-vhk` to migrate the legacy forcing-derivation fallbacks into `icemodel.surface`: `LONGIN -> incoming_longwave_radiation`, `PRESSURE -> atmospheric_pressure_from_elevation`, `SOLARIN -> incoming_shortwave_radiation`, and `SOLARRAD -> terrain_adjusted_shortwave_radiation`. The helpers now carry proper docstrings, code comments, namespaced call sites, and focused tests. `incoming_shortwave_radiation` also fixes the long-standing undefined `Qsi_tmp` bug by explicitly averaging hourly terrain-adjusted flux samples.
- [x] (2026-04-10) Reran focused native MATLAB validation for the forcing-fallback migration with `/Applications/MATLAB_R2025b.app/bin/matlab`. The following suites all passed: `test_flux_and_thermo_kernels`, `test_spectral_kernels`, and `test_met_contracts` (39/39 tests passed).
- [x] (2026-04-10) Approved the explicit phase-2 rename map before editing, then completed the column/surface move pass with `git mv`: `ICEMF -> icemodel.column.budget_surface_mass_balance`, `INFILTRATION -> icemodel.column.infiltration`, `LIQFLUX -> icemodel.column.liquid_flux`, `VAPORTRANSFER -> icemodel.column.vapor_mass_transfer`, `LAYERINDS -> icemodel.column.merge_layer_indices`, `COMBINEHEAT -> icemodel.column.merge_layers`, `VOLBAL -> icemodel.column.enforce_control_volume_balance`, `LIQAVAIL -> icemodel.column.available_liquid_water`, `CVMESH -> icemodel.column.control_volume_mesh`, `ICERUNOFF -> icemodel.column.diagnose_column_runoff`, `ICESUBL -> icemodel.surface.apply_surface_vapor_mass_change`, `ICEABLATION -> icemodel.surface.diagnose_surface_ablation`, and `SRFRUNOFF -> icemodel.surface.diagnose_surface_runoff`.
- [x] (2026-04-10) Canonicalized the remaining SEB duplication inside the turbulence helpers without widening the public solver contracts: `bulk_richardson/evaluate_surface_flux` now uses `icemodel.numerics.complex_step_derivative` plus the canonical `net_shortwave_radiation`, `net_longwave_radiation`, and `evaluate_surface_energy_balance` helpers; `bulk_richardson/surface_flux_linearization` now computes its current-state flux through the canonical surface-budget helpers while retaining its analytic Jacobian; and `monin_obukhov/surface_flux_linearization` now reuses the shared numerics derivative helper instead of keeping a private complex-step implementation.
- [x] (2026-04-10) Tightened the bulk-Richardson numerical-reference path and longwave helper ownership after the SEB review: `bulk_richardson/evaluate_surface_flux` is now a thin complex-step wrapper around `icemodel.surface.surface_energy_balance_residual` using raw conductive inputs instead of a frozen `Qc`; `solve_surface_temperature` now reuses `incoming_longwave_radiation`, `outgoing_longwave_radiation`, `evaluate_surface_energy_balance`, and the new optional derivative outputs from `bulk_richardson.latent_heat_flux` / `sensible_heat_flux`; and the old forcing fallback `incoming_longwave_radiation(Tair, ea_atm, ...)` was renamed to `empirical_incoming_longwave_radiation` so `incoming_longwave_radiation(Qli)` can be the canonical absorbed-incoming-longwave SEB helper.
- [x] (2026-04-10) Split the old `ICEMF` responsibilities so `icemodel.column.budget_surface_mass_balance` now owns only the mass-balance bookkeeping, while the new `icemodel.column.merge_thin_layers` helper performs the thin-layer merge/remesh step explicitly in the main `icemodel.m` flow immediately afterward.
- [x] (2026-04-10) Reran targeted native MATLAB validation for the phase-2 rename/move pass using `/Applications/MATLAB_R2025b.app/bin/matlab` with explicit source-path setup (`addpath(fullfile(pwd,'icemodel')); addpath(fullfile(pwd,'test'));`). The following suites all passed: `test_surface_phase_partitioning`, `test_mesh_and_timestep_kernels`, `test_phase_and_state_kernels`, and `test_reduced_solver_runs` (32/32 tests passed).
- [x] (2026-04-11) Internalized the remaining physical-constant and parameter lookups that were still being threaded from `icemodel.m` / `skinmodel.m` into downstream column, coupler, radiation, timestep, and surface-mass helpers. The active main-model callers no longer pass `ro_*`, latent heats, phase-curve parameters, or spectral `ro_ice` into downstream helpers; the affected functions now own those lookups persistently via `physicalConstant`, `parameterLookup`, or the new `icemodel.column.meltzone_bounds` helper.
- [x] (2026-04-11) Added `icemodel.water_fraction`, `icemodel.column.meltzone_bounds`, `icemodel.column.surface_linearization_error`, and `icemodel.column.subsurface_linearization_error`; updated `liquid_fraction_function`, `liquid_fraction_derivative`, `meltzone_transform`, `total_enthalpy`, `total_enthalpy_2`, `assemble_enthalpy_system`, `solve_column_enthalpy`, `solve_column_temperature`, `budget_surface_mass_balance`, `merge_thin_layers`, `apply_surface_vapor_mass_change`, `checksubstep`, `updatesubstep`, and `initialize_spectral_model` to use the new canonical ownership boundaries; and replaced inlined `f_wat = f_liq + f_ice * ro_ice / ro_liq` expressions with `icemodel.water_fraction(...)` where those calculations were still duplicated.
- [x] (2026-04-11) Reran focused native MATLAB validation for the constant-lookup and caller-contract cleanup with `/Applications/MATLAB_R2025b.app/bin/matlab` using explicit source-path setup. The following suites all passed: `test_phase_and_state_kernels`, `test_mesh_and_timestep_kernels`, `test_surface_phase_partitioning`, `test_surface_solver_kernels`, and `test_flux_and_thermo_kernels` (58/58 tests passed).
- [ ] Close or retitle the still-`IN_PROGRESS` Beads `icemodel-bwq`, `icemodel-hrx`, and `icemodel-2ex` if their scope is already fully satisfied by the landed work.

## Surprises & Discoveries

- Observation: the local Beads status-update step is currently unreliable even though reads and creates succeed.
  Evidence: `bd update icemodel-lbk --status in_progress` and `bd update icemodel-bwq --status in_progress` hit a local Dolt circuit breaker on 2026-04-07 after the issues were created successfully.

- Observation: the Beads state is readable and updatable when run outside the read-only sandbox, but not from the current read-only sandbox without escalation.
  Evidence: on 2026-04-09, `bd ready` failed in-sandbox because it could not create `.beads/dolt-server.lock`, then succeeded immediately with escalation; `bd show` also confirmed the current namespace-reorg Bead states once run with escalation.

- Observation: Octave can resolve the new namespaced entrypoints with `which(...)`, but it is not a reliable substitute for executing package-qualified `icemodel.*` calls in this repo because the root `icemodel.m` entrypoint collides with Octave package resolution.
  Evidence: `octave --no-init-file --quiet --eval "addpath(fullfile(pwd,'icemodel')); disp(which('icemodel.surface.solve_surface_energy_balance'))"` resolves the moved files successfully on 2026-04-07, while attempting direct package-qualified execution routes into the root `icemodel.m` entrypoint and fails before reaching the namespaced function bodies.

- Observation: the user wants namespaced function names to follow the existing MATLAB-style naming already present in the repo rather than a blanket snake_case rule.
  Evidence: after the first-pass commit on 2026-04-07, the user explicitly restored preferred names such as `initialize_timesteps`, `checksubstep`, `aitkenscalar`, and `fsearchzero`, and directed that existing camelCase names should remain unchanged.

- Observation: current values in `setopts` are not reliable evidence of long-term project defaults because they are often changed interactively during local development.
  Evidence: on 2026-04-09, the user explicitly clarified that periodically changing `solver` and similar `setopts` values is part of their interactive workflow and should not be recorded in the ExecPlan as an authoritative design decision.

- Observation: the ExecPlan was only partially maintained after the two surface revision requests.
  Evidence: the current file already contained genuine 2026-04-08 surface-migration updates when inspected on 2026-04-09, but it still described the pre-rename vapor/radiation names, omitted the `METINIT` / `LOADMETDATA` rename, omitted the later SEB-contract cleanup, and still carried stale "matlab not on PATH" language despite later commit history recording native MATLAB test runs.

- Observation: direct `matlab -batch "runtests(...)"` from repo root is not sufficient for this repo because the `icemodel/` source root is not automatically on MATLAB's path.
  Evidence: the first 2026-04-10 phase-2 validation attempt failed with `undefinedVarOrClass` errors for `icemodel.*` packages until the rerun added `addpath(fullfile(pwd,'icemodel')); addpath(fullfile(pwd,'test'));` explicitly.

## Decision Log

- Decision: keep the public surface call sites short by using `icemodel.surface.diagnose_turbulent_heat_fluxes(...)` and `icemodel.surface.diagnose_surface_energy_balance(...)` as the stable top-level surface contracts.
  Rationale: the user wants deep scheme namespaces to remain internal details rather than long production call sites, while still keeping the function names descriptive enough to stand alone semantically in a future Python translation.
  Date/Author: 2026-04-07 / Codex

- Decision: use `icemodel.surface.turbulence.monin_obukhov` as the namespace and `'monin_obukhov'` as the `opts.turbulent_flux_scheme` option string (replacing the prior `'bulk_mo'` abbreviation).
  Rationale: both the namespace and the option string should use the same physically descriptive name. The prior note "the namespace should describe the closure family rather than the current option string abbreviation" was incomplete — the option string abbreviation also needed to be updated to match. Updated 2026-04-07.
  Date/Author: 2026-04-07 / Codex

- Decision: function naming follows a meaning-first rule rather than blanket snake_case or blanket lowercase fusion. Short established helpers can remain single-word lowercase (`checksubstep`, `aitkenscalar`, `fsearchzero`), while opaque or abbreviated legacy names should be expanded to descriptive names even when that means snake_case (`initialize_vapor_model`, `bulk_extinction_coefficients`, `relative_humidity_from_vapor_pressure`).
  Date/Author: 2026-04-07 / user correction

- Decision: move the timestep/substep control helpers into `icemodel.timestepping` applying the syllable-count naming rule above.
  Rationale: these functions coordinate full-step and substep runtime control rather than column physics, even when they transport column state through the accepted/rejected-state workflow.
  Date/Author: 2026-04-07 / Codex

- Decision: move the Ambaum/Romps vapor thermodynamics and atmosphere-surface vapor conversion helpers into `icemodel.vapor`, while treating surface vapor mass application and column vapor transport as distinct contracts.
  Rationale: these functions are a coherent shared thermodynamic family used directly by surface, column, and test code, while `icemodel.surface.apply_surface_vapor_mass_change` owns the surface closure and `icemodel.column.vapor_mass_transfer` owns the column transport / grain-growth scaffold rather than a reusable vapor property helper.
  Date/Author: 2026-04-07 / Codex

- Decision: move the optical-property and two-stream helper stack into `icemodel.radiation`, but keep the absorbed-shortwave thermal-grid assembly at `icemodel.column.shortwave_source_term`.
  Rationale: the extinction/scattering tables, spectral model initialization, and two-stream machinery are general radiation services; the mapping from spectral net flux to the thermal-grid source term remains column-owned.
  Date/Author: 2026-04-07 / Codex

- Decision: keep the Dirichlet and Robin conduction treatments distinct in the public documentation and helper contracts.
  Rationale: in the standalone Dirichlet `solve_surface_temperature` path, `Qc` and `dQc/dT_sfc` belong directly in the Newton residual and Jacobian; in the Robin path, conduction enters through the top-node equation / `icemodel.column.assemble_enthalpy_system` linearization instead of as an explicit surface Jacobian term. The April 9 surface revision pass codified that distinction in code comments and signatures.
  Date/Author: 2026-04-09 / derived from revision-request implementation

- Decision: the second revision-request pass established that `icemodel.m` should pass current-step state into `solve_surface_column_dirichlet` rather than stale `x*` state from the external call site.
  Rationale: the coupler owns substep-state bookkeeping internally; requiring stale `x*` names at the external call boundary was unnecessary and obscured the actual state flow. This does not imply that whatever `setopts` currently contains should be treated as a stable project default.
  Date/Author: 2026-04-09 / derived from revision-request implementation

- Decision: the end state for this task includes zero permanent compatibility shims.
  Rationale: the user explicitly wants a tight end-state architecture. Temporary aliases are allowed only while a migration slice is in progress and must be removed before `icemodel-lbk` closes.
  Date/Author: 2026-04-07 / Codex

- Decision: use `icemodel.column.control_volume_mesh` as the explicit mesh-construction owner and `icemodel.column.budget_surface_mass_balance` as the phase-2 mass-balance owner.
  Rationale: `control_volume_mesh` preserves the control-volume specificity without the verbosity of `build_control_volume_mesh`, and `budget_surface_mass_balance` matches the scientific ownership of the old `ICEMF` bookkeeping step better than a generic `budget_mass_balance` name. The layer-remesh loop was intentionally split into `icemodel.column.merge_thin_layers` so the bookkeeping and remeshing steps are explicit in the main model flow.
  Date/Author: 2026-04-10 / Codex + user approval

- Decision: canonicalize SEB residual assembly and complex-step differentiation rather than letting each turbulence helper keep a private implementation.
  Rationale: `bulk_richardson/evaluate_surface_flux` is only a numerical-reference helper used by tests, so it should reduce to a thin wrapper around the canonical `surface_energy_balance_residual` instead of carrying its own residual contract. A shared `icemodel.numerics.complex_step_derivative` helper keeps the complex-step logic consistent between `complexstep`, the Monin-Obukhov Robin linearization, and the bulk-Richardson derivative reference while leaving the analytic Jacobian ownership in `solve_surface_temperature`.
  Date/Author: 2026-04-10 / Codex + user discussion

- Decision: make `incoming_longwave_radiation` the canonical absorbed-longwave SEB helper and rename the meteorological fallback to `empirical_incoming_longwave_radiation`.
  Rationale: the old name was overloaded between two different concerns: a forcing-estimation formula and the actual SEB term. Splitting them lets the surface-energy code use a clean canonical helper while preserving the translated-Fortran empirical fallback for sparse station forcing datasets.
  Date/Author: 2026-04-10 / Codex + user discussion

## Outcomes & Retrospective

This section will summarize the completed migration slices, the final namespace ownership achieved in this pass, the remaining open Beads if the full repo-wide move is not completed, and the exact validation evidence. It must also confirm that any temporary aliases introduced during the work were deleted before closeout.

Current completed slices:

- Surface public API: `diagnose_turbulent_heat_fluxes`, `diagnose_surface_energy_balance`, `solve_surface_energy_balance`, `surface_energy_balance_residual`, `surface_flux_linearization`, `physical_surface_temperature`, `conductive_heat_flux`, `advective_heat_flux`, `evaluate_surface_energy_balance`, `diagnose_melt_freeze_energy`, `solve_surface_temperature`, `surface_energy_balance_terms`, `initialize_surface_forcings`, `incoming_longwave_radiation`, `incoming_shortwave_radiation`, `terrain_adjusted_shortwave_radiation`, and `atmospheric_pressure_from_elevation` now live under `icemodel.surface`.
- Surface radiation bookkeeping: `net_shortwave_radiation` and `net_longwave_radiation` now live under `icemodel.surface`, and the old root `LONGOUT` helper has been removed.
- Surface turbulence internals: the bulk-Richardson helpers now live under `icemodel.surface.turbulence.bulk_richardson` (including `evaluate_surface_flux` and `richardson_number`), and the Monin-Obukhov helpers now live under `icemodel.surface.turbulence.monin_obukhov`.
- Couplers: the active names are now `solve_surface_column_dirichlet` and `solve_surface_column_robin`.
- Timestepping/numerics: the active runtime helpers now live under `icemodel.timestepping` and `icemodel.numerics`, including `getforcings`, `checksubstep`, `resetsubstep`, `updatesubstep`, `initialize_timesteps`, `newtimestep`, `nexttimestep`, `update_force_advance_guard`, `trisolve`, `aitkenscalar`, `fsearchzero`, `complexstep`, and `sign_or_one`.
- Vapor thermodynamics: the Ambaum/Romps vapor-property and atmosphere-surface conversion helpers now live under `icemodel.vapor` using the current descriptive names (`saturation_vapor_pressure`, `saturation_vapor_density`, `vapor_diffusivity`, `vapor_thermal_diffusion_coefficient`, `dew_point_temperature`, `wet_bulb_temperature`, `initialize_vapor_model`, and `relative_humidity_from_vapor_pressure`).
- Radiation / column shortwave: the spectral-model and two-stream helpers now live under `icemodel.radiation` using the current descriptive names (`initialize_spectral_model`, `get_scattering_coefficients`, `get_solar_spectrum`, `update_extinction_coefficients`, `spectral_extinction_coefficients`, `rescale_spectral_extinction_coefficients`, `bulk_extinction_coefficients`, `bulk_extinction_coefficients_lookup`, `solvetwostream`, and `smooth_twostream_fluxes`), and the absorbed-shortwave thermal-grid assembler now lives at `icemodel.column.shortwave_source_term`.
- Column phase 1: the state-initialization, enthalpy-solve, temperature-solve, melt-zone, and thermodynamics helpers now live under `icemodel.column` as `initialize_column_state`, `solve_column_enthalpy`, `solve_column_temperature`, `assemble_enthalpy_system`, `meltzone_transform`, `liquid_fraction_derivative`, `liquid_fraction_function`, `total_enthalpy`, `total_enthalpy_2`, `updatestate`, `bulk_thermal_conductivity`, and `firn_thermal_conductivity`.
- Column phase 2: the mass-balance, liquid-transport, mesh, and runoff helpers now live under `icemodel.column` as `budget_surface_mass_balance`, `merge_thin_layers`, `infiltration`, `liquid_flux`, `vapor_mass_transfer`, `merge_layer_indices`, `merge_layers`, `enforce_control_volume_balance`, `available_liquid_water`, `control_volume_mesh`, and `diagnose_column_runoff`.
- Column support helpers: `icemodel.water_fraction`, `icemodel.column.meltzone_bounds`, `icemodel.column.surface_linearization_error`, and `icemodel.column.subsurface_linearization_error` now own the canonical water-fraction transform, melt-zone bounds, and enthalpy linearization diagnostics used by the moved column solver stack.
- Surface mass-diagnostic slice: the surface-owned mass helpers now live under `icemodel.surface` as `apply_surface_vapor_mass_change`, `diagnose_surface_ablation`, and `diagnose_surface_runoff`.
- Helpers: `icemodel.helpers.rmttleapinds` now exists and is available to namespaced validation / plotting code.
- Validation / analysis: the THF scheme-comparison infrastructure for KANM/KANL now exists in `test/interactive/`.
- Remaining root production files are now limited to unresolved support / special-case files:
  - still-root special cases / support: `METSUB`, `SKINEBSOLVE`, `SPECINIT`, `SAVEOUTPUT`, `WRITEOUTPUT`

Current Bead state relevant to this umbrella task:

- `icemodel-lbk` is `IN_PROGRESS`
- closed: `icemodel-bto`, `icemodel-yvb`
- closed: `icemodel-bto`, `icemodel-vik`, `icemodel-x6z`, `icemodel-yvb`
- closed: `icemodel-vhk`
- still `IN_PROGRESS`: `icemodel-bwq`, `icemodel-hrx`, `icemodel-2ex`
- still `OPEN`: `icemodel-npa`, `icemodel-r43`, `icemodel-2w6`, `icemodel-gdk`

Current validation evidence:

- Latest explicit recorded native-MATLAB validation in git history:
  - commit `776fe4c` recorded all 157 unit tests passing after the radiation extraction / rename pass
  - commit `3916fec` recorded all 157 unit tests passing after the vapor rename and THF scheme-identifier rename pass
- Later commits through 2026-04-09 continued updating production callers and tests (`ffb0413`, `7373c6e`) to the revised SEB / turbulence contracts, but the reconstructed plan should still be treated as needing a fresh targeted rerun before closing any remaining namespace-reorg Beads.
- On 2026-04-09 during the column phase-1 continuation pass, `git diff --check` completed cleanly and Octave `which(...)` resolved the new `icemodel.column` entry points on disk (`initialize_column_state`, `solve_column_enthalpy`, `solve_column_temperature`, `assemble_enthalpy_system`, `bulk_thermal_conductivity`, `firn_thermal_conductivity`). MATLAB is available from the shell via `/Applications/MATLAB_R2025b.app/bin/matlab`, so the next verification step should use that explicit binary instead of assuming `PATH` coverage.
- On 2026-04-09 after the descriptive rename cleanup, focused native MATLAB validation passed via `/Applications/MATLAB_R2025b.app/bin/matlab` for `test_phase_and_state_kernels`, `test_flux_and_thermo_kernels`, `test_surface_solver_kernels`, `test_core_thermo_and_timestep`, and `test_reduced_solver_runs` (48 tests total).
- On 2026-04-09 after the SEB-suite cleanup, focused native MATLAB validation passed via `/Applications/MATLAB_R2025b.app/bin/matlab` for `test_flux_and_thermo_kernels`, `test_surface_solver_kernels`, `test_turbulent_flux_schemes`, and `test_reduced_solver_runs` (43 tests total).
- On 2026-04-10 after the `diagnose_surface_energy_balance` rename / simplification pass, focused native MATLAB validation again passed via `/Applications/MATLAB_R2025b.app/bin/matlab` for `test_flux_and_thermo_kernels`, `test_surface_solver_kernels`, `test_turbulent_flux_schemes`, and `test_reduced_solver_runs` (43 tests total).
- On 2026-04-10 after migrating the legacy forcing-derivation fallbacks into `icemodel.surface`, focused native MATLAB validation passed via `/Applications/MATLAB_R2025b.app/bin/matlab` for `test_flux_and_thermo_kernels`, `test_spectral_kernels`, and `test_met_contracts` (39 tests total).
- On 2026-04-10 after the column phase-2 rename / move pass, focused native MATLAB validation passed via `/Applications/MATLAB_R2025b.app/bin/matlab` with explicit source-path setup (`addpath(fullfile(pwd,'icemodel')); addpath(fullfile(pwd,'test'));`) for `test_surface_phase_partitioning`, `test_mesh_and_timestep_kernels`, `test_phase_and_state_kernels`, and `test_reduced_solver_runs` (32 tests total).
- On 2026-04-11 after the constant-lookup / caller-contract cleanup, focused native MATLAB validation passed via `/Applications/MATLAB_R2025b.app/bin/matlab` with explicit source-path setup (`addpath(fullfile(pwd,'icemodel')); addpath(fullfile(pwd,'test'));`) for `test_phase_and_state_kernels`, `test_mesh_and_timestep_kernels`, `test_surface_phase_partitioning`, `test_surface_solver_kernels`, and `test_flux_and_thermo_kernels` (58 tests total).
- On 2026-04-09, `bd ready` reported the next ready work as `icemodel-2w6`, `icemodel-gdk`, `icemodel-r43`, `icemodel-npa`, and `icemodel-x6z`.
- On 2026-04-10, `bd show` confirmed:
  - `icemodel-lbk` is `IN_PROGRESS`
  - `icemodel-bto`, `icemodel-vik`, `icemodel-x6z`, and `icemodel-yvb` are `CLOSED`
  - `icemodel-npa`, `icemodel-r43`, and `icemodel-2w6` are still open as the next namespace-reorg / follow-up slices.

## Context and Orientation

The repo already contains several stable support namespaces under `icemodel/+icemodel/`, including `run`, `netcdf`, `plot`, `test`, `internal`, `namelists`, and `validators`. The unsettled area is the production physics stack, which still relies heavily on root-level all-caps files in `icemodel/`.

The surface, turbulence, vapor, radiation, timestepping, numerics, and coupler slices are now materially moved. Another agent should not restart those migrations from scratch. The remaining heavy lift is the vertical-column root stack plus the final root-cleanup / residual-shim sweep.

The current public namespace owners on disk are:

- `icemodel.surface`
- `icemodel.surface.turbulence.bulk_richardson`
- `icemodel.surface.turbulence.monin_obukhov`
- `icemodel.couplers`
- `icemodel.timestepping`
- `icemodel.numerics`
- `icemodel.vapor`
- `icemodel.radiation`
- `icemodel.column`
  Current migrated entry points:
- `icemodel.column.initialize_column_state`
- `icemodel.column.solve_column_enthalpy`
- `icemodel.column.solve_column_temperature`
- `icemodel.column.assemble_enthalpy_system`
- `icemodel.column.meltzone_transform`
- `icemodel.column.liquid_fraction_derivative`
- `icemodel.column.liquid_fraction_function`
- `icemodel.column.total_enthalpy`
- `icemodel.column.total_enthalpy_2`
- `icemodel.column.updatestate`
- `icemodel.column.bulk_thermal_conductivity`
- `icemodel.column.firn_thermal_conductivity`
- `icemodel.column.budget_surface_mass_balance`
- `icemodel.column.merge_thin_layers`
- `icemodel.column.infiltration`
- `icemodel.column.liquid_flux`
- `icemodel.column.vapor_mass_transfer`
- `icemodel.column.merge_layer_indices`
- `icemodel.column.merge_layers`
- `icemodel.column.enforce_control_volume_balance`
- `icemodel.column.available_liquid_water`
- `icemodel.column.control_volume_mesh`
- `icemodel.column.diagnose_column_runoff`
- `icemodel.surface.apply_surface_vapor_mass_change`
- `icemodel.surface.diagnose_surface_ablation`
- `icemodel.surface.diagnose_surface_runoff`

The current tests already reflect the intended subsystem seams:

- `test/unit/test_turbulent_flux_schemes.m`
- `test/unit/test_surface_solver_kernels.m`
- `test/unit/test_flux_and_thermo_kernels.m`
- `test/unit/test_mesh_and_timestep_kernels.m`
- `test/unit/test_spectral_kernels.m`
- `test/unit/test_run_contracts.m`

The desired target layout for this task is:

- `icemodel.surface` for public surface contracts
- `icemodel.surface.turbulence.bulk_richardson` and `icemodel.surface.turbulence.monin_obukhov` for internal scheme implementations
- `icemodel.column` for vertical-column physics
- `icemodel.couplers` for surface-column coupling
- `icemodel.vapor` for vapor thermodynamics and atmospheric vapor conversions
- `icemodel.radiation` for optical and radiative-transfer machinery
- `icemodel.timestepping` for accepted/rejected-state timestep control
- `icemodel.numerics` for reusable math algorithms

The active Beads for this work are:

- `icemodel-lbk` umbrella reorganization task
- `icemodel-bwq` THF scheme namespace completion
- `icemodel-hrx` timestepping extraction
- `icemodel-npa` column phase 2
- `icemodel-bto` radiation extraction
- `icemodel-yvb` vapor extraction
- `icemodel-r43` compatibility-shim deletion
- `icemodel-2ex` docs/tests/function-signature updates
- `icemodel-gdk` diagnostics/measurements follow-up evaluation
- `icemodel-2w6` scalar-exchange-diagnostics completion follow-up from the first revision-request pass

## Plan of Work

The next agent should treat the surface/turbulence, vapor, radiation, timestepping, numerics, and coupler slices as already substantially complete, then verify and clean up their Bead states rather than redoing them.

The critical remaining engineering work is now:

1. Finish the root cleanup and compatibility sweep (`icemodel-r43`).
   The old column phase-2 all-caps roots are now gone. The next pass should
   decide final ownership for the remaining root support/special-case files:
   - `METSUB`
   - `SKINEBSOLVE`
   - `SPECINIT`
   - `SAVEOUTPUT`
   - `WRITEOUTPUT`

2. Reconcile Bead state with reality.
   `icemodel-npa` should now be closable after a final Beads update. The
   still-`IN_PROGRESS` Beads for completed slices (`icemodel-bwq`,
   `icemodel-hrx`, `icemodel-2ex`) should be reviewed and either closed or
   explicitly narrowed.

3. Keep `icemodel-2w6` and `icemodel-gdk` as follow-up work, not blockers.
   `scalar_exchange_diagnostics` and the diagnostics/measurements namespace discussion were surfaced by the revision-request passes and belong in the continuation context, but neither should block the column/root cleanup sequence.

## Concrete Steps

All commands are run from `/Users/mattcooper/MATLAB/projects/icemodel`.

1. Re-sync Beads and the ExecPlan before writing code.

      bd ready
      bd show icemodel-lbk
      bd show icemodel-x6z
      bd show icemodel-npa
      bd show icemodel-r43
      bd show icemodel-2w6

   Then explicitly decide which still-`IN_PROGRESS` namespace-reorg Beads can be closed immediately and record that in this file.

2. Audit the remaining root support files and choose the next cleanup slice.

      find icemodel -maxdepth 1 -type f | sort
      rg -n "METSUB|SKINEBSOLVE|SPECINIT|SAVEOUTPUT|WRITEOUTPUT" icemodel test doc

   Prefer a cleanup cut that keeps the repo runnable and testable rather than moving every remaining support file at once.

3. Move the chosen remaining support slice into the appropriate namespace, update callers/tests/docs/comments, and remove the root aliases immediately for the completed slice.

4. Run targeted validation after each cleanup slice, starting with the unit tests already closest to the touched files:

      /Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/mattcooper/MATLAB/projects/icemodel'); addpath(fullfile(pwd,'icemodel')); addpath(fullfile(pwd,'test')); results = runtests({'test/unit/test_phase_and_state_kernels.m','test/unit/test_core_thermo_and_timestep.m','test/unit/test_reduced_solver_runs.m'}); disp(table({results.Name}',[results.Passed]',[results.Failed]','VariableNames',{'Name','Passed','Failed'})); assert(all([results.Passed]));"

   Add `test/unit/test_spectral_kernels.m` when touching the remaining radiation/shortwave integration points, and keep `test/unit/test_surface_solver_kernels.m` in the loop whenever a surface/column boundary changes.

5. Once the remaining cleanup slices are in place, do a final `rg` pass over `icemodel/` root files and close `icemodel-r43` only when the namespace-reorg leftovers are either moved or explicitly documented as intentional exceptions.

## Validation and Acceptance

Acceptance for continuation from this point requires:

1. Another agent can read this file plus `bd show` for the open Beads and understand that the next core implementation work is `icemodel.column`, not a restart of surface/turbulence/vapor/radiation.
2. The current completed namespaces are described using the actual current names on disk:
   - surface / turbulence
   - couplers
   - timestepping / numerics
   - vapor
   - radiation
   - `column.shortwave_source_term`
   - `column` phase-1 state / enthalpy / thermodynamics helpers
3. The still-open Beads are clearly separated into:
   - remaining core namespace work (`x6z`, `npa`, `r43`)
   - follow-up / design work (`2w6`, `gdk`)
4. The latest recorded validation evidence and its limits are stated accurately:
   - native MATLAB test runs were recorded in git history
   - this column phase-1 continuation pass ran `git diff --check`, Octave `which(...)` resolution checks, and a focused native-MATLAB rerun
   - use `/Applications/MATLAB_R2025b.app/bin/matlab` for in-shell native MATLAB verification when `PATH` is incomplete
5. No future agent should infer from this file that the old short vapor / radiation names (`vappress`, `getscattercoefs`, etc.) are still current. The descriptive names are now canonical.

## Idempotence and Recovery

The Beads and ExecPlan steps are safe to repeat. Namespace moves should be performed in slices that each keep the repo runnable and testable. If a partial move breaks callers, restore the last passing namespace ownership for that slice and retry with a narrower move rather than leaving duplicate permanent APIs. Temporary aliases are acceptable only while a slice is incomplete; once a slice is complete, delete them immediately.

Because the worktree may evolve during this long-running task, read any touched file again before applying a follow-up patch. Do not revert unrelated user changes.

## Artifacts and Notes

Initial Bead creation transcript:

    bd ready
    ✨ No open issues

    bd create "Repo-wide namespace reorganization umbrella"
    icemodel-lbk

    bd create "Finish THF scheme namespace completion"
    icemodel-bwq

Original status-update blocker:

    bd update icemodel-lbk --status in_progress
    Error resolving icemodel-lbk: dolt circuit breaker is open: server appears down

This blocker was transient. By 2026-04-09, `bd ready` and `bd show` succeeded again when run outside the read-only sandbox.

## Interfaces and Dependencies

At the end of the active surface/turbulence migration, the production-facing surface interface must include:

- `icemodel.surface.diagnose_turbulent_heat_fluxes(...)`
- `icemodel.surface.diagnose_surface_energy_balance(...)`
- `icemodel.surface.solve_surface_energy_balance(...)`
- `icemodel.surface.surface_energy_balance_residual(...)`
- `icemodel.surface.surface_flux_linearization(...)`

The deep turbulence implementations must live under:

- `icemodel.surface.turbulence.bulk_richardson`
- `icemodel.surface.turbulence.monin_obukhov`

The completed active slices in this task must not retain old root aliases such as `SEBSOLVE`, `SEBFLUX`, `STABLEFN`, `LATENT`, `SENSIBLE`, or `CHECKSUBSTEP`.

Revision note: created at task start for Bead `icemodel-lbk` and updated from the final approved user plan to reflect the locked namespace decisions, the Bead breakdown, and the explicit no-permanent-shims rule.
