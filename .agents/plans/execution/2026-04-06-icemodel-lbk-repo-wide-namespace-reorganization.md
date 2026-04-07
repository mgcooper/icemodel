# icemodel-lbk Repo-wide namespace reorganization

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This document must be maintained in accordance with `.agents/PLANS.md` and `.agents/plans/PLANS.local.md`.

## Purpose / Big Picture

After this change, the production physics code will no longer be concentrated in the root `icemodel/` folder behind legacy all-caps entrypoints. Instead, the codebase will expose semantic namespaces for surface physics, surface turbulence, column physics, vapor thermodynamics, radiation, timestepping, numerics, and coupling. A contributor should be able to navigate directly to the owning namespace for a behavior, and the repo should no longer rely on compatibility shims or thin wrappers once the migration is complete. The first user-visible proof is that the active turbulent-flux and surface-energy-balance stack on `feat/bulk-mo-thf-scheme` moves into the new `icemodel.surface` and `icemodel.surface.turbulence` layout while continuing to pass the focused unit suites.

## Progress

- [x] (2026-04-07 03:41Z) Reviewed `bd ready`; confirmed there was no existing ready issue for the namespace reorganization and that the work must be opened explicitly.
- [x] (2026-04-07 03:42Z) Created the umbrella Bead `icemodel-lbk` plus follow-up Beads `icemodel-bwq`, `icemodel-vik`, `icemodel-hrx`, `icemodel-x6z`, `icemodel-npa`, `icemodel-bto`, `icemodel-yvb`, `icemodel-r43`, `icemodel-2ex`, and `icemodel-gdk` so another agent can resume without thread history.
- [x] (2026-04-07 03:44Z) Wrote the authoritative ExecPlan in `.agents/plans/execution/`.
- [ ] Mark `icemodel-lbk` and the active implementation slice `icemodel-bwq` in progress once the transient local Beads/Dolt circuit-breaker clears.
- [x] (2026-04-07 06:05Z) Completed the active surface and surface-turbulence migration: the public entry points now live at `icemodel.surface.diagnose_turbulent_heat_fluxes`, `icemodel.surface.diagnose_surface_energy_fluxes`, `icemodel.surface.solve_surface_energy_balance`, and `icemodel.surface.physical_surface_temperature`; the scheme implementations moved under `icemodel.surface.turbulence.bulk_richardson` and `icemodel.surface.turbulence.monin_obukhov`; the old root THF/SEB files and the legacy `fSEB.m` / `ENBALANCE.m` shims were removed.
- [x] (2026-04-07 05:20Z) Extracted the active timestep/substep helpers and generic math helpers into `icemodel.timestepping` and `icemodel.numerics`, updating the direct runtime callers and removing the root copies.
- [x] (2026-04-07 06:22Z) Normalized the active `timestepping` and `numerics` entry points to the current package naming style (`checksubstep`, `resetsubstep`, `updatesubstep`, `initialize_timesteps`, `newtimestep`, `nexttimestep`, `update_force_advance_guard`, `trisolve`, `aitkenscalar`, `fsearchzero`, `complexstep`, `sign_or_one`) and updated the declarations, callers, tests, comments, and README references accordingly.
- [x] (2026-04-07 05:45Z) Renamed the active coupling entry points to `icemodel.couplers.solve_surface_column_dirichlet` and `icemodel.couplers.solve_surface_column_robin`, updating production and test call sites without leaving compatibility wrappers.
- [x] (2026-04-07 06:34Z) Extracted the vapor thermodynamics slice into `icemodel.vapor`: `vappress`, `vappress2rh`, `vapordensity`, `vapordiffusivity`, `vapork`, `tdewpoint`, `twetbulb`, `vaporinit`, `moist_air_density`, and `specific_humidity_from_vapor_pressure` now live in the namespace, and production callers, tests, demos, and comments were updated to use the new names directly.
- [ ] Decide how much of `column`, `vapor`, and `radiation` can be completed safely in this pass on the current feature branch; if any phase is deferred, keep its Bead open and update this ExecPlan accordingly.
- [x] (2026-04-07 06:05Z) Removed the temporary root aliases for the completed active slices; there are no remaining compatibility wrappers for the finished THF/SEB, coupler, timestepping, or numerics migrations.
- [ ] Run the focused unit suites for the touched namespaces and record the results here. Native MATLAB is still unavailable in this shell, so only Octave/lookup smoke checks have run so far.

## Surprises & Discoveries

- Observation: the local Beads status-update step is currently unreliable even though reads and creates succeed.
  Evidence: `bd update icemodel-lbk --status in_progress` and `bd update icemodel-bwq --status in_progress` hit a local Dolt circuit breaker on 2026-04-07 after the issues were created successfully.

- Observation: Octave can resolve the new namespaced entrypoints with `which(...)`, but it is not a reliable substitute for executing package-qualified `icemodel.*` calls in this repo because the root `icemodel.m` entrypoint collides with Octave package resolution.
  Evidence: `octave --no-init-file --quiet --eval "addpath(fullfile(pwd,'icemodel')); disp(which('icemodel.surface.solve_surface_energy_balance'))"` resolves the moved files successfully on 2026-04-07, while attempting direct package-qualified execution routes into the root `icemodel.m` entrypoint and fails before reaching the namespaced function bodies.

- Observation: the user wants namespaced function names to follow the existing MATLAB-style naming already present in the repo rather than a blanket snake_case rule.
  Evidence: after the first-pass commit on 2026-04-07, the user explicitly restored preferred names such as `initialize_timesteps`, `checksubstep`, `aitkenscalar`, and `fsearchzero`, and directed that existing camelCase names should remain unchanged.

## Decision Log

- Decision: keep the public surface call sites short by using `icemodel.surface.diagnose_turbulent_heat_fluxes(...)` and `icemodel.surface.diagnose_surface_energy_fluxes(...)` as the stable top-level surface contracts.
  Rationale: the user wants deep scheme namespaces to remain internal details rather than long production call sites, while still keeping the function names descriptive enough to stand alone semantically in a future Python translation.
  Date/Author: 2026-04-07 / Codex

- Decision: use `icemodel.surface.turbulence.monin_obukhov` rather than `icemodel.surface.turbulence.bulk_mo`.
  Rationale: the user prefers the physically descriptive Monin-Obukhov name, and the namespace should describe the closure family rather than the current option string abbreviation.
  Date/Author: 2026-04-07 / Codex

- Decision: move `check_substep`, `reset_substep`, `update_substep`, `init_timesteps`, `new_timestep`, `next_step`, and `update_force_advance_guard` into `icemodel.timestepping`.
- Decision: move the timestep/substep control helpers into `icemodel.timestepping` using the repo's existing MATLAB-style lowercase naming rather than enforcing blanket snake_case.
  Rationale: these functions coordinate full-step and substep runtime control rather than column physics, even when they transport column state through the accepted/rejected-state workflow.
  Date/Author: 2026-04-07 / Codex

- Decision: move the Ambaum/Romps vapor thermodynamics and atmosphere-surface vapor conversion helpers into `icemodel.vapor`, keeping column microstructure transport such as `VAPORTRANSFER` outside the namespace for now.
  Rationale: these functions are a coherent shared thermodynamic family used directly by surface, column, and test code, while `VAPORTRANSFER` remains column-owned transport physics rather than a reusable vapor property helper.
  Date/Author: 2026-04-07 / Codex

- Decision: the end state for this task includes zero permanent compatibility shims.
  Rationale: the user explicitly wants a tight end-state architecture. Temporary aliases are allowed only while a migration slice is in progress and must be removed before `icemodel-lbk` closes.
  Date/Author: 2026-04-07 / Codex

## Outcomes & Retrospective

This section will summarize the completed migration slices, the final namespace ownership achieved in this pass, the remaining open Beads if the full repo-wide move is not completed, and the exact validation evidence. It must also confirm that any temporary aliases introduced during the work were deleted before closeout.

Current completed slices:

- Surface public API: `diagnose_turbulent_heat_fluxes`, `diagnose_surface_energy_fluxes`, `solve_surface_energy_balance`, `surface_energy_balance_residual`, `surface_flux_linearization`, and `physical_surface_temperature` now live under `icemodel.surface`.
- Surface turbulence internals: the bulk-Richardson helpers now live under `icemodel.surface.turbulence.bulk_richardson`, and the Monin-Obukhov helpers now live under `icemodel.surface.turbulence.monin_obukhov`.
- Couplers: the active names are now `solve_surface_column_dirichlet` and `solve_surface_column_robin`.
- Timestepping/numerics: the active runtime helpers now live under `icemodel.timestepping` and `icemodel.numerics`.
- Vapor thermodynamics: the Ambaum/Romps vapor-property and atmosphere-surface conversion helpers now live under `icemodel.vapor`.

Current validation evidence:

- `matlab` is not on `PATH` in this shell as of 2026-04-07, so the native MATLAB unit suites listed below could not be run here.
- `octave --no-init-file --quiet --eval "addpath(fullfile(pwd,'icemodel')); disp(which('icemodel.surface.diagnose_turbulent_heat_fluxes')); disp(which('icemodel.surface.solve_surface_energy_balance')); disp(which('icemodel.couplers.solve_surface_column_robin')); disp(which('icemodel.timestepping.checksubstep')); disp(which('icemodel.numerics.trisolve')); disp(which('icemodel.vapor.vappress')); disp(which('icemodel.vapor.vapordensity'));"` resolved the moved namespace entrypoints successfully on 2026-04-07.
- Faraday completed Octave smoke checks for `trisolve`, `aitkenscalar`, `fsearchzero`, `resetsubstep`, `nexttimestep`, `update_force_advance_guard`, and `checksubstep`; `updatesubstep` remained blocked there by Octave package-resolution behavior rather than by the namespace move itself.
- `rg -n --glob '!test/baselines/**' --glob '!**/*.html' "\\b(VAPPRESS2RH|VAPPRESS|VAPORDENSITY|VAPORDIFFUSIVITY|VAPORK|TDEWPOINT|TWETBULB|VAPORINIT)\\b|icemodel\\.kernels\\.(moist_air_density|specific_humidity_from_vapor_pressure)" icemodel test doc` returned no matches after the vapor migration on 2026-04-07, confirming that production, tests, and docs no longer reference the old root vapor entry points or the old `icemodel.kernels` humidity-helper paths.

## Context and Orientation

The repo already contains several stable support namespaces under `icemodel/+icemodel/`, including `run`, `netcdf`, `plot`, `test`, `internal`, `namelists`, and `validators`. The unsettled area is the production physics stack, which still relies heavily on root-level all-caps files in `icemodel/`.

The active production surface stack already has a partial package migration. `icemodel/+icemodel/+surface/` contains helper contracts for turbulent heat flux, SEB term assembly, forcing snow-depth resolution, roughness selection, and Robin linearization. However, the root files `SEBSOLVE.m`, `SEBFLUX.m`, `SFCFLUX.m`, `SFCFLIN.m`, `SFCTEMP.m`, `STABLEFN.m`, `LATENT.m`, `SENSIBLE.m`, and `WINDCOEF.m` still own the bulk of the production surface/turbulence behavior. The current public THF dispatch helper is `icemodel.surface.turbulent_heat_flux(...)`, and the scheme-specific implementations are `turbulent_heat_flux_bulk_richardson.m` and `turbulent_heat_flux_bulk_mo.m`.

The current coupling layer already lives in `icemodel/+icemodel/+couplers/`, with `solve_surface_subsurface_dirichlet.m` and `solve_surface_subsurface_robin.m`. The current tests already reflect several desired subsystem boundaries: `test/unit/test_turbulent_flux_schemes.m`, `test_surface_solver_kernels.m`, `test_flux_and_thermo_kernels.m`, `test_mesh_and_timestep_kernels.m`, and `test_spectral_kernels.m`.

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
- `icemodel-vik` SEB public surface contracts
- `icemodel-hrx` timestepping extraction
- `icemodel-x6z` column phase 1
- `icemodel-npa` column phase 2
- `icemodel-bto` radiation extraction
- `icemodel-yvb` vapor extraction
- `icemodel-r43` compatibility-shim deletion
- `icemodel-2ex` docs/tests/function-signature updates
- `icemodel-gdk` diagnostics/measurements follow-up evaluation

## Plan of Work

First, complete the active branch-critical surface and turbulence slice. Create the new namespace scaffolding under `icemodel/+icemodel/+surface/+turbulence/+bulk_richardson/` and `+monin_obukhov/`. Move the scheme-local helpers there, rename them to the final semantic names, and replace the current public THF entrypoint with `icemodel.surface.diagnose_turbulent_heat_fluxes(...)`. Update every touched caller so production code no longer reaches the old root helper names.

Second, finish the public surface contract cleanup. Rename `SEBSOLVE` to `icemodel.surface.solve_surface_energy_balance`, rename `SEBFLUX` to `icemodel.surface.diagnose_surface_energy_fluxes`, keep `surface_energy_balance_residual` as the nonlinear root-find contract, and remove `ENBALANCE` and `fSEB` once all callers and tests use the namespaced surface contracts directly.

Third, extract the active control-flow helpers into `timestepping` and the generic math helpers into `numerics`. This includes the root timestep/substep controller functions, `update_force_advance_guard`, and generic algorithms such as Brent root finding, Aitken acceleration, and the tridiagonal solver. Update the tests that currently describe these as “timestep kernels” or “mesh and timestep kernels” so they point at the new namespaces.

Fourth, continue outward only after the active surface/turbulence stack is stable: migrate the coupler names from `surface_subsurface` to `surface_column`, then move the owned column, radiation, and vapor files into their final semantic homes. If a later phase cannot be completed safely in this pass, keep the affected Bead open, record the exact state here, and do not leave any long-term compatibility wrappers behind.

Fifth, run focused unit validation after each major slice, then do a final pass to remove any temporary aliases or root legacy files that were only used during migration. The umbrella task cannot close while those remain.

## Concrete Steps

All commands are run from `/Users/mattcooper/MATLAB/projects/icemodel`.

1. Task-tracking setup and retry on transient Beads connectivity failures:

      bd ready
      bd create "Repo-wide namespace reorganization umbrella"
      bd create "Finish THF scheme namespace completion"
      bd create "Migrate SEB public surface contracts"
      bd create "Extract timestepping namespace"
      bd create "Extract column namespace phase 1"
      bd create "Extract column namespace phase 2"
      bd create "Extract radiation namespace"
      bd create "Extract vapor namespace"
      bd create "Delete compatibility shims and root legacy aliases"
      bd create "Update docs tests and function signatures for namespace reorg"
      bd create "Evaluate diagnostics and measurements namespace design"
      bd update icemodel-lbk --status in_progress
      bd update icemodel-bwq --status in_progress

   If the Dolt circuit breaker trips again, wait for the cooldown and retry rather than skipping the status update.

2. Surface/turbulence migration:

      rg -n "turbulent_heat_flux|SEBSOLVE|SEBFLUX|SFCFLUX|SFCFLIN|SFCTEMP|STABLEFN|LATENT|SENSIBLE|WINDCOEF" icemodel test

   Then move the touched files into the new `surface` and `surface.turbulence` namespaces, update the function names, and update all callers.

3. Timestepping/numerics migration:

      rg -n "check_substep|reset_substep|update_substep|init_timesteps|new_timestep|next_step|aitken_scalar|search_zero|trisolve" icemodel test

   Move the helpers, update tests, and remove the root copies once callers are switched.

4. Focused validation after the touched migrations:

      matlab -batch "addpath(genpath(pwd)); results = runtests({'test/unit/test_turbulent_flux_schemes.m','test/unit/test_surface_solver_kernels.m','test/unit/test_flux_and_thermo_kernels.m','test/unit/test_run_contracts.m','test/unit/test_mesh_and_timestep_kernels.m'}); disp(table({results.Name}',[results.Passed]',[results.Failed]','VariableNames',{'Name','Passed','Failed'})); assert(all([results.Passed]));"

   Expand the suite if column, vapor, or radiation files are moved in the same pass.

## Validation and Acceptance

Acceptance for the active branch-critical slice requires:

1. The surface/turbulence call chain uses the new public API names:
   - `icemodel.surface.diagnose_turbulent_heat_fluxes`
   - `icemodel.surface.diagnose_surface_energy_fluxes`
   - `icemodel.surface.solve_surface_energy_balance`
2. The deep scheme implementations live under:
   - `icemodel.surface.turbulence.bulk_richardson`
   - `icemodel.surface.turbulence.monin_obukhov`
3. The timestep/substep controller helpers live under `icemodel.timestepping`.
4. The generic math helpers moved in this pass live under `icemodel.numerics`.
5. Focused unit tests for the touched areas pass.
6. No permanent compatibility shims or thin wrappers remain for the completed migration slices.

If later phases such as `column`, `radiation`, or `vapor` remain open at a pause point, this plan must record exactly which Beads remain open and which namespaces are already complete.

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

Current status-update blocker:

    bd update icemodel-lbk --status in_progress
    Error resolving icemodel-lbk: dolt circuit breaker is open: server appears down

This blocker should be retried until the umbrella and active slice are marked in progress.

## Interfaces and Dependencies

At the end of the active surface/turbulence migration, the production-facing surface interface must include:

- `icemodel.surface.diagnose_turbulent_heat_fluxes(...)`
- `icemodel.surface.diagnose_surface_energy_fluxes(...)`
- `icemodel.surface.solve_surface_energy_balance(...)`
- `icemodel.surface.surface_energy_balance_residual(...)`
- `icemodel.surface.surface_flux_linearization(...)`

The deep turbulence implementations must live under:

- `icemodel.surface.turbulence.bulk_richardson`
- `icemodel.surface.turbulence.monin_obukhov`

The completed active slices in this task must not retain old root aliases such as `SEBSOLVE`, `SEBFLUX`, `STABLEFN`, `LATENT`, `SENSIBLE`, or `check_substep`.

Revision note: created at task start for Bead `icemodel-lbk` and updated from the final approved user plan to reflect the locked namespace decisions, the Bead breakdown, and the explicit no-permanent-shims rule.
