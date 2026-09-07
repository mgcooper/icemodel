# icemodel.couplers

Purpose: surface-column couplers and helpers.

Contents:

- `solve_surface_column_dirichlet`
  - Couples the surface energy balance to the column with a Dirichlet surface
    temperature.
- `solve_surface_column_robin`
  - Couples the surface energy balance to the column with a Robin surface
    condition.
- `solve_skin_surface_column`
  - Couples the skin surface temperature to the subsurface column.

Functions used by all couplers:

- `accelerate_coupler_iterate`
  - Accelerates one Picard step on T_sfc. It applies Aitken with relaxation as
    the fallback, then a safeguarded secant step when the last two residuals
    bracket a root. All three solvers call it.
- `initialize_coupler_history`
  - Returns the empty iterate history used by `accelerate_coupler_iterate`.
- `initialize_solver_settings`
  - Returns the settings used by the couplers and timestep controls. These
    settings include `dt_full_step`, `maxsubstep`, and `debug`.
- `initialize_solver_diag`
  - Returns the forcing-step diagnostics record. `diag.substep` is the default
    record for one solve attempt.
- `update_solver_diag`
  - Copies an accepted `diag.substep` into the forcing-step record and counts
    substeps accepted with recovery settings.

Recovery mode:

- Each coupler makes one attempt with its input settings.
  `icemodel.timestepping.checksubstep` handles a failed attempt. An outer-loop
  failure after successful inner solves gets one retry at the same `dt` when
  the input settings do not already match recovery mode. The retry disables
  acceleration and sets `cpl_alpha` to `cpl_recovery_alpha`. A failed recovery
  attempt shortens `dt`. After `maxsubstep` failures, `checksubstep` forces an
  advance.
- `update_solver_diag` adds each accepted recovery substep to
  `cpl_recovery_count`. Failure dumps use a session sequence number.
