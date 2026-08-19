# icemodel.couplers

Purpose: surface-column coupling workflows and convergence control.

Public entrypoints:
- `solve_surface_column_dirichlet`
- `solve_surface_column_robin`
- `solve_skin_surface_column`

Shared:
- `accelerate_coupler_iterate` accelerates one Picard step on T_sfc, and all
  three solvers call it. It applies Aitken with relaxation as the fallback,
  then a safeguarded secant step when the last two residuals bracket a root.
- `initialize_coupler_history` returns the empty iterate history that the
  accelerator expects.

Production recovery:
- `icemodel` uses the configured coupling relaxation and acceleration for the
  primary attempt. Solver 3 alone retries a healthy-inner Robin outer failure
  from the exact prognostic checkpoint with acceleration disabled and
  relaxation no larger than the central conservative cap. A successful retry
  latches that mode for the rest of the run. Diagnostic output records each
  successful latch in the per-forcing-step `cpl_recovery_count` channel.
  Failure dumps remain useful for the latest failure, but a later failure can
  replace an earlier dump.

Rules:
- own Picard/Aitken and cross-domain convergence logic here
- call surface and column contracts; do not absorb their physics policy
