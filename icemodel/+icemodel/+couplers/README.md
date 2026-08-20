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
- `initialize_solver_diag` returns the fixed-schema observability record
  the two column couplers (`solve_surface_column_dirichlet`,
  `solve_surface_column_robin`) fill and return as their final output.
  `solve_skin_surface_column` returns its plain `n_iters` count instead.

Outputs:
- Both column couplers return `[T_sfc, T_ice, f_ice, f_liq, U_vap, L_vap,
  k_eff, ok_seb, ok_ieb, ok_cpl, diag]`. `U_vap` is the accepted face vapor
  mass flux. `L_vap` is its face donor latent heat;
  `icemodel.column.couple_vapor_step` uses it to route mass to the phase the
  solve's energy carried.

Recovery:
- The Robin coupler owns recovery as default behavior. When the inner solve
  is healthy but the outer loop exhausts its iterations, the coupler reruns
  once from the entry state. The rerun disables acceleration and sets
  relaxation to the conservative cap `cpl_alpha_min`. It self-suppresses
  when the primary policy is already that conservative pair. `diag.cpl_phase` and
  `diag.cpl_recovered` record the path; the driver counts accepted
  recoveries in the per-forcing-step `cpl_recovery_count` channel. Failure
  dumps are sequence-numbered per session, so successive failures do not
  overwrite each other.

Rules:
- own Picard/Aitken and cross-domain convergence logic here
- call surface and column contracts; do not absorb their physics policy
