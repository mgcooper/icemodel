# icemodel.couplers

Purpose: surface-column coupling workflows and convergence control.

Public entrypoints:
- `solve_surface_column_dirichlet`
- `solve_surface_column_robin`
- `solve_skin_surface_column`

Shared:
- `accelerate_coupler_iterate` accelerates one Picard step on T_sfc. All three
  solvers call it, so they cannot drift apart on how the iterate is advanced.
  It applies Aitken with relaxation as the fallback, then a safeguarded secant
  step when the last two residuals bracket a root.
- `initialize_coupler_history` returns the empty iterate history that accelerator
  expects.

Rules:
- own Picard/Aitken and cross-domain convergence logic here
- call surface and column contracts; do not absorb their physics policy
- keep old `surface_subsurface` names deleted once callers are migrated

Migration status: active coupling names have been migrated to `surface_column_*`.
