# icemodel.couplers

Purpose: surface-column coupling workflows and convergence control.

Public entrypoints:
- `solve_surface_column_dirichlet`
- `solve_surface_column_robin`

Rules:
- own Picard/Aitken and cross-domain convergence logic here
- call surface and column contracts; do not absorb their physics policy
- keep old `surface_subsurface` names deleted once callers are migrated

Migration status: active coupling names have been migrated to `surface_column_*`.
