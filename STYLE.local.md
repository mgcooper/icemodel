# Project-specific code style — icemodel

Conventions specific to icemodel, extending the canonical `STYLE.md` (general +
MATLAB). This file is project-owned — `--update` never overwrites it.

## Function naming (override of the canonical MATLAB rule)

The canonical short-lowercase / camelCase-when-longer rule applies, with one
scale-driven exception: in the **core physics namespaces** —
`+icemodel/+column`, `+couplers`, `+surface`, `+vapor`, `+kernels`, `+radiation` —
use `snake_case` for longer descriptive functions, because the project's size makes
long camelCase names unreadable. Short core-physics names may still be all-lowercase
one word if readable.

- Examples: `solve_surface_temperature`, `evaluate_surface_energy_balance`,
  `solvetwostream`.
- Non-core namespaces (helper, verification, I/O, plotting, setup, utility, support)
  follow the canonical rule: all-lowercase for short names, camelCase for longer ones
  (e.g. `snowDataRoot`). Short five-syllable exceptions may remain one word for local
  consistency, e.g. `getsubstepforcings` next to `getforcings`.
- Rationale: variables here are commonly lowercase or snake_case, so camelCase in
  non-core helpers improves visual distinction between functions and variables, e.g.
  `snow_data_root = icemodel.verification.helpers.snowDataRoot()`.

Tests:

- Do not infer global naming rules from the `test/` folder.
- Existing test-suite files may use snake_case because long descriptive test entry points are already deeply embedded.
- The `+test` namespace follows the canonical lowercase/camelCase convention.
- New non-test namespaces should not copy `test/` snake_case unless they are core physics or long descriptive test-suite entry points.

## Variable prefixes

- `f_` for fractions.
- `d_` for deltas/increments.
- `x`-prefixed variables for checkpoint/prior-substep state.
- `cpl_` for coupler options/diagnostics.
- `seb_` for surface-energy-balance options/diagnostics.
- `ok` / `ok_*` for success flags.
- `n_` for counts.
- `use_` for boolean feature toggles.

## Sandbox new physics in production namespaces

When developing a new feature that touches core `icemodel.m` physics,
prefer adding the new kernels under the relevant namespace
(`+icemodel/+column/`, `+icemodel/+vapor/`, etc.) or as opt-in code
paths inside existing functions guarded by `opts` flags that nothing
currently sets. Do not wire the new code into `icemodel.m` until the
physics is verified and the wiring is itself the deliverable.

This keeps the production solver and regression baselines stable
while the repo accumulates physics surface that future features can
target. Verification kernels and analytical references can ship
before the model that uses them. Later wiring becomes a focused PR
of its own rather than a hidden side effect of the physics work.

When in doubt, file the "wire into `icemodel.m`" work as a separate
issue and stop at the namespace-only deliverable.

## Kernel conventions

- Keep core numerical kernels minimalist, without argument blocks or strict input parsing. Enforce contracts at the model/helper boundary instead.
- Argument blocks are acceptable in namespace functions and test-suite code.
- Use reusable validators under `+validators` instead of hard-coded member lists when the choices are part of a stable repo contract.
- Helper layers such as namespaced loaders, postprocessing, config functions, and setup functions may still use old-style parsing when that better matches surrounding code.
- Maintain codegen compatibility in kernel functions.
- Preserve established input/output ordering and naming. For example, use `T`, `f_ice`, `f_liq`, then related state/diagnostics.
- Surface namespace functions follow a preferred ordering schema:
  1. state variables (`T_sfc`, `T_ice`, `f_ice`, `f_liq`)
  2. source terms and grid (`Sc`, `Sp`, `dz`, `delz`, `fn`)
  3. timestep (`dt`)
  4. primary forcings (`tair`, `swd`, `lwd`, `albedo`, `wspd`, etc.)
  5. derived atmospheric variables (`ea_atm`, `ro_atm`, `cv_atm`, `nu_air`)
  6. transport prefactors (`H_h`, `H_e`, `hv_atm`)
  7. surface state (`br_coefs`, `liqflag`, `chi`)
  8. conduction terms (`T_ice`, `k_eff`, `dz`) in SEB-level functions
  9. surface properties (`ro_sfc`, `snow_depth`)
  10. solver options, then `opts` last
- This surface ordering schema is not enforced repo-wide; other namespaces may differ.
- When editing kernels, keep the code shape close to surrounding kernels unless there is a strong reason to refactor more broadly.
