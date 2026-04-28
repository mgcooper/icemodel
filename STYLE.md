# Code Style

## Function naming protocol

 Default:

- Prefer MATLAB-style all-lowercase, single-word function names when the name is short and readable.
- Use this for names of four syllables or fewer.
- Examples: loadcases, plotcase, getforcings, concatoutput.

Core physics namespaces:

- Use snake_case for longer descriptive functions in:
  +icemodel/+column
  +icemodel/+couplers
  +icemodel/+surface
  +icemodel/+vapor
  +icemodel/+kernels
  +icemodel/+radiation
- Short core-physics names may still be all-lowercase one word if readable.
- Examples: solve_surface_temperature, evaluate_surface_energy_balance, solvetwostream.

Non-core-physics namespaces:

- Use all-lowercase one word for short names.
- Use camelCase for longer names where single-word lowercase becomes hard to read.
- This includes helper, verification, I/O, plotting, setup, utility, and support namespaces.
- Short, punchy five-syllable exceptions may remain one word when this preserves local consistency, e.g. `getsubstepforcings` next to `getforcings`.
- Examples: snowDataRoot, loadcases, plotcase.

Tests:

- Do not infer global naming rules from the test/ folder.
- Existing test-suite files may use snake_case because long descriptive test entry points are already deeply embedded.
- The +test namespace follows the preferred lowercase/camelCase convention.
- New non-test namespaces should not copy test/ snake_case unless they are core physics or long descriptive test-suite entry points.

Rationale:

- Variable names are commonly lowercase or snake_case.
- camelCase in non-core helper/utility functions improves visual distinction between functions and variables, e.g.
  snow_data_root = icemodel.verification.helpers.snowDataRoot().

## Variable prefixes

- `f_` for fractions.
- `d_` for deltas/increments.
- `x`-prefixed variables for checkpoint/prior-substep state.
- `cpl_` for coupler options/diagnostics.
- `seb_` for surface-energy-balance options/diagnostics.
- `ok` / `ok_*` for success flags.
- `n_` for counts.
- `use_` for boolean feature toggles.

## Function files and reuse

- Do not repeat logic across functions; if logic is shared, define one helper and call it.
- Define canonical contract/helper methods as their own function files.
- Do not leave shared contract logic as local subfunctions or inline blocks.
- Local subfunctions are acceptable only for truly file-local glue such as small parsing helpers, formatting helpers, and one-off local table/display helpers.
- Before adding a new helper, check the existing `+icemodel` namespace for an equivalent function.
- Organize reusable functions into the relevant namespaces instead of leaving them at arbitrary top level.
- When moving or renaming a function, update all call sites, tests, and docs/comments that describe the public entry point.

## Documentation and comments

- Add docstrings when creating new functions and comment on all code.
- If an existing function is edited materially, update or add its docstring.
- Group operations logically when writing new code.
- Preserve information in existing comments while ensuring accuracy; if a comment is stale, update it rather than deleting it casually.

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
