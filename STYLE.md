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
