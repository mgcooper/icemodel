# icemodel.helpers

Purpose: shared runtime utilities used by two or more namespaces.

Contents:

- `copyFields`
  - Copies one scalar struct's fields onto another. By default, it adds all
    source fields. Set `only_matching = true` to copy only fields already in
    the target. `icemodel.couplers.update_solver_diag` uses this option to
    promote the accepted `diag.substep` record.
- `ensureDirExists`
  - Creates a directory when needed.
- `rmttleapinds`
  - Removes leap-day rows from a timetable.
