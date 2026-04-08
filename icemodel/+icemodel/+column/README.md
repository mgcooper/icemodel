# icemodel.column

Purpose: vertical-column physics, state evolution, phase transforms, hydrology, and mesh ownership.

Planned public areas:
- column-state initialization
- enthalpy/temperature solves
- shortwave source-term assembly
- mesh/state/hydrology helpers

Rules:
- own column physics here rather than in root all-caps files
- keep runtime control in `icemodel.timestepping`, not here

Migration status: planned namespace scaffold only; major file moves still open.
Current migrated entry points:
- `icemodel.column.shortwave_source_term`
