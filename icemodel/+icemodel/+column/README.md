# icemodel.column

Purpose: vertical-column physics, state evolution, phase transforms,
hydrology, vapor transport, mass-balance bookkeeping, and mesh ownership.

Planned public areas:
- column-state initialization
- enthalpy/temperature solves
- shortwave source-term assembly
- mesh/layer-merge helpers
- liquid transport / infiltration
- vapor transport
- runoff diagnostics
- mass-balance bookkeeping

Rules:
- own column physics here rather than in root all-caps files
- keep runtime control in `icemodel.timestepping`, not here

Migration status: phase-1 state / enthalpy / thermodynamics migration is
complete; phase 2 mass-transfer, mesh, and runoff helpers are active.
Current migrated entry points:
- `icemodel.column.shortwave_source_term`
- `icemodel.column.initialize_column_state`
- `icemodel.column.solve_column_enthalpy`
- `icemodel.column.solve_column_temperature`
- `icemodel.column.assemble_enthalpy_system`
- `icemodel.column.meltzone_transform`
- `icemodel.column.liquid_fraction_derivative`
- `icemodel.column.liquid_fraction_function`
- `icemodel.column.bulk_enthalpy`
- `icemodel.column.updatestate`
- `icemodel.column.bulk_thermal_conductivity`
- `icemodel.column.firn_thermal_conductivity`
- `icemodel.column.budget_surface_mass_balance`
- `icemodel.column.merge_thin_layers`
- `icemodel.column.infiltration`
- `icemodel.column.liquid_flux`
- `icemodel.column.vapor_mass_transfer`
- `icemodel.column.merge_layer_indices`
- `icemodel.column.merge_layers`
- `icemodel.column.enforce_control_volume_balance`
- `icemodel.column.available_liquid_water`
- `icemodel.column.control_volume_mesh`
- `icemodel.column.diagnose_column_runoff`
- `icemodel.column.water_fraction`
