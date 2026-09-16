% +COLUMN
%
%   Contents file for +COLUMN and its subfolders.
%
%   +COLUMN
%   icemodel.column.accumulate_phase_budget          - Add one substep's phase-change storage increments
%   icemodel.column.accumulate_remesh_budget         - Add one substep's remesh events to the budget
%   icemodel.column.accumulate_vapor_exchange        - Budget one substep's surface vapor exchange
%   icemodel.column.accumulate_vapor_transport       - Budget one substep's interior vapor transport
%   icemodel.column.apply_vapor_transfer             - Apply vapor phase change increments to the column
%   icemodel.column.assemble_enthalpy_system         - Compute the general equation coefficients
%   icemodel.column.assert_max_water                 - Assert that water fraction does not exceed the maximum
%   icemodel.column.available_liquid_water           - Compute available liquid water in one control volume
%   icemodel.column.budget_surface_mass_balance      - Apply and budget the surface mass balance
%   icemodel.column.bulk_density                     - Bulk density of ice, liquid, and air mixture
%   icemodel.column.bulk_enthalpy                    - Compute the solver-state bulk enthalpy [J m-3]
%   icemodel.column.bulk_specific_heat_capacity      - Bulk cp of ice, liquid, and air mixture
%   icemodel.column.bulk_thermal_conductivity        - Compute bulk effective thermal conductivity
%   icemodel.column.control_volume_mesh              - Compute cell edges and nodes for the column mesh
%   icemodel.column.couple_vapor_step                - Apply subsurface vapor transport across cells
%   icemodel.column.diagnose_column_runoff           - Diagnose cumulative runoff from column mass changes
%   icemodel.column.enforce_control_volume_balance   - Enforce the total-volume constraint
%   icemodel.column.finalize_budget_state            - Record the storage end endpoints for one forcing step
%   icemodel.column.firn_thermal_conductivity        - Compute porous ice thermal conductivity
%   icemodel.column.infiltration                     - Snow column liquid mass + cold-content + conduction update
%   icemodel.column.initialize_budget_state          - Return the zeroed budget for one forcing step
%   icemodel.column.initialize_column_state          - Initialize the 1-d ice column state
%   icemodel.column.initialize_remesh_ledger         - Zeroed remesh event ledger
%   icemodel.column.integrate_column_budget          - Return column-integrated mass and enthalpy storage
%   icemodel.column.liquid_flux                      - Compute the liquid water flux between snowpack layers
%   icemodel.column.liquid_fraction_derivative       - Liquid fraction derivative wrt temperature
%   icemodel.column.liquid_fraction_function         - Project state onto the liquid-fraction function
%   icemodel.column.max_liquid_fraction_change       - Largest f_liq increase a control volume takes
%   icemodel.column.meltzone_bounds                  - Return the canonical mushy-zone liquid fraction bounds
%   icemodel.column.meltzone_transform               - Apply the melt-zone temperature-enthalpy transform
%   icemodel.column.merge_layer_indices              - Choose the pair of layers to merge
%   icemodel.column.merge_layers                     - Combine two control volumes conserving state and sources
%   icemodel.column.merge_thin_layers                - Merge layers that fall below the minimum ice fraction
%   icemodel.column.potential_sublimation            - Convert surface vapor demand to ice fraction
%   README.md
%   icemodel.column.residual_water_fraction          - Volumetric residual-water floor for control volumes
%   icemodel.column.residual_water_pore_fraction     - Residual liquid fraction per pore volume
%   icemodel.column.saturated_hydraulic_conductivity - Snow saturated hydraulic conductivity
%   icemodel.column.shortwave_source_term            - Solve the spectral shortwave source term
%   icemodel.column.solve_column_enthalpy            - Solve the column enthalpy balance
%   icemodel.column.solve_column_temperature         - Solve the 1-dimensional column conduction equation
%   icemodel.column.subsurface_linearization_error   - Diagnose the top-node enthalpy
%   icemodel.column.surface_linearization_error      - Diagnose the Robin surface linearization error
%   icemodel.column.update_grain_radius              - Grow thermal grains from the substep vapor exchange
%   icemodel.column.updatestate                      - Update column thermodynamic state variables
%   icemodel.column.vapor_exchange_is_wet            - Decide which phase a cell exchanges vapor with
%   icemodel.column.vapor_transport_terms            - Build the coupled vapor face transport terms
%   icemodel.column.water_fraction                   - Compute the total volumetric water fraction
%
%   updatecontents.m generated this file on 14 Sep 2026 at 20:13:56.
