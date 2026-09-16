% +SURFACE
%
%   Contents file for +SURFACE and its subfolders.
%
%   +SURFACE
%   icemodel.surface.advective_heat_flux                                    - Compute heat advected to the surface by rainfall
%   icemodel.surface.apply_surface_vapor_exchange                           - Apply surface vapor energy demand
%   icemodel.surface.atmospheric_pressure_from_elevation                    - Estimate pressure from elevation
%   icemodel.surface.atmospheric_vapor_pressure                             - Relative humidity to atmospheric vapor pressure
%   icemodel.surface.conductive_heat_flux                                   - Conductive heat flux into the surface and derivative
%   icemodel.surface.diagnose_melt_freeze_energy                            - Diagnose surplus/deficit energy relative to Tf
%   icemodel.surface.diagnose_surface_ablation                              - Diagnose cumulative surface ablation terms
%   icemodel.surface.diagnose_surface_energy_balance                        - Diagnose the full surface energy budget
%   icemodel.surface.diagnose_surface_runoff                                - Diagnose cumulative runoff from surface fluxes
%   icemodel.surface.diagnose_turbulent_heat_fluxes                         - Run the configured THF scheme
%   icemodel.surface.dump_turbulent_heat_flux_debug_state                   - Save THF/SEB failure diagnostics
%   icemodel.surface.empirical_incoming_longwave_radiation                  - Estimate downwelling longwave
%   icemodel.surface.evaluate_surface_energy_balance                        - Evaluate the SEB from known flux terms
%   icemodel.surface.incoming_shortwave_radiation                           - Estimate downwelling shortwave radiation
%   icemodel.surface.initialize_surface_forcings                            - Load the meteorological forcing vectors
%   icemodel.surface.initialize_surface_state                               - Precompute forcing-derived surface state vectors
%   icemodel.surface.net_longwave_radiation                                 - Net surface longwave radiation and T_sfc derivative
%   icemodel.surface.net_shortwave_radiation                                - Compute net absorbed shortwave radiation
%   icemodel.surface.numerical_surface_flux                                 - Evaluate the SEB residual and derivative numerically
%   icemodel.surface.outgoing_longwave_radiation                            - Outgoing longwave radiation and T_sfc derivative
%   icemodel.surface.physical_surface_temperature                           - Cap surface temperature at the melting point
%   icemodel.surface.potential_surface_vapor_demand                         - Diagnose top-cell vapor energy demand
%   icemodel.surface.potential_surface_vapor_exchange                       - Partition surface vapor demand
%   README.md
%   icemodel.surface.resolve_forcing_snow_depth                             - Resolve scalar snow-depth for the THF scheme
%   icemodel.surface.solve_surface_energy_balance                           - Solve the nonlinear surface energy balance
%   icemodel.surface.solve_surface_temperature                              - Solve the explicit bulk-Richardson SEB for T_sfc
%   icemodel.surface.step_observation_heights                               - Select scalar observation heights for one step
%   icemodel.surface.surface_bulk_density                                   - Compute the bulk density of the top model layer
%   icemodel.surface.surface_energy_balance_residual                        - Return the SEB residual at T_sfc
%   icemodel.surface.surface_energy_balance_terms                           - Evaluate the SEB term set at T_sfc
%   icemodel.surface.surface_flux_linearization                             - Linearize the non-conductive surface flux
%   icemodel.surface.surface_roughness_length                               - Select the momentum roughness length z0m
%   icemodel.surface.surface_vapor_mass_flux                                - Convert a surface vapor fraction to a mass flux
%   icemodel.surface.terrain_adjusted_shortwave_radiation                   - Estimate terrain-adjusted shortwave
%   icemodel.surface.update_surface_state                                   - Update the surface state at substep entry
%
%   +SURFACE/+TURBULENCE
%   README.md
%
%   +SURFACE/+TURBULENCE/+BULK_RICHARDSON
%   icemodel.surface.turbulence.bulk_richardson.bulk_richardson_diagnostics - Assemble the full bulk-Richardson diagnostics
%   icemodel.surface.turbulence.bulk_richardson.exchange_coefficients       - Compute bulk-Richardson exchange coefficients
%   icemodel.surface.turbulence.bulk_richardson.latent_heat_flux            - Compute the turbulent latent heat flux
%   icemodel.surface.turbulence.bulk_richardson.richardson_number           - Compute the bulk Richardson number
%   icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux          - Compute the turbulent sensible heat flux
%   icemodel.surface.turbulence.bulk_richardson.stability_factor            - Compute the stability function and derivative wrt T_sfc
%   icemodel.surface.turbulence.bulk_richardson.surface_flux_linearization  - Linearize the surface energy balance equation
%   icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux         - Evaluate the bulk-Richardson THF scheme
%
%   +SURFACE/+TURBULENCE/+MONIN_OBUKHOV
%   icemodel.surface.turbulence.monin_obukhov.monin_obukhov_length          - Return the Monin-Obukhov stability length
%   icemodel.surface.turbulence.monin_obukhov.psi_h_paulson                 - Dyer/Paulson unstable scalar profile correction
%   icemodel.surface.turbulence.monin_obukhov.psi_holtslag                  - Holtslag and de Bruin stable profile correction
%   icemodel.surface.turbulence.monin_obukhov.psi_m_paulson                 - Paulson unstable momentum profile correction
%   icemodel.surface.turbulence.monin_obukhov.scalar_roughness_lengths      - Return scalar roughness lengths for bulk-MO
%   icemodel.surface.turbulence.monin_obukhov.stability_corrections         - Return Monin-Obukhov profile corrections
%   icemodel.surface.turbulence.monin_obukhov.surface_flux_linearization    - Linearize the bulk-MO surface flux
%   icemodel.surface.turbulence.monin_obukhov.turbulent_heat_flux           - Evaluate the Monin-Obukhov THF scheme
%
%   updatecontents.m generated this file on 14 Sep 2026 at 20:13:59.
