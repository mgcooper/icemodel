function [d_pevp, pevp, Qe, T_sfc_phys] = ...
      potential_surface_vapor_tendency(T_sfc, tair, wspd, psfc, ea_atm, ...
      ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
      f_ice, f_liq, dt, dz, snow_depth, opts)
   %POTENTIAL_SURFACE_VAPOR_TENDENCY Top-layer potential vapor tendency.
   %
   %  [d_pevp, pevp, Qe, Ts_phys] = ...
   %     icemodel.surface.potential_surface_vapor_tendency(...)
   %
   % This helper is the physical-flux contract used by icemodel mass
   % balance updates. Ts may be an internal solver temperature that exceeds
   % Tf; the turbulent latent-heat flux and derived vapor tendency must use
   % the physical diagnosed surface temperature Ts_phys = min(Ts, Tf).
   %
   % Delegates to icemodel.kernels.potential_surface_vapor_tendency
   %
   % See also: icemodel.surface.apply_surface_vapor_mass_change,
   %  icemodel.kernels.potential_surface_vapor_tendency
   %
   %#codegen

   % Update surface density for the surface turbulent heat flux scheme.
   ro_sfc = icemodel.surface.surface_bulk_density(f_ice, f_liq);

   % Update latent heat flux.
   T_sfc_phys = icemodel.surface.physical_surface_temperature(T_sfc);
   [Qe, ~] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc_phys, tair, wspd, psfc, ea_atm, ro_atm, cv_atm, nu_air, ...
      H_h, H_e, hv_atm, br_coefs, liqflag, ro_sfc, snow_depth, opts);

   % Compute potential vapor tendency and derivative.
   [d_pevp, pevp] = ...
      icemodel.kernels.potential_surface_vapor_tendency(Qe, dt, dz);
end
