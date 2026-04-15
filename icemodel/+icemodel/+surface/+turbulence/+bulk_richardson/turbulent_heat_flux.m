function [Qe, Qh, diag] = turbulent_heat_flux(T_sfc, tair, wspd, psfc, ...
      ea_atm, De, br_coefs, ro_air_Lv, liqflag, z0_bulk)
   %TURBULENT_HEAT_FLUX Evaluate the bulk-Richardson THF scheme.
   %
   %  [Qe, Qh] = ...
   %     icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux(...)
   %  [Qe, Qh, diag] = ...
   %     icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux(...)
   %
   % This is the existing bulk-Richardson surface-transfer scheme used by the
   % production SEB path. The formulation uses one aerodynamic roughness, one
   % exchange coefficient, and the Louis/Liston stability factor.
   %
   %#codegen

   % Saturation vapor pressure at the surface.
   es_sfc = icemodel.vapor.saturation_vapor_pressure(T_sfc, liqflag);

   % Bulk-Richardson stability factor.
   stability = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc, tair, wspd, br_coefs);

   % Latent heat flux at the surface.
   Qe = icemodel.surface.turbulence.bulk_richardson.latent_heat_flux( ...
      es_sfc, ea_atm, De, stability, psfc, ro_air_Lv);

   % Sensible heat flux at the surface. Keep as an optional second output for
   % callers like potential_surface_vapor_tendency that only require Qe.
   if nargout > 1
      Qh = icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux( ...
         T_sfc, tair, De, stability);
   end

   % Diagnostics for callers that request it.
   if nargout > 2
      diag = ...
         icemodel.surface.turbulence.bulk_richardson.bulk_richardson_diagnostics( ...
         T_sfc, es_sfc, tair, wspd, psfc, De, ea_atm, stability, ro_air_Lv, z0_bulk);
   end
end
