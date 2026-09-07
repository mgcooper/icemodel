function [Qh, dQh_dT_sfc] = sensible_heat_flux(T_sfc, tair, H_h, stability, ...
      dstability_dT_sfc)
   %SENSIBLE_HEAT_FLUX Compute the turbulent sensible heat flux.
   %
   %  Qh = icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux(...)
   %  [Qh, dQh_dT_sfc] = ...
   %     icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux(...)
   %
   % Qh = H_h * stability * (Ta - Ts);
   % [W m-2] = [W m-2 K-1] * [-] * [K]
   %
   % where H_h = cv_atm * De = ro_atm * cp_air * De is the sensible heat
   % transport coefficient precomputed at model initialization.
   %
   % When the derivative is requested, provide the temperature derivative of
   % the stability factor as the fifth input. The returned derivative is then
   % the full dQh/dT_sfc that the newton solve uses, not the fixed-stability
   % partial that the linearization uses.
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.latent_heat_flux
   %
   %#codegen

   Qh = H_h .* stability .* (tair - T_sfc);

   if nargout > 1
      dQh_dT_sfc = H_h .* (dstability_dT_sfc .* (tair - T_sfc) - stability);
   end
end
