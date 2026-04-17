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
   % transport prefactor precomputed at model initialization.
   %
   % When the derivative is requested, provide the temperature derivative of
   % the stability factor as the fifth input so the returned derivative is
   % the full dQh/dT_sfc used in the newton solve rather than only the
   % fixed-stability partial used in the linearization.
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.latent_heat_flux
   %
   %#codegen

   Qh = H_h .* stability .* (tair - T_sfc);

   if nargout > 1
      dQh_dT_sfc = H_h .* (dstability_dT_sfc .* (tair - T_sfc) - stability);
   end
end
