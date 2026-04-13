function [Qh, dQh_dT_sfc] = sensible_heat_flux(T_sfc, tair, De, stability, ...
      dstability_dT_sfc)
   %SENSIBLE_HEAT_FLUX Compute the turbulent sensible heat flux.
   %
   %  Qh = icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux(...)
   %  [Qh, dQh_dT_sfc] = ...
   %     icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux(...)
   %
   % Qh = ro_air * Cp_air * De_h * stability * (Ta - Ts);
   % [W m-2] = [kg m-3] * [J kg-1 K-1] * [m s-1] * [-] * [K]
   %
   % When the derivative is requested, provide the temperature derivative of
   % the stability factor as the fifth input so the returned derivative is
   % the full dQh/dT_sfc rather than only the fixed-stability partial.
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.latent_heat_flux
   %
   %#codegen
   persistent cv_air
   if isempty(cv_air)
      cv_air = icemodel.physicalConstant('cv_air');
   end

   Qh = cv_air * (De .* stability .* (tair - T_sfc));

   if nargout > 1
      dQh_dT_sfc = cv_air * ...
         (De .* ((tair - T_sfc) .* dstability_dT_sfc - stability));
   end
end
