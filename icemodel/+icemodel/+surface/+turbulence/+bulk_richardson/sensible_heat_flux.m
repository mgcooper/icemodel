function Qh = sensible_heat_flux(T_sfc, tair, De, stability)
   %SENSIBLE_HEAT_FLUX Compute the turbulent sensible heat flux.
   %
   %  Qh = icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux(...)
   %
   % Qh = ro_air * Cp_air * De_h * stability * (Ta - Ts);
   % [W m-2] = [kg m-3] * [J kg-1 K-1] * [m s-1] * [-] * [K]
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.latent_heat_flux
   %
   %#codegen
   persistent cv_air
   if isempty(cv_air)
      cv_air = icemodel.physicalConstant('cv_air');
   end

   Qh = cv_air * De * stability * (tair - T_sfc);
end
