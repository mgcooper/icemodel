function Qh = SENSIBLE(T_sfc, tair, De, stability)
   %SENSIBLE Compute the sensible heat flux
   %
   %  Qh = SENSIBLE(T_sfc, tair, De, stability)
   %
   % Qh = ro_air * Cp_air * De_h * stability * (Ta - Ts);
   % [W m-2] = [kg m-3] * [J kg-1 K-1] * [m s-1] * [-] * [K]
   %
   % see also: LATENT
   %
   %#codegen
   persistent cv_air
   if isempty(cv_air)
      cv_air = icemodel.physicalConstant('cv_air');
   end

   Qh = cv_air * De * stability * (tair - T_sfc);
end
