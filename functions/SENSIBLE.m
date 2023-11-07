function Qh = SENSIBLE(De, S, Tair, Tsfc, cv_air)
   %SENSIBLE compute the sensible heat flux
   %
   % Qh = ro_air * Cp_air * De_h * stability * (Tair - Tsfc);
   % [W m-2] = [kg m-3] * [J kg-1 K-1] * [m s-1] * [-] * [K]
   %
   % see also: LATENT
   Qh = cv_air * De * S * (Tair - Tsfc);
end
