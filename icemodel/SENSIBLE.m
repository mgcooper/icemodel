function Qh = SENSIBLE(De, S, Ta, Ts, cv_air)
   %SENSIBLE Compute the sensible heat flux
   %
   %  Qh = SENSIBLE(De, S, Ta, Ts, cv_air)
   %
   % Qh = ro_air * Cp_air * De_h * stability * (Ta - Ts);
   % [W m-2] = [kg m-3] * [J kg-1 K-1] * [m s-1] * [-] * [K]
   %
   % see also: LATENT
   %
   %#codegen
   Qh = cv_air * De * S * (Ta - Ts);
end
