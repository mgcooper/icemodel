function Qe = LATENT(De, S, ea, es0, roL, epsilon, Pa)
   %LATENT compute latent heat flux
   %
   %  Qe = ro_air * L * De_h * stability * (0.622/Pa * (ea - es0));
   %  [W m-2] = [kg m-3] * [J kg-1] * [m s-1] * [-] * [Pa-1 * Pa]
   %
   % epsilon = 0.622 = Rd/Rv where Rd = 287 J/K/kg is gas constant for dry air
   % and Rv is gas constant for water Pa is reference pressure, ea vapor
   % pressure, es0 saturation vapor pressure
   %
   % See also: SENSIBLE
   Qe = roL * De * S * (epsilon/Pa * (ea - es0));
end
