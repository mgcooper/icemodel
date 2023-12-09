function Qe = LATENT(De, S, ea, es, roL, epsilon, Pa)
   %LATENT Compute latent heat flux
   %
   %  Qe = LATENT(De, S, ea, es, roL, epsilon, Pa)
   %
   %  Qe = ro_air * L * De_h * stability * (0.622/Pa * (ea - es));
   %  [W m-2] = [kg m-3] * [J kg-1] * [m s-1] * [-] * [Pa-1 * Pa]
   %
   % epsilon = 0.622 = Rd/Rv where Rd = 287 J/K/kg is gas constant for dry air
   % and Rv is gas constant for water. Pa is reference pressure, ea atmospheric
   % (2-m) vapor pressure, es surface saturation vapor pressure.
   %
   % See also: SENSIBLE
   Qe = roL * De * S * (epsilon / Pa * (ea - es));
end
