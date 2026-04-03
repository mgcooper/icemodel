function Qe = LATENT(es_sfc, ea_atm, De, stability, psfc, roL)
   %LATENT Compute latent heat flux
   %
   %  Qe = LATENT(es_sfc, ea_atm, De, stability, psfc, roL)
   %
   %  Qe = ro_air * L * De_h * stability * (epsilon/Pa * (ea_atm - es_sfc));
   %  [W m-2] = [kg m-3] * [J kg-1] * [m s-1] * [-] * [Pa-1 * Pa]
   %
   % epsilon = Rd/Rv where Rd = 287 J/K/kg is gas constant for dry air
   % and Rv is gas constant for water. psfc is atmospheric pressure, ea_atm
   % atmospheric (2-m) vapor pressure, es_sfc surface saturation vapor pressure.
   %
   % See also: SENSIBLE
   %
   %#codegen
   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end
   Qe = roL * De * stability * (epsilon / psfc * (ea_atm - es_sfc));
end
