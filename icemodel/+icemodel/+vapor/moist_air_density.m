function ro_atm = moist_air_density(Pa, ea, Ta)
   %MOIST_AIR_DENSITY Return moist-air density from partial pressures.
   %
   % rho = (p_d / (R_d T)) + (p_v / (R_v T))
   %
   % where p_d = p - e is dry-air pressure and p_v = e is vapor pressure.

   %#codegen

   persistent Rd Rv
   if isempty(Rd)
      [Rd, Rv] = icemodel.physicalConstant('Rd', 'Rv');
   end

   ro_atm = (Pa - ea) / (Rd * Ta) + ea / (Rv * Ta);
end
