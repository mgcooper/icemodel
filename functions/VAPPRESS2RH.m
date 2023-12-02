function rh = VAPPRESS2RH(ea, Ta, liqflag)
   %VAPPRESS2RH Relative humidity from atmospheric vapor pressure
   %
   %  Coeffs for saturation vapor pressure over water (Buck 1981).
   %  Note: temperatures for Buck's equations are in deg C, and
   %  vapor pressures are in mb.  Do the adjustments so that the
   %  calculations are done with temperatures in K, and vapor
   %  pressures in Pa. (1 mb = 100 Pa)
   %
   %  Magnus formula: ea = A*exp(B*T/C+T)
   %
   % See also: VAPPRESS

   if liqflag == true
      % Over water.
      rh = 100.0.*ea./(611.21.*exp((17.502.*(Ta-273.16))./(Ta-32.19)));
   else
      % Over ice.
      rh = 100.0.*ea./(611.15.*exp((22.452.*(Ta-273.16))./(Ta-0.61)));
   end
end
