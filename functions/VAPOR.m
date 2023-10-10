function es0 = VAPOR(Tsfc,Tf,liqflag)
   %VAPOR compute the saturation vapor pressure at the surface
   %
   % capital P in Jordan
   % Pv0_sat = saturation water vapor pressure at T = 0oC = 6.1360 mb
   % Pvk_sat = es0 = saturation water vapor pressure wrt phase k
   % note: Pvk_sat = esi for ice i.e., within the ice matrix, here, es0 is
   %
   % Coeffs for saturation vapor pressure over water (Buck 1981).
   %   Note: temperatures for Buck's equations are in deg C, and vapor
   %   pressures are in mb.  Do the adjustments so that the calculations are
   %   done with temperatures in K, and vapor pressures in Pa.

   if liqflag == true   % Over water.
      A = 611.21;
      B = 17.502;
      C = 240.97;
   else                 % Over ice.
      A = 611.15;
      B = 22.452;
      C = 272.55;
   end

   % Compute the saturated water vapor pressure at the surface.
   es0 = A * exp((B * (Tsfc - Tf))/(C + (Tsfc - Tf))); 			% [Pa]
end
