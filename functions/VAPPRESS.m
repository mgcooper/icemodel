function ea = VAPPRESS(rh,Tair,liqflag)
   %VAPPRESS compute atmospheric vapor pressure from relative humidity data.
   %
   % Compute atmospheric vapor pressure from relative humidity data using the
   % magnus formula: ea = A*exp(B*T/C+T) and coeffs for saturation vapor
   % pressure over water (Buck 1981). Note: temperatures for Buck's equations
   % are in deg C, and vapor pressures are in mb.  Do the adjustments so that
   % the calculations are done with temperatures in K, and vapor pressures in
   % Pa. (1 mb = 100 Pa)
   %
   % Formula:
   % ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)));
   %
   % Units:
   % [Pa = J m-3 = N m-2 = kg m-1 s-2]
   %
   % See also: VAPOR

   if liqflag == true
      % Over water.
      ea = rh/100.0*611.21*exp((17.502*(Tair-273.16))/(Tair-32.19));
   else
      % Over ice.
      ea = rh/100.0*611.15*exp((22.452*(Tair-273.16))/(Tair-0.61));
   end

   % if Tsfc >= Tf && liqflag == true
   %    % Over water.
   %    A = 611.21;     % reference vapor pressure over water [Pa]
   %    B = 17.502;     % unitless coefficient
   %    C = 240.97;     % unitless coefficient
   % else
   %    % Over ice.
   %    A = 611.15;
   %    B = 22.452;
   %    C = 272.55;
   % end
end
