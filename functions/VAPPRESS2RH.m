function rh = VAPPRESS2RH(ea,Tair,liqflag)
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
   % See also: VAPOR

   if liqflag == true
      % Over water.
      rh = 100.0.*ea./(611.21.*exp((17.502.*(Tair-273.16))./(Tair-32.19)));
   else
      % Over ice.
      rh = 100.0.*ea./(611.15.*exp((22.452.*(Tair-273.16))./(Tair-0.61)));
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
   %
   % % Atmospheric vapor pressure from relative humidity data.
   %    ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)));

   % [Pa = J m-3 = N m-2 = kg m-1 s-2]
end
