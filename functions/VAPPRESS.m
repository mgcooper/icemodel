function [e_sat, e_air] = VAPPRESS(T, Tf, liqflag, rh)
   %VAPPRESS Compute saturation vapor pressure over water or ice.
   %
   %  E_SAT = VAPPRESS(T, TF, LIQFLAG) Computes saturation vapor pressure using
   %  the Magnus formula:
   %
   %     E_SAT = A * exp(B * (T - Tf) / (C + T - Tf))
   % 
   %  Where A, B, C are coefficients from Buck (1981), and Tf = 273.16 is the
   %  triple point temperature. Note that temperatures for Buck's equations are
   %  in deg C, and vapor pressures are in mb. Here they have been converted so
   %  the calculations are done with temperatures in K, and vapor pressures in
   %  Pa. (1 mb = 100 Pa). Curves ew1 and ei2 are adopted for water and ice,
   %  respectively (See Buck Eq. 3, Table 2, and Eq. 8).
   % 
   %  [E_SAT, E_AIR] = VAPPRESS(T, TF, LIQFLAG, RH) Also computes vapor pressure
   %  from relative humidity: E_AIR = E_SAT * RH / 100.
   % 
   % Units:
   % [Pa = J m-3 = N m-2 = kg m-1 s-2]
   % 
   % See also: VAPORHEAT
   
   if liqflag == true || T >= Tf
      % Over water:
      % A = 100 * 6.1121; % reference vapor pressure [Pa]
      % B = 17.502;
      % C = 240.97;
      e_sat = 611.21 * exp(17.502 * (T - Tf) / (240.97 + T - Tf));
   else
      % Over ice:
      % A = 100 * 6.1115; % reference vapor pressure [Pa]
      % B = 22.452;
      % C = 272.55;
      e_sat = 611.15 * exp(22.452 * (T - Tf) / (272.55 + T - Tf));
   end
   
   % If requested, compute vapor pressure from relative humidity:
   if nargout == 2 && nargin == 4
      e_air = e_sat * rh / 100;
   end
   
   % To convert from e_sat to T:
   % T = Tf + C * log(e_sat / A) / (B - log(e_sat / A));
   
   % Buck suggests ew1, ei2, fw3, fi3. Below are expressions for f, not used
   % here, but for reference:
   % fw3 = 1.0007 + 3.46e-6 * P; % note: P is in mb, need to convert.
   % fi3 = 1.0003 + 4.18e-6 * P;
end
