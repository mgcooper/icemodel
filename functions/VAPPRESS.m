function [e_sat, e_air, T_dew] = VAPPRESS(T, Tf, liqflag, rh)
   %VAPPRESS Compute saturation vapor pressure over water or ice.
   %
   %  E_SAT = VAPPRESS(T, TF, LIQFLAG) Computes saturation vapor pressure 
   %  using the Magnus formula:
   %
   %     E_SAT = ES0 * exp(B * (T - Tf) / (C + T - Tf))
   % 
   %  Where A, B, C are coefficients from Buck (1981), and Tf = 273.16 is 
   %  the triple point temperature. Note that temperatures for Buck's 
   %  equations are in deg C, and vapor pressures are in mb. Here they have
   %  been converted so the calculations are done with temperatures in K, 
   %  and vapor pressures in Pa. (1 mb = 100 Pa). Curves ew1 and ei2 are 
   %  adopted for water and ice, respectively (See Buck Eq. 3,8 and Table 2).
   % 
   %  [E_SAT, E_AIR] = VAPPRESS(T, TF, LIQFLAG, RH) Also computes vapor 
   %  pressure from relative humidity: E_AIR = E_SAT * RH / 100.
   %  
   %  [E_SAT, E_AIR, T_DEW] = VAPPRESS(T, TF, LIQFLAG, RH) Also computes the 
   %  dew point temperature in Kelvins.
   % 
   %  Note that A is the reference vapor pressure, es0 [Pa]
   % 
   % Units: [Pa = J m-3 = N m-2 = kg m-1 s-2]
   % 
   % See also: VAPORHEAT
   
   if liqflag == true || T >= Tf
      % Over water
      % a = 611.21
      % b = 17.502
      % c = 240.97 (c - Tf = -32.19)
      e_sat = 611.21 * exp(17.502 * (T - Tf) / (T - 32.19));
   else
      % Over ice
      % a = 611.15
      % b = 22.452
      % c = 272.55 (c - Tf = -0.61)
      e_sat = 611.15 * exp(22.452 * (T - Tf) / (T - 0.61));
   end
   
   % If requested, compute vapor pressure from relative humidity:
   if nargout > 1 && nargin == 4
      e_air = e_sat * rh / 100;

      % To convert from e_sat to T:
      % T = Tf + c * log(e_sat / a) / (b - log(e_sat / a));

      % To compute Tdew, substitute e_air for e_sat:
      if nargout == 3
         if liqflag == true || T >= Tf
            T_dew = 514.13 * log(e_air/611.21) / (17.502 - log(e_air/611.21));
         else
            T_dew = 545.71 * log(e_air/611.15) / (22.452 - log(e_air/611.15));
         end
      end
   end
   
   % Buck suggests ew1, ei2, fw3, fi3. Below are expressions for f, not used
   % here, but for reference:
   % fw3 = 1.0007 + 3.46e-6 * P; % note: P is in mb, need to convert.
   % fi3 = 1.0003 + 4.18e-6 * P;
end
