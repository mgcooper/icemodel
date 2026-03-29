function rh = VAPPRESS2RH(ea, Ta, liqflag)
   %VAPPRESS2RH Relative humidity from atmospheric vapor pressure.
   %
   %  rh = VAPPRESS2RH(ea, Ta, liqflag) computes relative humidity [%] from
   %  atmospheric vapor pressure ea [Pa] and air temperature Ta [K] using
   %  VAPPRESS as the single source for saturation vapor pressure.
   %
   % See also: VAPPRESS
   %
   %#codegen

   rh = 100 * ea ./ VAPPRESS(Ta, liqflag);
end
