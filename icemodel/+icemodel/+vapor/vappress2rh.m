function rh = vappress2rh(ea, Ta, liqflag)
   %vappress2rh Relative humidity from atmospheric vapor pressure.
   %
   %  rh = icemodel.vapor.vappress2rh(ea, Ta, liqflag) computes relative humidity [%] from
   %  atmospheric vapor pressure ea [Pa] and air temperature Ta [K] using
   %  icemodel.vapor.vappress as the single source for saturation vapor pressure.
   %
   % See also: icemodel.vapor.vappress
   %
   %#codegen

   rh = 100 * ea ./ icemodel.vapor.vappress(Ta, liqflag);
end
