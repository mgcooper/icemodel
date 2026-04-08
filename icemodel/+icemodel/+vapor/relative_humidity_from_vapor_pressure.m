function rh = relative_humidity_from_vapor_pressure(e, tair, liqflag)
   %relative_humidity_from_vapor_pressure Relative humidity from vapor pressure.
   %
   %  rh = icemodel.vapor.relative_humidity_from_vapor_pressure(ea, Ta, liqflag)
   %  computes relative humidity [%] from vapor pressure ea [Pa] and air
   %  temperature tair [K] using icemodel.vapor.saturation_vapor_pressure as the
   %  single source for saturation vapor pressure.
   %
   % See also: icemodel.vapor.saturation_vapor_pressure
   %
   %#codegen

   rh = 100 * e ./ icemodel.vapor.saturation_vapor_pressure(tair, liqflag);
end
