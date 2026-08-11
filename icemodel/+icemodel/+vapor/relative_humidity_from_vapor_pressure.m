function rh = relative_humidity_from_vapor_pressure(e, tair, liqflag)
   %RELATIVE_HUMIDITY_FROM_VAPOR_PRESSURE Relative humidity from vapor pressure.
   %
   %  rh = icemodel.vapor.relative_humidity_from_vapor_pressure(ea, Ta, liqflag)
   %  computes relative humidity [%] from vapor pressure ea [Pa] and air
   %  temperature tair [K]. It calls
   %  icemodel.vapor.saturation_vapor_pressure for the saturation vapor
   %  pressure.
   %
   % See also: icemodel.vapor.saturation_vapor_pressure
   %
   %#codegen

   rh = 100 * e ./ icemodel.vapor.saturation_vapor_pressure(tair, liqflag);
end
