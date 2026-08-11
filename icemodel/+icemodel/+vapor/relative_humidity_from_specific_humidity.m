function rh = relative_humidity_from_specific_humidity(qair, psfc, tair)
   %RELATIVE_HUMIDITY_FROM_SPECIFIC_HUMIDITY Specific humidity to RH [%].
   %
   %  rh = icemodel.vapor.relative_humidity_from_specific_humidity( ...
   %     qair, psfc, tair)
   %
   % This helper converts specific humidity into vapor pressure with the
   % moist-air mixing-ratio identity. It then maps that vapor pressure onto
   % relative humidity with the repo's canonical saturation-vapor routines.

   % Convert specific humidity q into vapor pressure ea.
   ea = icemodel.vapor.vapor_pressure_from_specific_humidity(qair, psfc);

   % Determine whether rh should be computed wrt ice or liquid.
   liqflag = tair >= icemodel.physicalConstant('Tf');

   % Convert vapor pressure to relative humidity.
   rh = icemodel.vapor.relative_humidity_from_vapor_pressure(ea, tair, liqflag);
   rh = min(max(rh, 0), 100);
end
