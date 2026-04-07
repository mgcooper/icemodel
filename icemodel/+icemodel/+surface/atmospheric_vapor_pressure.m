function ea = atmospheric_vapor_pressure(Tair, rh, liqflag)
   %ATMOSPHERIC_VAPOR_PRESSURE Convert RH to atmospheric vapor pressure.
   %
   %  ea = icemodel.surface.atmospheric_vapor_pressure(Tair, rh, liqflag)
   %
   % This preserves the current repo contract:
   %
   %   ea = icemodel.vapor.vappress(Tair, liqflag) * rh / 100
   %
   % where RH is supplied in percent and ea is returned in Pa.
   %
   %#codegen

   ea = icemodel.vapor.vappress(Tair, liqflag) .* rh ./ 100;
end
