function [Tair,Qsi,Qli,albedo,wspd,Pa,De,ea,snow_depth] = LOADMETDATA(met,metstep,liqflag)
   %LOADMETDATA load the met data for this timestep
   %
   %#codegen

   Tair = met.tair(metstep);
   wspd = met.wspd(metstep);
   Qsi = met.swd(metstep);
   Qli = met.lwd(metstep);
   Pa = met.psfc(metstep);
   De = met.De(metstep);
   ea = icemodel.surface.atmospheric_vapor_pressure(Tair, met.rh(metstep), liqflag);
   albedo = met.albedo(metstep);
   if isvariable('snow_depth', met)
      snow_depth = met.snow_depth(metstep);
   else
      snow_depth = nan;
   end
end
