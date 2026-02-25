function [Tair,Qsi,Qli,albedo,wspd,Pa,De,ea] = LOADMETDATA(met,metstep,liqflag)
   %LOADMETDATA load the met data for this timestep
   %
   %#codegen

   Tair = met.tair(metstep);
   wspd = met.wspd(metstep);
   Qsi = met.swd(metstep);
   Qli = met.lwd(metstep);
   Pa = met.psfc(metstep);
   De = met.De(metstep);
   ea = VAPPRESS(Tair, 273.16, liqflag) * met.rh(metstep) / 100;
   albedo = met.albedo(metstep);
end
