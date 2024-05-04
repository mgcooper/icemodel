function [Tair,Qsi,Qli,albedo,wspd,Pa,De,ea] = LOADMETDATA(met,iter,liqflag)
   %LOADMETDATA load the met data for this timestep
   %
   %#codegen

   Tair = met.tair(iter);
   wspd = met.wspd(iter);
   Qsi = met.swd(iter);
   Qli = met.lwd(iter);
   Pa = met.psfc(iter);
   De = met.De(iter);
   ea = VAPPRESS(Tair, 273.16, liqflag) * met.rh(iter) / 100;
   albedo = met.albedo(iter);
end
