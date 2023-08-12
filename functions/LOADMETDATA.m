function [Tair, rh, wspd, Qsi, Qli, Pa, albedo, De, ea] = LOADMETDATA(met, ...
   iter, wcoef, liqflag)
%LOADMETDATA load the met data for this timestep

Tair = met.tair(iter);
wspd = met.wspd(iter);
Qsi = met.swd(iter);
Qli = met.lwd(iter);
Pa = met.psfc(iter);
rh = met.rh(iter);
De = wspd*wcoef;
ea = VAPPRESS(rh, Tair, liqflag);
albedo = met.albedo(iter);
