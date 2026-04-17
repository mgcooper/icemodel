function [tair, swd, lwd, albedo, wspd, psfc, ea, De, snow_depth] = ...
      getforcings(met, metstep, liqflag, opts)
   %GETFORCINGS load the forcing data for this timestep
   %
   %#codegen

   % Extract the core forcings
   swd = met.swd(metstep);
   lwd = met.lwd(metstep);
   tair = met.tair(metstep);
   wspd = met.wspd(metstep);
   psfc = met.psfc(metstep);
   albedo = met.albedo(metstep);

   % Compute atmospheric vapor pressure
   ea = icemodel.surface.atmospheric_vapor_pressure( ...
      tair, met.rh(metstep), liqflag);

   % Compute the bulk-richardson exchange coefficient
   De = icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, opts.z0_bulk, opts.z_tair, opts.z_wind);

   % Extract snow depth if it exists in the met data
   if isvariable('snow_depth', met)
      snow_depth = met.snow_depth(metstep);
   else
      snow_depth = nan;
   end
end
