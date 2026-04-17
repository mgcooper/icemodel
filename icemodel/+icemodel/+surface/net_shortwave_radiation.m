function Qsn = net_shortwave_radiation(Qsi, albedo, chi)
   %NET_SHORTWAVE_RADIATION Compute net absorbed shortwave radiation.
   %
   %  Qsn = icemodel.surface.net_shortwave_radiation(Qsi, albedo, chi)
   %
   % The sign convention is positive toward the surface:
   %
   %   Qsn = chi * Qsi * (1 - albedo)
   %
   %#codegen

   Qsn = chi .* Qsi .* (1.0 - albedo);
end
