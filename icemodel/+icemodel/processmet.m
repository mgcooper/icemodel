function met = processmet(met)
   %PROCESSMET post process an icemodel met file

   [Tf, emiss, SB] = icemodel.physicalConstant('Tf', 'emiss', 'SB');

   met = retime(met, 'hourly', 'mean');
   met = met(~(month(met.Time) == 2 & day(met.Time) == 29), :);

   if ~isvariable('tsfc', met)
      met.tsfc = nan.*met.tair;
   end

   met.swu = met.swd .* met.albedo;
   met.swn = met.swd .* (1-met.albedo);
   met.lwu = emiss * SB * met.tsfc.^4;
   met.lwd = emiss * met.lwd;
   met.lwn = met.lwd - met.lwu;
   met.netr = met.swn + met.lwn;
   met.tsfc = met.tsfc - Tf; % do this after computing lwu etc
   met.tair = met.tair - Tf; % do this after computing lwu etc
end
