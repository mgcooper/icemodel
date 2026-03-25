function met = processmet(met, kwargs)
   %PROCESSMET Post-process an icemodel met timetable.
   %
   %  met = icemodel.processmet(met)
   %  met = icemodel.processmet(met, newTimeStep="native")
   %  met = icemodel.processmet(met, newTimeStep="hourly")
   %
   % Default is "hourly"
   %
   % See also: icemodel.loadmet, icemodel.postprocess

   arguments
      met timetable
      kwargs.newTimeStep (1, :) string {mustBeMember(kwargs.newTimeStep, ...
         ["native", "hourly"])} = "hourly"
   end

   [Tf, SB] = icemodel.physicalConstant('Tf', 'SB');
   emiss = icemodel.parameterLookup('emiss');

   if kwargs.newTimeStep == "hourly"
      met = retime(met, 'hourly', 'mean');
   end
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
