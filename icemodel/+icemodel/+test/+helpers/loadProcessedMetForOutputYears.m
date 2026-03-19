function met = loadProcessedMetForOutputYears(opts, kwargs)
   %LOADPROCESSEDMETFOROUTPUTYEARS Load processed met limited to output years.
   %
   %  met = icemodel.test.helpers.loadProcessedMetForOutputYears(opts)
   %  met = icemodel.test.helpers.loadProcessedMetForOutputYears(opts, newTimeStep="native")

   arguments
      opts struct
      kwargs.newTimeStep (1, :) string {mustBeMember(kwargs.newTimeStep, ...
         ["native", "hourly"])} = "hourly"
   end

   % Reload the saved forcing and trim it to the retained output years.
   met = icemodel.loadmet(opts);
   met = met(ismember(year(met.Time), opts.output_years), :);
   % Process with the requested cadence so regression metrics stay aligned
   % with the saved-output contract under test.
   met = icemodel.processmet(met, newTimeStep=kwargs.newTimeStep);
end
