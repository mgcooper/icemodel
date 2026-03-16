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

   met = icemodel.loadmet(opts);
   met = met(ismember(year(met.Time), opts.output_years), :);
   met = icemodel.processmet(met, newTimeStep=kwargs.newTimeStep);
end
