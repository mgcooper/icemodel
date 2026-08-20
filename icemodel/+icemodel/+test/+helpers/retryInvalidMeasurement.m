function [perf_data, valid, reason, dispersion, attempts] = ...
      retryInvalidMeasurement(measure_fcn)
   %RETRYINVALIDMEASUREMENT Measure once; re-measure once if invalid.
   %
   %  [perf_data, valid, reason, dispersion, attempts] = ...
   %     icemodel.test.helpers.retryInvalidMeasurement(measure_fcn)
   %
   % MEASURE_FCN(attempt) performs one measurement and returns a struct
   % with sample_times and valid fields (the runPerfCase shape). An
   % invalid sample set (framework error or dispersion above the gate)
   % earns one automatic re-measure, because the common cause is one-off
   % interference. A second invalid set is returned as-is with
   % valid=false, so the verdict fails as "measurement invalid" rather
   % than pass or fail on contaminated numbers.
   %
   % A transient subprocess launch failure
   % (icemodel:test:perf:subprocessFailed) also earns one retry; a second
   % launch failure, or any other error, is rethrown.
   %
   % The function-handle seam exists so the retry policy is testable with
   % stub measurements; production callers close over
   % icemodel.test.helpers.measurePerfCase.
   %
   % See also: icemodel.test.helpers.measurePerfCase,
   %  icemodel.test.helpers.perfSampleValidity

   attempts = 0;
   [perf_data, valid, reason, dispersion] = deal([], false, "", nan);
   while attempts < 2 && ~valid
      attempts = attempts + 1;
      try
         perf_data = measure_fcn(attempts);
      catch launch_err
         % A subprocess launch can fail transiently. One retry covers
         % that; a second failure is real and must stop the run.
         if launch_err.identifier ~= "icemodel:test:perf:subprocessFailed" ...
               || attempts >= 2
            rethrow(launch_err)
         end
         continue
      end
      [valid, reason, dispersion] = ...
         icemodel.test.helpers.perfSampleValidity( ...
         perf_data.sample_times, perf_data.valid);
   end
end
