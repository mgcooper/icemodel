function assertPerfBuildQuality(quality, ambient_drift_accepted)
   %ASSERTPERFBUILDQUALITY Require measurement quality before a managed build.
   %
   %  icemodel.test.helpers.assertPerfBuildQuality(quality, ...
   %     ambient_drift_accepted)
   %
   % build_perf_baseline accepts a managed rolling file on the five quality
   % conditions of perfMeasurementQuality and on nothing else. The accepted
   % ambient-drift override is the one exception: when the anchor drifted
   % and the caller accepted that drift, the file is written and records
   % ambient_drift_accepted=true, so snapshot_perf_baseline can refuse it as
   % a release source. Every other failed condition refuses the build.
   %
   % Input
   %  quality                 Output of perfMeasurementQuality.
   %  ambient_drift_accepted  True when the caller accepted a drifted anchor.
   %
   % See also: icemodel.test.helpers.perfMeasurementQuality,
   %  icemodel.test.helpers.assertAmbientBaselineAcceptance,
   %  build_perf_baseline

   arguments
      quality (1, 1) struct
      ambient_drift_accepted (1, 1) logical
   end

   if quality.passed
      return
   end
   names = string(fieldnames(quality.conditions));
   failed = names(~cellfun(@(name) quality.conditions.(name), ...
      cellstr(names)));
   % The override covers the ambient condition alone, so remove it from the
   % failed set only when the caller accepted the drift.
   if ambient_drift_accepted
      failed = failed(failed ~= "ambient_stable");
   end
   if isempty(failed)
      return
   end
   error('icemodel:test:perf:buildQuality', ...
      ['The measured candidate fails the measurement quality ', ...
      'conditions (%s): %s. The managed rolling baseline was not ', ...
      'written.'], char(strjoin(failed, ", ")), ...
      char(strjoin(quality.reasons, "; ")))
end
