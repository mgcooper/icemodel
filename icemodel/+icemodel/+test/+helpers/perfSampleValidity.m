function [valid, reason, dispersion] = perfSampleValidity( ...
      sample_times, framework_valid, max_dispersion)
   %PERFSAMPLEVALIDITY Decide whether one case's timing samples are usable.
   %
   %  [valid, reason, dispersion] = ...
   %     icemodel.test.helpers.perfSampleValidity(sample_times, ...
   %     framework_valid, max_dispersion)
   %
   % A formal performance verdict must never come from a contaminated
   % sample set. The perf framework's own Valid flag catches measurement
   % errors, but not one-off interference (another process, a paused
   % machine, JIT state left by an earlier case). Interference shows up as
   % dispersion: one sample far above the median. This gate marks such a
   % set invalid so the runner re-measures once, then fails the case as
   % "measurement invalid" rather than returning a phantom verdict.
   %
   % Inputs
   %  sample_times    - measured wall times [s], one per sample
   %  framework_valid - the perf framework's Valid flag for the result
   %  max_dispersion  - largest accepted max/median ratio (default from
   %                    icemodel.test.helpers.perfMeasurementPolicy)
   %
   % Outputs
   %  valid      - true when the sample set supports a formal verdict
   %  reason     - "" when valid; one-sentence cause otherwise
   %  dispersion - max(sample_times) / median(sample_times)
   %
   % See also: icemodel.test.helpers.formalPerformanceVerdict,
   %  icemodel.test.helpers.runPerfCase

   % The default gate comes from the one formal timing policy.
   if nargin < 3
      max_dispersion = ...
         icemodel.test.helpers.perfMeasurementPolicy().max_dispersion;
   end

   % The dispersion ratio is reported even when an earlier check fails, so
   % the artifact always records what the samples looked like.
   dispersion = nan;
   if ~isempty(sample_times)
      dispersion = max(sample_times) / median(sample_times, 'omitnan');
   end

   % The framework flag overrides the later sample checks: a protocol
   % error invalidates clean-looking samples.
   valid = false;
   if ~framework_valid
      reason = "perf framework marked the samples invalid";
      return
   end

   % Timings must be finite and positive to mean anything.
   if isempty(sample_times) || any(~isfinite(sample_times)) ...
         || any(sample_times <= 0)
      reason = "samples contain nonfinite or nonpositive times";
      return
   end

   % One slow outlier marks interference, not code performance.
   if dispersion > max_dispersion
      reason = sprintf( ...
         "sample dispersion max/median = %.3f exceeds %.3f", ...
         dispersion, max_dispersion);
      return
   end

   valid = true;
   reason = "";
end
