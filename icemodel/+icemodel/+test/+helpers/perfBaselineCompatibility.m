function [compatible, reason] = perfBaselineCompatibility( ...
      baseline_meta, isolation)
   %PERFBASELINECOMPATIBILITY Decide whether wall-time comparison is fair.
   %
   %  [compatible, reason] = ...
   %     icemodel.test.helpers.perfBaselineCompatibility(baseline_meta, ...
   %     isolation)
   %
   % A timing gate is fair only when the baseline and the current run
   % measured under the same environment and the same isolation protocol.
   % When COMPATIBLE is false, the caller reports validity-only results
   % and REASON states the mismatch.
   %
   % See also: icemodel.test.helpers.formalPerformanceVerdict,
   %  run_perf_suite, build_perf_baseline

   compatible = false;
   reason = "";

   if ~isstruct(baseline_meta) || isempty(fieldnames(baseline_meta))
      reason = "perf baseline metadata not found";
      return
   end

   if ~isfield(baseline_meta, 'matlab_version') || ...
         ~isfield(baseline_meta, 'host') || ...
         ~isfield(baseline_meta, 'hostname') || ...
         isblanktext(baseline_meta.matlab_version) || ...
         isblanktext(baseline_meta.host) || ...
         isblanktext(baseline_meta.hostname)
      reason = "perf baseline predates environment metadata";
      return
   end

   % Timings measured under different isolation protocols are not
   % comparable: process-isolated cases pay uniform cold-ish starts,
   % while in-session cases inherit shared warm state (measured here as
   % a systematic ~25 percent offset). A baseline saved before the
   % isolation field exists counts as "session".
   baseline_isolation = "session";
   if isfield(baseline_meta, 'isolation') ...
         && ~isblanktext(baseline_meta.isolation)
      baseline_isolation = string(baseline_meta.isolation);
   end
   if baseline_isolation ~= isolation
      reason = sprintf([ ...
         'baseline measured under isolation="%s"; this run uses ', ...
         'isolation="%s"; timings are not comparable'], ...
         char(baseline_isolation), char(isolation));
      return
   end

   current_version = string(version);
   current_host = string(computer);
   current_hostname = icemodel.test.helpers.machineHostname();
   baseline_version = string(baseline_meta.matlab_version);
   baseline_host = string(baseline_meta.host);
   baseline_hostname = string(baseline_meta.hostname);
   compatible = current_version == baseline_version && ...
      current_host == baseline_host && current_hostname == baseline_hostname;

   if ~compatible
      reason = sprintf([ ...
         'baseline built under MATLAB %s on %s (%s); current environment ', ...
         'is MATLAB %s on %s (%s)'], char(baseline_version), ...
         char(baseline_host), char(baseline_hostname), ...
         char(current_version), char(current_host), char(current_hostname));
   end
end
