function [passed, ref_wall, floor_wall, gate_wall, reason] = ...
      formalPerformanceVerdict(valid, current_wall, baseline, rows, ...
      baseline_compatible, tol_perf, compatibility_reason)
   %FORMALPERFORMANCEVERDICT Evaluate one formal timing comparison row.
   %
   % Compatible baselines fail closed unless exactly one finite positive
   % reference row exists. Incompatible environments return a validity-only
   % result because their timings are not comparable.

   [ref_wall, floor_wall, gate_wall] = deal(nan);
   reason = string(compatibility_reason);
   if ~baseline_compatible
      passed = logical(valid);
      if ~passed
         reason = "performance samples are invalid";
      end
      return
   end

   % A compatible comparison needs one unique row and valid current samples.
   passed = false;
   if ~valid
      reason = "performance samples are invalid";
      return
   end
   if ~isscalar(current_wall) || ~isfinite(current_wall) || current_wall <= 0
      reason = "performance sample median is not finite and positive";
      return
   end
   if numel(rows) ~= 1
      reason = "case is missing or duplicated in compatible perf baseline";
      return
   end
   if ~ismember('median_wall_s', baseline.Properties.VariableNames)
      reason = "compatible perf baseline lacks median_wall_s";
      return
   end
   ref_wall = baseline.median_wall_s(rows);
   if ~isscalar(ref_wall) || ~isfinite(ref_wall) || ref_wall <= 0
      reason = "compatible perf baseline reference is not finite and positive";
      return
   end

   % A saved positive case tolerance overrides the runner default.
   tol_case = tol_perf;
   if ismember('tol_perf', baseline.Properties.VariableNames) ...
         && isfinite(baseline.tol_perf(rows)) ...
         && baseline.tol_perf(rows) > 0
      tol_case = baseline.tol_perf(rows);
   end
   [passed, floor_wall, gate_wall, gate_reason] = ...
      icemodel.test.helpers.performanceGate( ...
      current_wall, ref_wall, tol_case);
   if ~passed
      reason = gate_reason;
   end
end
