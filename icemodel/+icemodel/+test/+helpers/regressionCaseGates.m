function [case_passed, failed_gates, checks] = regressionCaseGates( ...
      case_id, S, baseline, bid, baseline_tag, tolerance)
   %REGRESSIONCASEGATES Evaluate every regression gate for one formal case.
   %
   %  [case_passed, failed_gates, checks] = ...
   %     icemodel.test.helpers.regressionCaseGates( ...
   %     case_id, S, baseline, bid, baseline_tag, tolerance)
   %
   % Inputs
   %  case_id      - formal case identifier, used in diagnostic text.
   %  S            - scalar struct of case metrics from summarizeIce1Metrics.
   %  baseline     - baseline table from loadBaseline.
   %  bid          - row index of the case in BASELINE, or [] when absent.
   %  baseline_tag - baseline selector for formalRegressionMetricEvidence.
   %  tolerance    - struct of tolerances that IcemodelRegressionTest
   %                 defines: the scalar pair rel_scalar and abs_scalar, and
   %                 the runoff volume pair rel_runoff_m3 and abs_runoff_m3.
   %
   % Outputs
   %  case_passed  - true when the baseline row exists and every gate passes.
   %  failed_gates - string scalar of comma-separated failed gate names, ""
   %                 when every gate passes. "baseline_row" means the baseline
   %                 has no row for the case. "<metric>:evidence" means the
   %                 metric is not comparable. "<metric>" means the compared
   %                 metric is outside its tolerance.
   %  checks       - struct array with one element per evaluated gate and the
   %                 fields gate, passed, and diagnostic. IcemodelRegressionTest
   %                 verifies each element, so the runner output keeps every
   %                 diagnostic.
   %
   % See also: icemodel.test.helpers.formalRegressionMetricEvidence,
   %  icemodel.test.helpers.regressionFailures

   % Size the check list for the row gate plus two gates per metric.
   metric_names = string(fieldnames(S));
   checks = repmat(makeCheck("", true, ""), 1, 1 + 2 * numel(metric_names));

   % A missing baseline row fails the case before any metric comparison.
   checks(1) = makeCheck("baseline_row", ~isempty(bid), ...
      "baseline missing case=" + string(case_id));
   n = 1;

   % Check that each metric is comparable, then compare it with the baseline
   % inside its metric-specific tolerance.
   if ~isempty(bid)
      for imetric = 1:numel(metric_names)
         metric = metric_names(imetric);
         actual = S.(metric);
         [compare_metric, evidence_passed, reason] = ...
            icemodel.test.helpers.formalRegressionMetricEvidence( ...
            actual, baseline, bid, metric, baseline_tag);
         n = n + 1;
         checks(n) = makeCheck(metric + ":evidence", evidence_passed, reason);

         if compare_metric && evidence_passed
            expected = baseline.(metric)(bid);
            tol = metricTolerance(metric, expected, tolerance);
            n = n + 1;
            checks(n) = makeCheck(metric, abs(actual - expected) <= tol, ...
               "baseline mismatch var=" + metric);
         end
      end
   end
   checks = checks(1:n);

   % A case passes only when every evaluated gate passes. The report row
   % stores the failed gate names as one comma-separated string.
   case_passed = all([checks.passed]);
   failed_gates = strjoin( ...
      [strings(1, 0), checks(~[checks.passed]).gate], ",");
end

function check = makeCheck(gate, passed, diagnostic)
   %MAKECHECK Build one gate outcome record.
   check = struct('gate', string(gate), 'passed', logical(passed), ...
      'diagnostic', string(diagnostic));
end

function tol = metricTolerance(varname, expected, tolerance)
   %METRICTOLERANCE Return one metric-specific scalar tolerance.

   % Runoff volumes use a volume tolerance.
   if endsWith(varname, "_m3")
      tol = max(tolerance.abs_runoff_m3, ...
         tolerance.rel_runoff_m3 * abs(expected));
      return
   end

   % Runoff and melt depths use a relative tolerance with an absolute floor.
   if any(startsWith(varname, ["runoff_", "melt_"]))
      tol = max(1e-4, 1e-4 * abs(expected));
      return
   end

   % Iteration statistics allow a half-iteration difference.
   if contains(varname, "numiter")
      tol = 0.5;
      return
   end

   % Non-converged counts allow a difference of one substep.
   if contains(varname, "not_converged")
      tol = 1.0;
      return
   end

   % Closure and goodness-of-fit statistics share one absolute tolerance.
   if startsWith(varname, "closure_") || startsWith(varname, "gof_")
      tol = 5e-3;
      return
   end

   % Every other scalar metric uses the scalar tolerance pair.
   tol = max(tolerance.abs_scalar, tolerance.rel_scalar * abs(expected));
end
