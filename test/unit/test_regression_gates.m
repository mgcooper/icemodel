classdef test_regression_gates < matlab.unittest.TestCase
   %TEST_REGRESSION_GATES Test per-case regression gate outcomes.
   %
   % The tests drive icemodel.test.helpers.regressionCaseGates and
   % icemodel.test.helpers.regressionFailures with synthetic metrics, so no
   % model run is needed. Expected tolerances restate the rules that
   % IcemodelRegressionTest documents for each metric family.

   properties (Constant)
      % Tolerance values from the IcemodelRegressionTest properties.
      Tolerance = struct('rel_scalar', 1e-6, 'abs_scalar', 1e-9, ...
         'rel_runoff_m3', 1e-4, 'abs_runoff_m3', 1.0)
   end

   properties (TestParameter)
      % Each metric family with its baseline value and expected tolerance.
      family = struct( ...
         'volume', struct('metric', "icemodel_eval_m3", 'value', 5e4, ...
         'tol', 5.0), ...
         'depth', struct('metric', "melt_final", 'value', 2.5, ...
         'tol', 2.5e-4), ...
         'depth_floor', struct('metric', "runoff_eval", 'value', 0.2, ...
         'tol', 1e-4), ...
         'iterations', struct('metric', "mean_Tice_numiter", 'value', 3.0, ...
         'tol', 0.5), ...
         'not_converged', struct('metric', "n_not_converged", 'value', 2.0, ...
         'tol', 1.0), ...
         'closure', struct('metric', "closure_seb_rmse", 'value', 0.6, ...
         'tol', 5e-3), ...
         'fit', struct('metric', "gof_tsfc_rmse", 'value', 1.1, ...
         'tol', 5e-3), ...
         'scalar', struct('metric', "albedo_mean", 'value', 0.5, ...
         'tol', 5e-7))
   end

   methods (Test)
      function tolerance_bounds_each_metric_family(testCase, family)
         % A metric inside its tolerance passes and one outside it fails.
         baseline = table("case_a", family.value, ...
            'VariableNames', {'case_id', char(family.metric)});
         S = struct(char(family.metric), family.value + 0.9 * family.tol);
         returned = icemodel.test.helpers.regressionCaseGates( ...
            "case_a", S, baseline, 1, "rolling", testCase.Tolerance);
         testCase.verifyTrue(returned)

         S.(char(family.metric)) = family.value + 1.1 * family.tol;
         [returned, gates] = icemodel.test.helpers.regressionCaseGates( ...
            "case_a", S, baseline, 1, "rolling", testCase.Tolerance);
         expected = family.metric;
         testCase.verifyFalse(returned)
         testCase.verifyEqual(gates, expected)
      end

      function matching_case_passes_every_gate(testCase)
         % A case equal to its baseline row passes one row gate and two
         % gates per metric.
         [baseline, S] = syntheticCase();
         [passed, returned, checks] = ...
            icemodel.test.helpers.regressionCaseGates( ...
            "case_a", S, baseline, 1, "rolling", testCase.Tolerance);
         expected = "";
         testCase.verifyTrue(passed)
         testCase.verifyEqual(returned, expected)
         testCase.verifyNumElements(checks, 1 + 2 * numel(fieldnames(S)))
      end

      function failed_gates_name_evidence_and_tolerance(testCase)
         % A non-finite metric fails its evidence gate and skips its
         % tolerance gate; an out-of-tolerance metric fails by name.
         [baseline, S] = syntheticCase();
         S.runoff_final = NaN;
         S.melt_final = S.melt_final + 1.0;
         [passed, returned, checks] = ...
            icemodel.test.helpers.regressionCaseGates( ...
            "case_a", S, baseline, 1, "rolling", testCase.Tolerance);
         expected = "runoff_final:evidence,melt_final";
         testCase.verifyFalse(passed)
         testCase.verifyEqual(returned, expected)
         testCase.verifyFalse(any([checks.gate] == "runoff_final"))
         mismatch = checks([checks.gate] == "melt_final");
         testCase.verifyEqual(mismatch.diagnostic, ...
            "baseline mismatch var=melt_final")
      end

      function frozen_release_skips_an_absent_metric(testCase)
         % A frozen release may omit a metric added after its release. The
         % metric passes its evidence gate and gets no tolerance gate, while a
         % rolling baseline fails the same absent metric.
         [baseline, S] = syntheticCase();
         S.new_metric = 1.0;
         [passed, returned, checks] = ...
            icemodel.test.helpers.regressionCaseGates( ...
            "case_a", S, baseline, 1, "v1.1", testCase.Tolerance);
         expected = "";
         testCase.verifyTrue(passed)
         testCase.verifyEqual(returned, expected)
         testCase.verifyTrue(any([checks.gate] == "new_metric:evidence"))
         testCase.verifyFalse(any([checks.gate] == "new_metric"))

         [passed, returned] = icemodel.test.helpers.regressionCaseGates( ...
            "case_a", S, baseline, 1, "rolling", testCase.Tolerance);
         expected = "new_metric:evidence";
         testCase.verifyFalse(passed)
         testCase.verifyEqual(returned, expected)
      end

      function missing_baseline_row_fails_the_row_gate(testCase)
         % A case without a baseline row evaluates only the row gate.
         [baseline, S] = syntheticCase();
         [passed, returned, checks] = ...
            icemodel.test.helpers.regressionCaseGates( ...
            "case_c", S, baseline, [], "rolling", testCase.Tolerance);
         expected = "baseline_row";
         testCase.verifyFalse(passed)
         testCase.verifyEqual(returned, expected)
         testCase.verifyEqual(checks.diagnostic, ...
            "baseline missing case=case_c")
      end

      function passing_run_reports_no_failure(testCase)
         % A passing run returns empty failure lists with fixed columns.
         report = syntheticReport([true; true], ["", ""]);
         [returned, gates] = ...
            icemodel.test.helpers.regressionFailures(report, true);
         expected = strings(0, 1);
         testCase.verifyEqual(returned, expected)
         testCase.verifyEqual(height(gates), 0)
         testCase.verifyEqual(string(gates.Properties.VariableNames), ...
            ["case_id", "failed_gates"])
      end

      function every_failed_case_is_reported(testCase)
         % Two failing cases in one run both appear with their gates.
         report = syntheticReport([true; false; false], ...
            ["", "runoff_final:evidence,melt_final", "baseline_row"]);
         [returned, gates] = ...
            icemodel.test.helpers.regressionFailures(report, false);
         expected = ["case_2"; "case_3"];
         testCase.verifyEqual(returned, expected)
         testCase.verifyEqual(gates.case_id, expected)
         testCase.verifyEqual(gates.failed_gates, ...
            ["runoff_final:evidence,melt_final"; "baseline_row"])
      end

      function framework_failure_marks_every_case(testCase)
         % A unittest failure with no failed row applies to every case.
         report = syntheticReport([true; true], ["", ""]);
         [returned, gates] = ...
            icemodel.test.helpers.regressionFailures(report, false);
         expected = ["case_1"; "case_2"];
         testCase.verifyEqual(returned, expected)
         testCase.verifyEqual(gates.failed_gates, ...
            ["test_framework"; "test_framework"])
      end
   end
end

function [baseline, S] = syntheticCase()
   %SYNTHETICCASE Build one baseline row and matching case metrics.
   S = struct('runoff_final', 1.9, 'melt_final', 2.3, ...
      'closure_seb_rmse', 0.5, 'n_not_converged', 0);
   baseline = [table("case_a", 'VariableNames', {'case_id'}), ...
      struct2table(S)];
end

function report = syntheticReport(passed, failed_gates)
   %SYNTHETICREPORT Build a saved-report table with gate outcome columns.
   n = numel(passed);
   report = table("case_" + string((1:n)'), passed(:), ...
      reshape(string(failed_gates), [], 1), ...
      'VariableNames', {'case_id', 'passed', 'failed_gates'});
end
