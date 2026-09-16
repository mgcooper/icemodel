function [failed_cases, failed_gates] = regressionFailures(report, test_passed)
   %REGRESSIONFAILURES List the failed cases and gates of one regression run.
   %
   %  [failed_cases, failed_gates] = ...
   %     icemodel.test.helpers.regressionFailures(report, test_passed)
   %
   % Inputs
   %  report      - saved regression report table with the columns case_id,
   %                passed, and failed_gates.
   %  test_passed - true when every unittest result of the run passed.
   %
   % Outputs
   %  failed_cases - string column of failed case identifiers.
   %  failed_gates - table with one row per failed case and the columns
   %                 case_id and failed_gates. failed_gates holds the
   %                 comma-separated gate names from regressionCaseGates.
   %
   % A unittest failure after the class saved its report, such as a teardown
   % error, belongs to no single report row. Every case then lists the gate
   % "test_framework". A failure before the save never reaches this function:
   % run_regression_suite raises icemodel:test:regressionArtifactMissing.
   %
   % See also: icemodel.test.helpers.regressionCaseGates, run_regression_suite

   % A passing run has no failed case.
   failed_cases = strings(0, 1);
   gates = strings(0, 1);

   % Report the rows whose gates failed. When no row failed, the unittest
   % failure came from outside the case loop and applies to every case.
   if ~test_passed
      failed = ~report.passed;
      if any(failed)
         failed_cases = string(report.case_id(failed));
         gates = string(report.failed_gates(failed));
      else
         failed_cases = string(report.case_id);
         gates = repmat("test_framework", size(failed_cases));
      end
   end

   failed_cases = failed_cases(:);
   failed_gates = table(failed_cases, gates(:), ...
      'VariableNames', {'case_id', 'failed_gates'});
end
