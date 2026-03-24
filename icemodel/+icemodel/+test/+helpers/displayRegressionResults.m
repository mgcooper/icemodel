function displayRegressionResults(results)
   %DISPLAYREGRESSIONRESULTS Display compact results from run_regression_suite.
   %
   %  icemodel.test.helpers.displayRegressionResults(results)
   %  icemodel.test.helpers.displayRegressionResults(artifact_file)
   %
   % The struct form displays the results returned by run_regression_suite.
   % The string form loads a saved regression artifact MAT file and
   % reconstructs the results struct for display.
   %
   % The results struct must contain at minimum:
   %   report - table with per-case regression comparison data
   %
   % Displays four baseline-vs-current comparison tables:
   %   1. runoff_final  vs baseline_runoff_final  (runoff_pct_delta)
   %   2. melt_final    vs baseline_melt_final    (melt_pct_delta)
   %   3. runoff_eval   vs baseline_runoff_eval   (runoff_eval_pct_delta)
   %   4. melt_eval     vs baseline_melt_eval     (melt_eval_pct_delta)
   %
   % This function is the runner-level display entry point. For displaying
   % baseline tables or report tables directly, use displayRegressionSummary
   % instead.

   % Support the one-argument file-path form for interactive use.
   if (ischar(results) || isStringScalar(results))
      results = loadFromFile(results);
   end

   report = results.report;

   if isempty(report)
      return
   end

   % Define the four comparison table specifications. Each row maps a
   % current metric, baseline metric, and percent-delta column to compact
   % display labels for the console summary table.
   pairs = {
      'runoff_final', 'baseline_runoff_final', 'runoff_pct_delta', ...
         'runoff', 'baseline', 'pct_delta'
      'melt_final', 'baseline_melt_final', 'melt_pct_delta', ...
         'melt', 'baseline', 'pct_delta'
      'runoff_eval', 'baseline_runoff_eval', 'runoff_eval_pct_delta', ...
         'runoff_eval', 'baseline_eval', 'pct_delta'
      'melt_eval', 'baseline_melt_eval', 'melt_eval_pct_delta', ...
         'melt_eval', 'baseline_eval', 'pct_delta'
      };

   % Display each comparison table, skipping any whose columns are not
   % present in the report (e.g., when the baseline is empty).
   for k = 1:size(pairs, 1)
      current_var = pairs{k, 1};
      baseline_var = pairs{k, 2};
      pct_var = pairs{k, 3};
      current_label = pairs{k, 4};
      baseline_label = pairs{k, 5};
      pct_label = pairs{k, 6};

      vars = string(report.Properties.VariableNames);
      if ~all(ismember({current_var, baseline_var, pct_var}, vars))
         continue
      end

      T = table(report.case_id, report.(current_var), ...
         report.(baseline_var), report.(pct_var), ...
         'VariableNames', ...
         {'case_id', current_label, baseline_label, pct_label});
      disp(T)
      if k < size(pairs, 1)
         fprintf('\n')
      end
   end
end

function results = loadFromFile(filepath)
   %LOADFROMFILE Load a regression artifact and reconstruct a results struct.

   S = load(char(filepath));

   results = struct();

   % Load the regression report table.
   if isfield(S, 'report')
      results.report = S.report;
   elseif isfield(S, 'RegressionBaseline')
      results.report = S.RegressionBaseline;
   else
      results.report = table();
   end

   % Load metadata and case opts when present.
   if isfield(S, 'meta')
      results.meta = S.meta;
   else
      results.meta = struct();
   end
   if isfield(S, 'case_opts')
      results.case_opts = S.case_opts;
   end

   results.artifact_file = string(filepath);
end
