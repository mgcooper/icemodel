function displayRegressionReport(report)
   %DISPLAYREGRESSIONREPORT Display compact regression compare summaries.
   %
   %  icemodel.test.helpers.displayRegressionReport(report)
   %  icemodel.test.helpers.displayRegressionReport(artifact_file)
   %
   % The one-argument table form displays the four baseline-vs-current compare
   % tables (runoff, melt, runoff_eval, melt_eval). The one-argument string
   % form loads the report from a saved regression artifact MAT file.

   if nargin == 1 && (ischar(report) || isStringScalar(report))
      report = loadFromFile(report);
   end

   if isempty(report)
      return
   end

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
         'VariableNames', {'case_id', current_label, baseline_label, pct_label});
      disp(T)
      if k < size(pairs, 1)
         fprintf('\n')
      end
   end
end

function report = loadFromFile(filepath)
   %LOADFROMFILE Load a regression report from a saved artifact MAT file.

   S = load(char(filepath));
   if isfield(S, 'report')
      report = S.report;
   elseif isfield(S, 'RegressionBaseline')
      report = S.RegressionBaseline;
   else
      report = table();
   end
end
