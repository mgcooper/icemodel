function displayPerfSummary(kwargs)
   %DISPLAYPERFSUMMARY Display a compact perf case summary and benchmark table.
   %
   %  icemodel.test.helpers.displayPerfSummary(case_summary, benchmark)
   %  icemodel.test.helpers.displayPerfSummary(artifact_file)
   %
   % The two-argument form displays in-memory tables from run_perf_suite or
   % build_perf_baseline. The one-argument form loads a saved perf artifact
   % or baseline MAT file and displays whatever it contains.
   %
   % If run_perf_suite is called with include_benchmarks=false, then the
   % benchmark struct will have an empty benchmark.summary table and results are
   % not reported.

   % Use an arguments block that allows
   arguments(Input)
      kwargs.results struct = struct.empty()
      kwargs.artifacts_file string = string.empty()
   end

   % --- Put all of this in a helper ---

   % Refactored so the output of run_perf_suite is passed in
   % TODO: unify parsing with loadFromFile logic
   [results, artifacts_file] = deal(kwargs.results, kwargs.artifacts_file);

   if isempty(results) && ~isempty(artifacts_file)
      [case_summary, benchmark] = loadFromFile(artifacts_file);

   elseif ~isempty(results) && isempty(artifacts_file)
      case_summary = results.case_summary;
      benchmark = results.benchmark;
   end
   has_case_summary = ~isempty(case_summary);
   has_benchmark_summary = ~isempty(benchmark.summary);
   has_benchmark_comparison = ~isempty(benchmark.comparison);
   % ------

   % Case summary
   if has_case_summary
      keep = intersect(["case_id", "median_wall_s", "ref_wall_s", ...
         "gate_wall_s", "baseline_compatible", "passed_perf"], ...
         string(case_summary.Properties.VariableNames), 'stable');
      if isempty(keep)
         keep = intersect(["case_id", "median_wall_s", "n_runs", ...
            "n_warmups"], ...
            string(case_summary.Properties.VariableNames), 'stable');
      end
      if ~isempty(keep)
         disp(case_summary(:, keep))
      end
   end

   % Benchmark summary
   if has_benchmark_summary
      fprintf('Benchmark summary:\n')
      keep = intersect(["Name", "SampleSize", "Mean", ...
         "StandardDeviation"], ...
         string(benchmark.summary.Properties.VariableNames), 'stable');
      if ~isempty(keep)
         disp(benchmark.summary(:, keep))
      end
   end

   % Benchmark comparison
   if has_benchmark_comparison
      if isfield(benchmark.meta, 'baseline_compatible') ...
            && benchmark.meta.baseline_compatible
         fprintf('Benchmark comparison:\n')
         keep = intersect(["Name", "Mean", "ref_mean", "pct_delta"], ...
            string(benchmark.comparison.Properties.VariableNames), 'stable');
         if ~isempty(keep)
            disp(benchmark.comparison(:, keep))
         end
      elseif isfield(benchmark.meta, 'compare_reason') ...
            && ~isblanktext(benchmark.meta.compare_reason)
         fprintf('Benchmark comparison skipped: %s\n', ...
            char(benchmark.meta.compare_reason));
      end
   end
end

function [case_summary, benchmark] = loadFromFile(filepath)
   %LOADFROMFILE Load perf data from an artifact or baseline MAT file.

   S = load(char(filepath));
   if isfield(S, 'case_summary')
      case_summary = S.case_summary;
   elseif isfield(S, 'PerfBaseline')
      case_summary = S.PerfBaseline;
   else
      case_summary = table();
   end

   benchmark = struct('summary', table(), 'comparison', table(), ...
      'meta', struct());
   if isfield(S, 'benchmark_summary')
      benchmark.summary = S.benchmark_summary;
   elseif isfield(S, 'BenchmarkBaseline')
      benchmark.summary = S.BenchmarkBaseline;
   end
   if isfield(S, 'benchmark_comparison')
      benchmark.comparison = S.benchmark_comparison;
   end
   if isfield(S, 'benchmark_meta')
      benchmark.meta = S.benchmark_meta;
   end
end
