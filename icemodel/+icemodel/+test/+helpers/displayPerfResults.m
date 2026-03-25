function displayPerfResults(results)
   %DISPLAYPERFRESULTS Display compact performance results from run_perf_suite.
   %
   %  icemodel.test.helpers.displayPerfResults(results)
   %  icemodel.test.helpers.displayPerfResults(artifact_file)
   %
   % The struct form displays the results returned by run_perf_suite. The
   % string form loads a saved perf artifact MAT file and reconstructs the
   % results struct for display.
   %
   % The results struct must contain at minimum:
   %   case_summary - table of per-case timing summaries
   %   benchmark    - struct with .summary, .comparison, .meta
   %   meta         - struct with run metadata
   %
   % This function is the runner-level display entry point. For displaying
   % baseline tables directly, use displayPerfSummary instead.

   % Support the one-argument file-path form for interactive use.
   if (ischar(results) || isStringScalar(results))
      results = loadFromFile(results);
   end

   % --- Case summary ---

   % Select display columns using intersect for forward-compatible column
   % selection: intersect returns only columns present in the table, so new
   % columns added to case_summary in the future do not break display.
   keep = intersect(["case_id", "median_wall_s", "ref_wall_s", ...
      "gate_wall_s", "baseline_compatible", "passed_perf"], ...
      string(results.case_summary.Properties.VariableNames), 'stable');

   % Fall back to a simpler column set when baseline comparison columns are
   % absent (e.g., when displaying a baseline-build artifact that lacks
   % ref_wall_s and gate_wall_s).
   if isempty(keep)
      keep = intersect(["case_id", "median_wall_s", "n_runs", ...
         "n_warmups"], ...
         string(results.case_summary.Properties.VariableNames), 'stable');
   end
   if ~isempty(keep)
      disp(results.case_summary(:, keep))
   end

   % Report if whole-model perf comparison was skipped due to environment
   % mismatch (e.g., different MATLAB version or host platform).
   if isfield(results.meta, 'baseline_compatible') ...
         && ~results.meta.baseline_compatible ...
         && isfield(results.meta, 'compare_reason') ...
         && ~isblanktext(results.meta.compare_reason)
      fprintf('Whole-model perf comparison skipped: %s\n', ...
         char(results.meta.compare_reason));
   end

   % --- Benchmark summary ---
   % Benchmarks are optional (include_benchmarks=false produces empty
   % benchmark.summary), so this is a legitimate emptiness check.
   if ~isempty(results.benchmark.summary)
      fprintf('Benchmark summary:\n')
      keep = intersect(["Name", "SampleSize", "Mean", ...
         "StandardDeviation"], ...
         string(results.benchmark.summary.Properties.VariableNames), ...
         'stable');
      if ~isempty(keep)
         disp(results.benchmark.summary(:, keep))
      end
   end

   % --- Benchmark comparison ---
   if ~isempty(results.benchmark.comparison)
      if isfield(results.benchmark.meta, 'baseline_compatible') ...
            && results.benchmark.meta.baseline_compatible
         fprintf('Benchmark comparison:\n')
         keep = intersect(["Name", "Mean", "ref_mean", "pct_delta"], ...
            string(results.benchmark.comparison.Properties.VariableNames), ...
            'stable');
         if ~isempty(keep)
            disp(results.benchmark.comparison(:, keep))
         end
      elseif isfield(results.benchmark.meta, 'compare_reason') ...
            && ~isblanktext(results.benchmark.meta.compare_reason)
         fprintf('Benchmark comparison skipped: %s\n', ...
            char(results.benchmark.meta.compare_reason));
      end
   end
end

function results = loadFromFile(filepath)
   %LOADFROMFILE Load a perf artifact and reconstruct a results struct.
   %
   % Handles both the new format (benchmark saved as a single struct) and
   % the old format (benchmark_summary, benchmark_comparison, benchmark_meta
   % saved as separate variables).

   S = load(char(filepath));

   results = struct();

   % Load case_summary, falling back to PerfBaseline for baseline files.
   if isfield(S, 'case_summary')
      results.case_summary = S.case_summary;
   elseif isfield(S, 'PerfBaseline')
      results.case_summary = S.PerfBaseline;
   else
      results.case_summary = table();
   end

   % Load benchmark: prefer the struct form, fall back to flat vars.
   if isfield(S, 'benchmark')
      results.benchmark = S.benchmark;
   else
      results.benchmark = struct('summary', table(), ...
         'comparison', table(), 'meta', struct());
      if isfield(S, 'benchmark_summary')
         results.benchmark.summary = S.benchmark_summary;
      elseif isfield(S, 'BenchmarkBaseline')
         results.benchmark.summary = S.BenchmarkBaseline;
      end
      if isfield(S, 'benchmark_comparison')
         results.benchmark.comparison = S.benchmark_comparison;
      end
      if isfield(S, 'benchmark_meta')
         results.benchmark.meta = S.benchmark_meta;
      end
   end

   % Load metadata.
   if isfield(S, 'meta')
      results.meta = S.meta;
   else
      results.meta = struct();
   end

   % Load remaining fields when present in the artifact.
   if isfield(S, 'sample_detail')
      results.sample_detail = S.sample_detail;
   end
   if isfield(S, 'activity_detail')
      results.activity_detail = S.activity_detail;
   end
   if isfield(S, 'case_opts')
      results.case_opts = S.case_opts;
   end

   results.artifact_file = string(filepath);
end
