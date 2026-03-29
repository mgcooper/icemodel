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
   for k = 1:numel(results.meta)
      meta = results.meta(k);
      if isfield(meta, 'baseline_compatible') ...
            && ~all(meta.baseline_compatible) ...
            && isfield(meta, 'compare_reason') ...
            && hasText(meta.compare_reason)

         fprintf('Whole-model perf comparison skipped%s: %s\n', ...
            char(resultLabelSuffix(results, k)), ...
            char(compactText(meta.compare_reason)));
      end
   end

   % --- Benchmark summary ---
   % Benchmarks are optional (include_benchmarks=false produces empty
   % benchmark.summary), so this is a legitimate emptiness check.
   for k = 1:numel(results.benchmark)
      benchmark = results.benchmark(k);
      if isfield(benchmark, 'summary') && ~isempty(benchmark.summary)
         fprintf('Benchmark summary%s:\n', ...
            char(resultLabelSuffix(results, k)))
         keep = intersect(["Name", "SampleSize", "Mean", ...
            "StandardDeviation"], ...
            string(benchmark.summary.Properties.VariableNames), ...
            'stable');
         if ~isempty(keep)
            disp(benchmark.summary(:, keep))
         end
      end
   end

   % --- Benchmark comparison ---
   for k = 1:numel(results.benchmark)
      benchmark = results.benchmark(k);
      if ~isfield(benchmark, 'comparison') || isempty(benchmark.comparison)
         continue
      end

      if isfield(benchmark, 'meta') ...
            && isfield(benchmark.meta, 'baseline_compatible') ...
            && all(benchmark.meta.baseline_compatible)
         fprintf('Benchmark comparison%s:\n', ...
            char(resultLabelSuffix(results, k)))
         keep = intersect(["Name", "Mean", "ref_mean", "pct_delta"], ...
            string(benchmark.comparison.Properties.VariableNames), ...
            'stable');
         if ~isempty(keep)
            disp(benchmark.comparison(:, keep))
         end
      elseif isfield(benchmark, 'meta') ...
            && isfield(benchmark.meta, 'compare_reason') ...
            && hasText(benchmark.meta.compare_reason)
         fprintf('Benchmark comparison skipped%s: %s\n', ...
            char(resultLabelSuffix(results, k)), ...
            char(compactText(benchmark.meta.compare_reason)));
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

function suffix = resultLabelSuffix(results, idx)
   %RESULTLABELSUFFIX Build a compact label for multi-result displays.

   suffix = "";
   parts = strings(0, 1);

   if isfield(results, 'meta') && idx <= numel(results.meta)
      meta = results.meta(idx);
      if isfield(meta, 'smbmodel_filter') && hasText(meta.smbmodel_filter)
         parts(end + 1) = string(meta.smbmodel_filter);
      end
      if isfield(meta, 'solver_filter') && hasText(meta.solver_filter)
         parts(end + 1) = "solver " + string(meta.solver_filter);
      end
   end

   if isempty(parts) && isfield(results, 'artifact_file')
      artifact_file = string(results.artifact_file(:));
      if idx <= numel(artifact_file) && hasText(artifact_file(idx))
         [~, name, ext] = fileparts(char(artifact_file(idx)));
         parts = string(name + ext);
      end
   end

   if ~isempty(parts)
      suffix = " (" + strjoin(parts, ", ") + ")";
   end
end

function tf = hasText(value)
   %HASTEXT True when the value contains any non-blank text.

   if isempty(value)
      tf = false;
      return
   end

   value = string(value(:));
   tf = any(~isblanktext(value));
end

function text = compactText(value)
   %COMPACTTEXT Collapse one-or-more text values into one display string.

   value = string(value(:));
   value = value(~isblanktext(value));
   text = strjoin(value, "; ");
end
