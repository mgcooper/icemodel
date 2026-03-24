function benchmark = runBenchmarkDiagnostics(simyear, baseline_tag, smbmodel, ...
      include_benchmarks, benchmark_sampling_profile)
   %RUNBENCHMARKDIAGNOSTICS Run and compare managed component benchmarks.
   %
   %  benchmark = icemodel.test.helpers.runBenchmarkDiagnostics( ...
   %     simyear, baseline_tag, smbmodel, include_benchmarks, ...
   %     benchmark_sampling_profile)
   %
   % Runs the managed benchmark suite, loads the embedded benchmark baseline
   % from the perf baseline file, and returns a struct with fields:
   %   summary    - table of current benchmark timings
   %   comparison - table joining current to baseline with pct_delta
   %   meta       - struct with compatibility flags and provenance

   % Initialize an empty benchmark struct for early return when benchmarks
   % are disabled.
   benchmark = struct();
   benchmark.summary = table();
   benchmark.comparison = table();
   benchmark.meta = struct();

   % Short-circuit when the caller disabled component benchmarks.
   if ~include_benchmarks
      return
   end

   % Run the managed benchmark suite with summary display suppressed
   % (display is handled by displayPerfResults at the runner level).
   bench_results = run_benchmark_suite( ...
      sampling_profile=benchmark_sampling_profile, show_summary=false);

   % Compute the current suite content hash for baseline compatibility.
   [current_signature, suite_files] = ...
      icemodel.test.helpers.benchmarkSuiteSignature();

   % Build the summary table from the benchmark MeasurementResult array.
   benchmark.summary = sampleSummary(bench_results);

   % Coerce benchmark names to string and ensure the Name column is
   % consistent across MATLAB versions.
   if ~isempty(benchmark.summary) ...
         && ismember('Name', benchmark.summary.Properties.VariableNames)
      benchmark.summary.Name = string(benchmark.summary.Name);
   end

   % Carry the Valid flag from each benchmark result for artifact
   % completeness.
   benchmark.summary.Valid = reshape(logical([bench_results.Valid]), [], 1);

   % Load the embedded benchmark baseline from the matching perf baseline.
   [baseline_summary, baseline_meta, source_file] = ...
      loadBenchmarkBaselineFromPerf(simyear, baseline_tag, smbmodel);

   % Compare current benchmark timings against the loaded baseline.
   benchmark.comparison = compareBenchmarkSummary( ...
      benchmark.summary, baseline_summary, baseline_meta, current_signature);

   % Populate benchmark metadata: sampling profile, suite signatures,
   % baseline provenance, and compatibility status.
   benchmark.meta = struct();
   benchmark.meta.sampling_profile = benchmark_sampling_profile;
   benchmark.meta.current_signature = current_signature;
   benchmark.meta.suite_files = suite_files;
   benchmark.meta.baseline_file = source_file;
   benchmark.meta.baseline_meta = baseline_meta;
   if isfield(baseline_meta, 'suite_signature')
      benchmark.meta.baseline_signature = string(baseline_meta.suite_signature);
   else
      benchmark.meta.baseline_signature = "";
   end

   % Determine baseline compatibility based on suite signature matching.
   if isempty(source_file)
      benchmark.meta.baseline_compatible = false;
      benchmark.meta.compare_reason = "no embedded benchmark baseline found";
   elseif benchmark.meta.baseline_signature == ""
      benchmark.meta.baseline_compatible = false;
      benchmark.meta.compare_reason = ...
         "embedded benchmark baseline predates suite signatures";
   else
      benchmark.meta.baseline_compatible = ...
         benchmark.meta.current_signature == benchmark.meta.baseline_signature;
      if benchmark.meta.baseline_compatible
         benchmark.meta.compare_reason = "";
      else
         benchmark.meta.compare_reason = ...
            "embedded benchmark baseline was built from a different benchmark suite";
      end
   end
   benchmark.meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
end

function [BenchmarkBaseline, meta, source_file] = loadBenchmarkBaselineFromPerf( ...
      simyear, baseline_tag, smbmodel)
   %LOADBENCHMARKBASELINEFROMPERF Load managed benchmark timing from perf files.

   % Resolve the perf baseline selector once before probing model files.
   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_tag);

   % Expand the model selector to probe one or more baseline files.
   if smbmodel == "all"
      models = icemodel.namelists.smbmodel("test");
   else
      models = string(smbmodel);
   end

   % Initialize empty returns for the case where no matching baseline is
   % found.
   BenchmarkBaseline = table();
   meta = struct();
   source_file = "";

   % Probe each candidate baseline file until one carries the managed
   % benchmark bundle.
   for i = 1:numel(models)
      candidate = icemodel.test.helpers.defaultBaselinePath( ...
         "perf", baseline_type, baseline_tag, models(i), simyear);
      if exist(char(candidate), 'file') ~= 2
         continue
      end

      % Probe the MAT variables first so older perf baselines that predate
      % benchmark bundles do not emit missing-variable warnings.
      matvars = string({whos('-file', char(candidate)).name});
      has_baseline = any(matvars == "BenchmarkBaseline");
      has_meta = any(matvars == "benchmark_meta");
      if ~has_baseline && ~has_meta
         continue
      end

      % Load only the variables that actually exist in this baseline file.
      if has_baseline && has_meta
         S = load(char(candidate), 'BenchmarkBaseline', 'benchmark_meta');
      elseif has_baseline
         S = load(char(candidate), 'BenchmarkBaseline');
      else
         S = load(char(candidate), 'benchmark_meta');
      end

      if has_baseline
         BenchmarkBaseline = S.BenchmarkBaseline;
      end
      if has_meta
         meta = S.benchmark_meta;
      end
      source_file = string(candidate);
      return
   end
end

function comparison = compareBenchmarkSummary(current_summary, baseline_summary, ...
      baseline_meta, current_signature)
   %COMPAREBENCHMARKSUMMARY Join current benchmark timings to a baseline.

   % Early return for empty current summary.
   if isempty(current_summary)
      comparison = table();
      return
   end

   % Retain only the timing columns meaningful for benchmark-to-baseline
   % comparison. intersect selects only columns present in the table for
   % forward-compatible column selection.
   comparison = current_summary(:, intersect(["Name", "SampleSize", "Mean", ...
      "StandardDeviation"], string(current_summary.Properties.VariableNames), ...
      'stable'));
   if ~ismember('Name', comparison.Properties.VariableNames)
      return
   end

   % Initialize baseline comparison columns with default NaN/false values.
   n_rows = height(comparison);
   comparison.baseline_compatible = false(n_rows, 1);
   comparison.ref_mean = nan(n_rows, 1);
   comparison.ref_std = nan(n_rows, 1);
   comparison.pct_delta = nan(n_rows, 1);

   % Only trust benchmark deltas when the embedded baseline was built from
   % the same benchmark-suite definition as the current checkout.
   if ~isstruct(baseline_meta) || ~isfield(baseline_meta, 'suite_signature') ...
         || isblanktext(baseline_meta.suite_signature) ...
         || string(baseline_meta.suite_signature) ~= string(current_signature)
      return
   end

   comparison.baseline_compatible(:) = true;

   % Guard against missing or incompatible baseline summary tables.
   if isempty(baseline_summary) || ...
         ~ismember('Name', baseline_summary.Properties.VariableNames)
      return
   end

   % Populate reference values where the baseline includes a matching
   % benchmark by name.
   for i = 1:n_rows
      hit = find(string(baseline_summary.Name) == string(comparison.Name(i)), ...
         1);
      if isempty(hit)
         continue
      end
      if ismember('Mean', baseline_summary.Properties.VariableNames)
         comparison.ref_mean(i) = baseline_summary.Mean(hit);
      end
      if ismember('StandardDeviation', baseline_summary.Properties.VariableNames)
         comparison.ref_std(i) = baseline_summary.StandardDeviation(hit);
      end
      if isfinite(comparison.ref_mean(i)) && comparison.ref_mean(i) > 0
         comparison.pct_delta(i) = ...
            100 * (comparison.Mean(i) - comparison.ref_mean(i)) ...
            / comparison.ref_mean(i);
      end
   end
end
