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

   benchmark = struct();
   benchmark.summary = table();
   benchmark.comparison = table();
   benchmark.meta = struct();

   if ~include_benchmarks
      return
   end

   bench_results = run_benchmark_suite( ...
      sampling_profile=benchmark_sampling_profile, show_summary=false);
   [current_signature, suite_files] = ...
      icemodel.test.helpers.benchmarkSuiteSignature();
   benchmark.summary = sampleSummary(bench_results);
   if ~isempty(benchmark.summary) ...
         && ismember('Name', benchmark.summary.Properties.VariableNames)
      benchmark.summary.Name = string(benchmark.summary.Name);
   end
   benchmark.summary.Valid = reshape(logical([bench_results.Valid]), [], 1);

   [baseline_summary, baseline_meta, source_file] = ...
      loadBenchmarkBaselineFromPerf(simyear, baseline_tag, smbmodel);
   benchmark.comparison = compareBenchmarkSummary( ...
      benchmark.summary, baseline_summary, baseline_meta, current_signature);
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

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_tag);

   if smbmodel == "all"
      models = icemodel.namelists.smbmodel("test");
   else
      models = string(smbmodel);
   end

   BenchmarkBaseline = table();
   meta = struct();
   source_file = "";

   for i = 1:numel(models)
      candidate = icemodel.test.helpers.defaultBaselinePath( ...
         "perf", baseline_type, baseline_tag, models(i), simyear);
      if exist(char(candidate), 'file') ~= 2
         continue
      end

      matvars = string({whos('-file', char(candidate)).name});
      has_baseline = any(matvars == "BenchmarkBaseline");
      has_meta = any(matvars == "benchmark_meta");
      if ~has_baseline && ~has_meta
         continue
      end

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

   if isempty(current_summary)
      comparison = table();
      return
   end

   comparison = current_summary(:, intersect(["Name", "SampleSize", "Mean", ...
      "StandardDeviation"], string(current_summary.Properties.VariableNames), ...
      'stable'));
   if ~ismember('Name', comparison.Properties.VariableNames)
      return
   end

   n_rows = height(comparison);
   comparison.baseline_compatible = false(n_rows, 1);
   comparison.ref_mean = nan(n_rows, 1);
   comparison.ref_std = nan(n_rows, 1);
   comparison.pct_delta = nan(n_rows, 1);

   if ~isstruct(baseline_meta) || ~isfield(baseline_meta, 'suite_signature') ...
         || isblanktext(baseline_meta.suite_signature) ...
         || string(baseline_meta.suite_signature) ~= string(current_signature)
      return
   end

   comparison.baseline_compatible(:) = true;

   if isempty(baseline_summary) || ...
         ~ismember('Name', baseline_summary.Properties.VariableNames)
      return
   end

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
