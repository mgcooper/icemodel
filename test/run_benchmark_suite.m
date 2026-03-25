function results = run_benchmark_suite(options)
   %RUN_BENCHMARK_SUITE Run component benchmarks under test/benchmarks.
   %
   %  results = run_benchmark_suite()
   %  results = run_benchmark_suite(testname="SebKernelPerfTest")
   %  results = run_benchmark_suite(testname="micro", include_subfolders=true)
   %  results = run_benchmark_suite(sampling_profile="fast")
   %
   % Use this runner for component benchmarks and one-off performance
   % experiments. The default suite uses only the top-level benchmark files
   % under test/benchmarks/. Opt-in microbenchmarks can live in subfolders.
   %
   % For formal performance regression against managed rolling or release
   % baselines, use `run_perf_suite(...)` instead.
   %
   % See RenameRoundTest for notes on results that informed postprocess.

   arguments (Input)

      options.testname (1, :) string ...
         = ""

      options.include_subfolders (1, 1) logical ...
         = false

      options.sampling_profile (1, :) string ...
         {icemodel.validators.mustBeBenchmarkSamplingProfileName( ...
         options.sampling_profile)} ...
         = "default"

      options.num_warmups (1, 1) double ...
         = NaN

      options.max_samples (1, 1) double ...
         = NaN

      options.relative_margin_of_error (1, 1) double ...
         = NaN

      options.confidence_level (1, 1) double ...
         = NaN

      options.show_summary (1, 1) logical ...
         = true
   end

   import matlab.perftest.TimeExperiment

   % Install the canonical formal-suite config and keep the cleanup handle
   % in scope so the caller's environment is restored on return.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Set the experiment sampling options.
   [num_warmups, max_samples, relative_margin_of_error, ...
      confidence_level] = resolveSamplingOptions(options);

   % Resolve the requested benchmark selector before constructing the
   % limiting-sampling-error experiment.
   suite = buildBenchmarkSuite(options);

   % Run the component benchmark experiment.
   experiment = TimeExperiment.limitingSamplingError( ...
      "NumWarmups", num_warmups, ...
      "MaxSamples", max_samples, ...
      "RelativeMarginOfError", relative_margin_of_error, ...
      "ConfidenceLevel", confidence_level);
   results = run(experiment, suite);

   if options.show_summary
      disp(compactBenchmarkSummary(results));
   end
end

function suite = buildBenchmarkSuite(options)
   %BUILDBENCHMARKSUITE Resolve the requested benchmark file/folder selector.

   % Resolve simple benchmark names to the formal benchmark file/folder.
   target = resolveBenchmarkTarget(options.testname);

   if exist(target, 'dir') == 7
      if options.include_subfolders
         suite = testsuite(target, 'IncludingSubfolders', true);
      else
         suite = testsuite(target);
      end
   else
      suite = testsuite(target);
   end
end

function [num_warmups, max_samples, relative_margin_of_error, ...
      confidence_level] = resolveSamplingOptions(options)
   %RESOLVESAMPLINGOPTIONS Merge profile defaults with explicit overrides.

   [num_warmups, max_samples, relative_margin_of_error, ...
      confidence_level] = profileDefaults(options.sampling_profile);

   % Explicit numeric overrides win over the named profile defaults.
   if ~isnan(options.num_warmups)
      num_warmups = options.num_warmups;
   end
   if ~isnan(options.max_samples)
      max_samples = options.max_samples;
   end
   if ~isnan(options.relative_margin_of_error)
      relative_margin_of_error = options.relative_margin_of_error;
   end
   if ~isnan(options.confidence_level)
      confidence_level = options.confidence_level;
   end

   mustBeInteger(num_warmups);
   mustBeNonnegative(num_warmups);
   mustBeInteger(max_samples);
   mustBePositive(max_samples);
   mustBePositive(relative_margin_of_error);
   mustBePositive(confidence_level);
   mustBeLessThanOrEqual(confidence_level, 1.0);
end

function target = resolveBenchmarkTarget(testname)
   %RESOLVEBENCHMARKTARGET Map a selector string to a benchmark file/folder.

   testdir = icemodel.getpath('test');
   benchdir = fullfile(testdir, 'benchmarks');

   % Default to the top-level core benchmark suite. Microbenchmarks live in
   % subfolders and are opt-in through explicit selectors.
   if isblanktext(testname)
      target = benchdir;
      return
   end

   if contains(testname, filesep) || endsWith(testname, ".m")
      target = char(testname);
      return
   end

   target = fullfile(benchdir, char(testname) + ".m");
   if exist(target, 'file') == 2
      return
   end

   target = fullfile(benchdir, char(testname));
   if exist(target, 'dir') == 7
      return
   end

   dirhits = dir(fullfile(benchdir, '**', char(testname)));
   dirhits = dirhits([dirhits.isdir]);
   if isscalar(dirhits)
      target = fullfile(dirhits(1).folder, dirhits(1).name);
      return
   elseif numel(dirhits) > 1
      error('benchmark selector is ambiguous: %s', testname)
   end

   filehits = dir(fullfile(benchdir, '**', ...
      char(testname) + ".m"));
   if isscalar(filehits)
      target = fullfile(filehits(1).folder, filehits(1).name);
      return
   elseif numel(filehits) > 1
      error('benchmark selector is ambiguous: %s', testname)
   end

   error('benchmark selector does not resolve to a file/folder: %s', testname)
end

function [num_warmups, max_samples, relative_margin_of_error, ...
      confidence_level] = profileDefaults(profile_name)
   %PROFILEDEFAULTS Return the named sampling policy for benchmark runs.

   switch char(profile_name)
      case 'fast'
         num_warmups = 2;
         max_samples = 20;
         relative_margin_of_error = 0.10;
         confidence_level = 0.90;

      case 'default'
         num_warmups = 4;
         max_samples = 300;
         relative_margin_of_error = 0.02;
         confidence_level = 0.98;

      case 'strict'
         num_warmups = 6;
         max_samples = 600;
         relative_margin_of_error = 0.01;
         confidence_level = 0.99;

      otherwise
         error('unsupported sampling profile: %s', profile_name)
   end
end

function summary = compactBenchmarkSummary(results)
   %COMPACTBENCHMARKSUMMARY Keep only the main timing columns for display.

   summary = sampleSummary(results);
   keep = ["Name", "SampleSize", "Mean", "StandardDeviation"];
   keep = keep(ismember(keep, string(summary.Properties.VariableNames)));
   summary = summary(:, cellstr(keep));
end
