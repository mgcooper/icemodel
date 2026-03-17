function results = run_benchmark_suite(options)
   %RUN_BENCHMARK_SUITE Run component benchmarks under test/benchmarks.
   %
   %  results = run_benchmark_suite()
   %  results = run_benchmark_suite(testname="SebKernelPerfTest")
   %  results = run_benchmark_suite(testname="micro", include_subfolders=true)
   %
   % Use this runner for component benchmarks and one-off performance
   % experiments. The default suite uses only the top-level benchmark files
   % under test/benchmarks/. Opt-in microbenchmarks can live in subfolders.
   %
   % For formal performance regression against managed rolling or release
   % baselines, use `run_perf_suite(...)` instead.
   %
   % See RenameRoundTest for notes on results that informed postprocess.

   arguments
      options.testname (1, :) string = ""
      options.include_subfolders (1, 1) logical = false
      options.num_warmups (1, 1) double {mustBeInteger, mustBeNonnegative} = 4
      options.max_samples (1, 1) double {mustBeInteger, mustBePositive} = 100
      options.relative_margin_of_error (1, 1) double {mustBePositive} = 0.02
      options.confidence_level (1, 1) double {mustBePositive, ...
         mustBeLessThanOrEqual(options.confidence_level, 1.0)} = 0.98
      options.show_summary (1, 1) logical = true
   end

   import matlab.perftest.TimeExperiment

   thisdir = fileparts(mfilename('fullpath'));
   addpath(fullfile(fileparts(thisdir), 'icemodel'))
   addpath(fullfile(fileparts(thisdir), 'icemodel', 'dependencies'))
   addpath(thisdir)

   % Resolve the simple benchmark names used in this helper to the formal
   % benchmark suite. Callers can still pass an explicit path or selector.
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
   experiment = TimeExperiment.limitingSamplingError( ...
      "NumWarmups", options.num_warmups, ...
      "MaxSamples", options.max_samples, ...
      "RelativeMarginOfError", options.relative_margin_of_error, ...
      "ConfidenceLevel", options.confidence_level);
   results = run(experiment, suite);

   if options.show_summary
      disp(sampleSummary(results));
   end
end

function target = resolveBenchmarkTarget(testname)
   rootdir = icemodel.internal.fullpath();

   if isempty(testname) || all(strlength(testname) == 0)
      target = fullfile(rootdir, 'test', 'benchmarks');
      return
   end

   testname = string(testname);
   if contains(testname, filesep) || endsWith(testname, ".m")
      target = char(testname);
      return
   end

   target = fullfile(rootdir, 'test', 'benchmarks', char(testname) + ".m");
   if exist(target, 'file') == 2
      return
   end

   target = fullfile(rootdir, 'test', 'benchmarks', char(testname));
   if exist(target, 'dir') == 7
      return
   end

   dirhits = dir(fullfile(rootdir, 'test', 'benchmarks', '**', char(testname)));
   dirhits = dirhits([dirhits.isdir]);
   if numel(dirhits) == 1
      target = fullfile(dirhits(1).folder, dirhits(1).name);
      return
   elseif numel(dirhits) > 1
      error('benchmark selector is ambiguous: %s', testname)
   end

   filehits = dir(fullfile(rootdir, 'test', 'benchmarks', '**', ...
      char(testname) + ".m"));
   if numel(filehits) == 1
      target = fullfile(filehits(1).folder, filehits(1).name);
      return
   elseif numel(filehits) > 1
      error('benchmark selector is ambiguous: %s', testname)
   end

   error('benchmark selector does not resolve to a file/folder: %s', testname)
end
