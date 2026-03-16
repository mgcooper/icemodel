function results = run_benchmark_suite(varargin)
   %RUN_BENCHMARK_SUITE Run component benchmarks under test/benchmarks.
   %
   %  results = run_benchmark_suite()
   %  results = run_benchmark_suite("SebSolverTest")
   %
   % Use this runner for component benchmarks and one-off performance
   % experiments. For formal performance regression against managed rolling or
   % release baselines, use `run_perf_suite(...)` instead.
   %
   % See RoundRenameTest for notes on results that informed postprocess.

   import matlab.perftest.TimeExperiment

   if nargin < 1
      testname = "";
   else
      testname = string(varargin{1});
   end

   thisdir = fileparts(mfilename('fullpath'));
   addpath(fullfile(fileparts(thisdir), 'icemodel'))
   addpath(thisdir)

   % Resolve the simple benchmark names used in this helper to the formal
   % benchmark suite. Callers can still pass an explicit path or selector.
   target = resolveBenchmarkTarget(testname);
   if exist(target, 'dir') == 7
      suite = testsuite(target, 'IncludingSubfolders', true);
   else
      suite = testsuite(target);
   end
   experiment = TimeExperiment.limitingSamplingError("NumWarmups", 4, ...
      "MaxSamples", 100, "RelativeMarginOfError", 0.02, "ConfidenceLevel", 0.98);
   results = run(experiment, suite);

   % Display the results
   sampleSummary(results)
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
end
