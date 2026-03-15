function results = runperf(varargin)

   import matlab.perftest.TimeExperiment

   if nargin < 1
      testname = "SebSolverTest";
   else
      testname = string(varargin{1});
   end

   % Resolve the simple benchmark names used in this helper to the component
   % benchmark folder. Callers can still pass an explicit path or selector.
   suite = testsuite(resolveBenchmarkTarget(testname));
   experiment = TimeExperiment.limitingSamplingError("NumWarmups", 4, ...
      "MaxSamples", 100, "RelativeMarginOfError", 0.02, "ConfidenceLevel", 0.98);
   results = run(experiment, suite);
   sampleSummary(results)
   
   % % This is the simplest way to run
   % results = runperf(testname);
   %
   % % Display the mean measured times for each test. To exclude data
   % % collected in the warm-up runs, use the values in the Samples property.
   % for n = 1:numel(results)
   %    sampleTimes = results(n).Samples.MeasuredTime;
   %    printf(mean(sampleTimes), 6)
   % end
   %
   % results.sampleSummary
   
   
   %% Use the default test suite with the script-based test

   % This test indicates that renameing w/ intersect is faster by almost 2x and
   % has lower standard deviation and tighter min/max range, but rounding is
   % generally faster with the ismember method. NOTE: I am not sure if the
   % script-based test framework is valid in this case b/c the data might get
   % modified i.e. I am not sure of the scope.

   
   % testname = "test_renameRoundScript";
   % results = runperf(testname);
   % sampleSummary(results)

   %% Use the default test suite with the function-based test

   % This test indicates that renameing w/ ismember is faster but they are nearly
   % identical. I prefer the ismember method so I will use that.
   %
   % Rounding is generally fastest with intersect then switch then ismember but
   % they are so close its immaterial. Maybe intersect has enough of an edge to
   % matter over many iterations.

   % results = runperf("test_renameRoundFunc");
   % sampleSummary(results)

   %% Use the example from the documentation
   % suite = testsuite(resolveBenchmarkTarget("test_renameRoundScript"));
   % experiment = TimeExperiment.limitingSamplingError("NumWarmups",2, ...
   %    "RelativeMarginOfError",0.04,"ConfidenceLevel",0.98);
   % resultsTE = run(experiment,suite);
   % sampleSummary(resultsTE)

   %% run the example test

end

function target = resolveBenchmarkTarget(testname)
   testname = string(testname);
   if contains(testname, filesep) || endsWith(testname, ".m")
      target = char(testname);
      return
   end
   rootdir = icemodel.internal.fullpath();
   target = fullfile(rootdir, 'test', 'benchmarks', char(testname) + ".m");
end
