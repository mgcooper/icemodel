import matlab.perftest.TimeExperiment

%% Use the default test suite with the script-based test

% This test indicates that renameing w/ intersect is faster by almost 2x and has
% lower standard deviation and tighter min/max range, but rounding is generally
% faster with the ismember method. NOTE: I am not sure if the script-based test
% framework is valid in this case b/c the data might get modified i.e. I am not
% sure of the scope.

results = runperf("test_renameRoundScript");
sampleSummary(results)

%% Use the default test suite with the function-based test

% This test indicates that renameing w/ ismember is faster but they are nearly
% identical. I prefer the ismember method so I will use that.
% 
% Rounding is generally fastest with intersect then switch then ismember but
% they are so close its immaterial. Maybe intersect has enough of an edge to
% matter over many iterations.

results = runperf("test_renameRoundFunc");
sampleSummary(results)

%% Use the example from the documentation
suite = testsuite("test_renameRoundScript");
experiment = TimeExperiment.limitingSamplingError("NumWarmups",2, ...
    "RelativeMarginOfError",0.04,"ConfidenceLevel",0.98);
resultsTE = run(experiment,suite);
sampleSummary(resultsTE)

%% run the example test

% This one is quite slow
% results = runperf('PreallocationTest');
% sampleSummary(results)

%% fprintfTest example

results = runperf("fprintfTest")
results(1)
results(1).Samples

% Display the mean measured time for the first test. To exclude data collected
% in the warm-up runs, use the values in the Samples property.
sampleTimes = results(1).Samples.MeasuredTime;
meanTest = mean(sampleTimes)

% Compute Statistics for All Test Elements
% To compare the different calls to fprintf, create a table of summary
% statistics from results. In this example, both test methods write the same
% amount of data to a file. Therefore, some of the difference between the
% statistical values is attributed to calling the fprintf function with an
% output argument.
sampleSummary(results)

% Change the statistical objectives defined by the runperf function by
% constructing and running a time experiment. Construct a time experiment with
% measurements that reach a sample mean with a 2% relative margin of error
% within a 98% confidence level. Collect 4 warm-up measurements and up to 16
% sample measurements.     

suite = testsuite("fprintfTest");

% Construct a time experiment with the specified requirements, and run the
% tests. In this example, the performance testing framework is not able to meet
% the stricter statistical objectives with the specified number of maximum
% samples. Your results might vary.
experiment = TimeExperiment.limitingSamplingError("NumWarmups",4, ...
    "MaxSamples",16,"RelativeMarginOfError",0.02,"ConfidenceLevel",0.98);
resultsTE = run(experiment,suite);

% Increase the maximum number of samples to 32 and rerun the time experiment.
experiment = TimeExperiment.limitingSamplingError("NumWarmups",4, ...
    "MaxSamples",32,"RelativeMarginOfError",0.02,"ConfidenceLevel",0.98);
resultsTE = run(experiment,suite);

% Compute the summary statistics for the test elements.
T1 = sampleSummary(resultsTE)

%% Measure First-Time Cost

% Start a new MATLAB® session. A new session ensures that MATLAB has not run the
% code contained in your tests. 

% Measure the first-time cost of your code by creating and running a fixed time
% experiment with zero warm-up measurements and one sample measurement. 

% Create a test suite. Because you are measuring the first-time cost of a
% function, run a single test. To run multiple tests, save the results and start
% a new MATLAB session between tests.  

suite = testsuite("fprintfTest/testPrintingToFile");

% Construct and run the time experiment.
experiment = TimeExperiment.withFixedSampleSize(1);
results = run(experiment,suite);

% Display the results. The TestActivity table shows that there were no warm-up
% measurements. 
fullTable = results.TestActivity




