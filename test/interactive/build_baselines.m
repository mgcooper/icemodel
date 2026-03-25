%BUILD_BASELINES Rebuild rolling baselines and optionally snapshot a release.
%
% Usage:
%   build_baselines  % Rebuild rolling baselines only.
%   build_baselines  % Then set do_snapshot = true below to tag a release.
%
% After building, smoke perf and regression suite runs verify the new
% baselines. Review the summary output before proceeding.

% Configuration
do_snapshot = false;
baseline_tag = ['v' icemodel.internal.version]; % used when do_snapshot is true

% Build rolling baselines. The standard rebuild is to run with default options.
RegressionBaseline = build_regression_baseline();
PerfBaseline = build_perf_baseline();

% Summaries are printed when build_* functions are run.
% For reference, this displays the same baseline results.
icemodel.test.helpers.displayPerfSummary(PerfBaseline)

% The file name can also be set and passed to the display function
filename = icemodel.test.helpers.baselineFilePath("perf");
icemodel.test.helpers.displayPerfSummary(filename)

% Verify with a smoke perf run
perf_results = run_perf_suite( ...
   tier="smoke", smbmodel="icemodel", solver=2, include_benchmarks=false);

% Display perf results
icemodel.test.helpers.displayPerfResults(perf_results)

% Verify with a smoke regression run
regression_results = run_regression_suite( ...
   tier="smoke", smbmodel="icemodel", solver=2);

% Display regression results
icemodel.test.helpers.displayRegressionResults(regression_results)

%% Scratch space to demo path building and results display

filename = icemodel.test.helpers.baselineFilePath("perf");
icemodel.test.helpers.displayPerfSummary(filename)

filename = icemodel.test.helpers.baselineFilePath("perf", ...
   baseline_type="release", baseline_tag="v1.1");
icemodel.test.helpers.displayPerfSummary(filename)

% This works but the file doesn't exist
icemodel.test.helpers.baselineFilePath("perf", ...
   smbmodel="skinmodel", ...
   baseline_type="release", ...
   baseline_tag="v1.0.1")
icemodel.test.helpers.displayPerfSummary(filename)

%% Snapshots

% Snapshot to release baseline (optional)
if do_snapshot == true

   models = unique(PerfBaseline.smbmodel);

   for i = 1:numel(models)

      % Snapshot perf
      icemodel.test.helpers.snapshotBaseline( ...
         "perf", baseline_tag, models(i), true);

      % Snapshot regression
      icemodel.test.helpers.snapshotBaseline( ...
         "regression", baseline_tag, models(i), true);
   end
   fprintf('Snapshot saved as %s\n', baseline_tag)
end
