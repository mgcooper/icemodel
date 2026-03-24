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
baseline_tag = "v1.0";  % used only when do_snapshot is true

% Build rolling baselines. The standard rebuild is to run with default options.
RegressionBaseline = build_regression_baseline();
PerfBaseline = build_perf_baseline();

% Display baseline results
icemodel.test.helpers.displayPerfSummary(PerfBaseline)

% Verify with a smoke perf run
perf_results = run_perf_suite( ...
   tier="smoke", smbmodel="icemodel", solver=2, include_benchmarks=false);

% Verify with a smoke regression run
regression_results = run_regression_suite( ...
   tier="smoke", smbmodel="icemodel", solver=2);

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
