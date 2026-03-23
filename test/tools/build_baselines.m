
% The standard rebuild is to run with default options.
RegressionBaseline = build_regression_baseline();
PerfBaseline = build_perf_baseline();

% Display the results
PerfBaseline(:, ["case_id", "median_wall_s", "n_runs", "n_warmups"])

disp([PerfBaseline.case_id, ...
   PerfBaseline.median_wall_s ./ (PerfBaseline.n_runs + PerfBaseline.n_warmups)])

% After building, run the perf suite
perf_results = run_perf_suite( ...
   tier="smoke", smbmodel="icemodel", solver=2, include_benchmarks=false);
