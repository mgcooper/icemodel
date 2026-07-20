function tests = test_suite_reports
   %TEST_SUITE_REPORTS Test Quarto source generation and suite summaries.

   tests = functiontests(localfunctions);
end

function test_regression_report_writes_plots_table_and_qmd(testCase)
   % A numerical report must expose physical and iteration changes visually.

   folder = temporaryFolder(testCase);
   results = regressionResults();
   report_file = icemodel.verification.report.buildTestSuiteReport( ...
      "regression", results, render=false, output_dir=folder);

   verifyEqual(testCase, report_file, ...
      fullfile(folder, "regression-suite-report.html"))
   verifyTrue(testCase, isfile(fullfile(folder, ...
      "regression-suite-summary.csv")))
   verifyTrue(testCase, isfile(fullfile(folder, ...
      "report-assets", "regression-percent-deltas.png")))
   verifyTrue(testCase, isfile(fullfile(folder, ...
      "report-assets", "regression-iteration-deltas.png")))
   verifyTrue(testCase, isfile(fullfile(folder, ...
      "report-assets", "regression-seb-closure.png")))
   source = fileread(fullfile(folder, "regression-suite-report.qmd"));
   verifySubstring(testCase, source, "**PASSED**")
   verifySubstring(testCase, source, "max iteration delta")
   verifySubstring(testCase, source, "closure seb max abs")
end

function test_performance_report_writes_current_ratio_and_history_plots(testCase)
   % A performance report must prioritize current and longitudinal timing plots.

   folder = temporaryFolder(testCase);
   baseline_root = fullfile(folder, "baselines");
   writePerfHistory(baseline_root)
   results = performanceResults();
   icemodel.verification.report.buildTestSuiteReport( ...
      "performance", results, render=false, output_dir=folder, ...
      baseline_root=baseline_root);

   verifyTrue(testCase, isfile(fullfile(folder, ...
      "report-assets", "performance-current-reference.png")))
   verifyTrue(testCase, isfile(fullfile(folder, ...
      "report-assets", "performance-runtime-ratio.png")))
   verifyTrue(testCase, isfile(fullfile(folder, ...
      "report-assets", "performance-history.png")))
   source = fileread(fullfile(folder, "performance-suite-report.qmd"));
   verifySubstring(testCase, source, "accepted rolling baselines")
   verifySubstring(testCase, source, "runtime ratio")
end

function test_report_requires_artifact_path_without_output_override(testCase)
   % Missing placement metadata must fail before writing ambiguous output.

   results = regressionResults();
   results = rmfield(results, "artifact_file");
   verifyError(testCase, @() ...
      icemodel.verification.report.buildTestSuiteReport( ...
      "regression", results, render=false), ...
      "icemodel:verification:report:missingArtifactPath")
end

function test_display_regression_summary_displays_accepted_baseline(testCase)
   % Baseline tables must not be silently skipped for lacking compare columns.

   fixture = regressionResults();
   baseline = fixture.report;
   keep = ~startsWith(string(baseline.Properties.VariableNames), "baseline_") ...
      & ~endsWith(string(baseline.Properties.VariableNames), "_pct_delta");
   baseline = baseline(:, keep); %#ok<NASGU>
   output = evalc( ...
      'icemodel.test.helpers.displayRegressionSummary(baseline)');
   verifySubstring(testCase, output, "icemodel_kanl_2016_solver2")
   verifySubstring(testCase, output, "max_Tice_numiter")
end

function test_perf_entrypoints_disable_default_profiling(testCase)
   % Formal timing entry points must keep profiling opt-in and disabled.

   test_root = icemodel.getpath("test");
   builder = fileread(fullfile(test_root, "tools", "build_perf_baseline.m"));
   runner = fileread(fullfile(test_root, "run_perf_suite.m"));
   start = strfind(builder, "kwargs.include_profile_artifacts");
   verifyNotEmpty(testCase, start)
   default_block = builder(start(1):(start(1) + 100));
   verifySubstring(testCase, default_block, "= false")
   verifyGreaterThanOrEqual(testCase, ...
      numel(strfind(builder, "profile off")), 2)
   verifySubstring(testCase, builder, ...
      "include_profile_artifacts=true requires one concrete smbmodel")
   verifySubstring(testCase, runner, "profile off")
end

function test_regression_entrypoint_resolves_one_batch_run_name(testCase)
   % Aggregate regression artifacts must share one directory before dispatch.

   runner = fileread(fullfile( ...
      icemodel.getpath("test"), "run_regression_suite.m"));
   resolve_at = strfind(runner, ...
      "icemodel.test.helpers.resolveRunStamp(run_name)");
   dispatch_at = strfind(runner, ...
      "per_model = arrayfun(@(mdl) runSingleModelRegression");
   verifyNumElements(testCase, resolve_at, 1)
   verifyNumElements(testCase, dispatch_at, 1)
   verifyLessThan(testCase, resolve_at, dispatch_at)
end

function folder = temporaryFolder(testCase)
   %TEMPORARYFOLDER Create one test-owned report directory.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   folder = string(fixture.Folder);
end

function results = regressionResults()
   %REGRESSIONRESULTS Return a two-case numerical comparison fixture.

   case_id = ["icemodel_kanl_2016_solver2"; "skinmodel_kanl_2016"];
   tier = repmat("full", 2, 1);
   smbmodel = ["icemodel"; "skinmodel"];
   sitename = repmat("kanl", 2, 1);
   solver = [2; NaN];
   runoff_pct_delta = [-1.2; -5.8];
   melt_pct_delta = [-0.8; -5.2];
   runoff_eval_pct_delta = [-1.0; -5.5];
   melt_eval_pct_delta = [-0.7; -5.0];
   mean_Tice_numiter = [2.5; 2.0];
   baseline_mean_Tice_numiter = [2.0; 1.5];
   max_Tice_numiter = [26.25; 14.5];
   baseline_max_Tice_numiter = [4.0; 4.0];
   n_not_converged = [0; 0];
   closure_seb_rmse = [3.1; 2.4];
   baseline_closure_seb_rmse = [3.0; 2.0];
   closure_seb_max_abs = [147.9; 42.0];
   baseline_closure_seb_max_abs = [140.0; 40.0];
   passed = [true; true];
   report = table(case_id, tier, smbmodel, sitename, solver, ...
      runoff_pct_delta, melt_pct_delta, runoff_eval_pct_delta, ...
      melt_eval_pct_delta, mean_Tice_numiter, ...
      baseline_mean_Tice_numiter, max_Tice_numiter, ...
      baseline_max_Tice_numiter, n_not_converged, closure_seb_rmse, ...
      baseline_closure_seb_rmse, closure_seb_max_abs, ...
      baseline_closure_seb_max_abs, passed);
   results = struct('report', report, 'passed', true, ...
      'artifact_file', "fixture.mat");
end

function results = performanceResults()
   %PERFORMANCERESULTS Return a two-case performance comparison fixture.

   case_id = ["icemodel_kanl_2016_solver2"; "skinmodel_kanl_2016"];
   tier = repmat("full", 2, 1);
   smbmodel = ["icemodel"; "skinmodel"];
   sitename = repmat("kanl", 2, 1);
   solver = [2; NaN];
   median_wall_s = [12; 8];
   ref_wall_s = [10; 10];
   floor_wall_s = [8; 8];
   gate_wall_s = [12; 12];
   baseline_compatible = [true; true];
   passed_perf = [true; true];
   compare_reason = [""; ""];
   case_summary = table(case_id, tier, smbmodel, sitename, solver, ...
      median_wall_s, ref_wall_s, floor_wall_s, gate_wall_s, ...
      baseline_compatible, passed_perf, compare_reason);
   results = struct('case_summary', case_summary, 'passed', true, ...
      'artifact_file', "fixture.mat");
end

function writePerfHistory(baseline_root)
   %WRITEPERFHISTORY Save two accepted baseline snapshots for plotting.

   case_id = ["icemodel_kanl_2016_solver2"; "skinmodel_kanl_2016"];
   for k = 1:2
      median_wall_s = [9 + k; 7 + k];
      last_updated_utc = repmat(datetime(2026, k, 1, TimeZone="UTC"), 2, 1);
      PerfBaseline = table(case_id, median_wall_s, last_updated_utc);
      folder = fullfile(baseline_root, "archive", "snapshot" + string(k));
      mkdir(folder)
      save(fullfile(folder, ...
         "perf_baseline_2016_rolling_fixture.mat"), "PerfBaseline")
   end
end
