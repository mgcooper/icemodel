function tests = test_baseline_contracts
   %TEST_BASELINE_CONTRACTS Verify baseline-selection and legacy-load helpers.
   tests = functiontests(localfunctions);
end

function test_resolveBaselineSelector_handles_rolling_and_release(testCase)
   % Cover the public selector contract for rolling and release baselines.

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector("rolling");
   testCase.verifyEqual(baseline_type, "rolling");
   testCase.verifyEqual(baseline_tag, "");

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector("v1.1");
   testCase.verifyEqual(baseline_type, "release");
   testCase.verifyEqual(baseline_tag, "v1.1");
end

function test_resolveBaselineBuild_defaults_blank_selector_to_rolling(testCase)
   % Baseline-build helpers should treat a fully blank selector as the
   % managed rolling baseline rather than forcing the caller to restate it.

   [baseline_type, baseline_tag, output_file] = ...
      icemodel.test.helpers.resolveBaselineBuild( ...
      "regression", "", "", "skinmodel", "", NaN);

   testCase.verifyEqual(baseline_type, "rolling");
   testCase.verifyEqual(baseline_tag, "");
   testCase.verifyTrue(contains(string(output_file), ...
      "regression_baseline_rolling_skinmodel.mat"));
end

function test_loadBaseline_regression_normalizes_legacy_schema(testCase)
   % Legacy baseline rows should be upgraded into the current schema on
   % load so older accepted files remain usable during the transition.

   filepath = [tempname '.mat'];
   cleanup = onCleanup(@() deleteIfExists(filepath));
   legacy_time = datenum(datetime(2026, 3, 9, 4, 15, 39, ...
      'TimeZone', 'UTC')); %#ok<DATNM>

   baseline = table( ...
      "smoke_icemodel_kanm_2016_bc2", ...
      1.5, ...
      legacy_time, ...
      'VariableNames', {'case_id', 'runoff_final', 'last_updated_utc'});
   save(filepath, 'baseline');

   loaded = icemodel.test.helpers.loadBaseline("regression", ...
      smbmodel="icemodel", filename=filepath);

   testCase.verifyEqual(loaded.case_id, "icemodel_kanm_2016_solver2");
   testCase.verifyEqual(loaded.baseline_tag, "");
   testCase.verifyEqual(loaded.baseline_type, "rolling");
   testCase.verifyEqual(loaded.smbmodel_filter, "icemodel");
   testCase.verifyTrue(isdatetime(loaded.last_updated_utc));
   clear cleanup
end

function test_loadBaseline_regression_prefers_explicit_path(testCase)
   % Explicit baseline paths should still load rolling files even when the
   % rolling selector has already been normalized to a blank tag.

   filepath = [tempname '.mat'];
   cleanup = onCleanup(@() deleteIfExists(filepath));

   RegressionBaseline = table( ...
      "icemodel_kanm_2016_solver2", ...
      1.5, ...
      'VariableNames', {'case_id', 'runoff_final'});
   save(filepath, 'RegressionBaseline');

   loaded = icemodel.test.helpers.loadBaseline("regression", ...
      smbmodel="icemodel", filename=filepath);

   testCase.verifyEqual(loaded.case_id, "icemodel_kanm_2016_solver2");
   testCase.verifyEqual(loaded.runoff_final, 1.5);
   clear cleanup
end

function test_getRegressionCaseMatrix_keeps_skinmodel_under_filter(testCase)
   % Solver filters should narrow only icemodel rows when the formal suite
   % is asked for smbmodel="all".

   cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
      tier="full", smbmodel="all", solver=2);

   testCase.verifyTrue(any(cases.smbmodel == "skinmodel"));
   testCase.verifyEqual(unique(cases.solver(cases.smbmodel == "skinmodel")), 1);
   testCase.verifyEqual(unique(cases.solver(cases.smbmodel == "icemodel")), 2);
end

function test_getPerfCaseMatrix_keeps_skinmodel_under_filter(testCase)
   % The same solver-filter contract should hold for the managed perf
   % matrix so build/run/bootstrap paths stay consistent.

   cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier="full", smbmodel="all", solver=2);

   testCase.verifyTrue(any(cases.smbmodel == "skinmodel"));
   testCase.verifyEqual(unique(cases.solver(cases.smbmodel == "skinmodel")), 1);
   testCase.verifyEqual(unique(cases.solver(cases.smbmodel == "icemodel")), 2);
end

function test_testSmbmodel_group_stays_limited_to_formal_models(testCase)
   % New smbmodel options should not enter the formal suite automatically
   % before they have accepted baselines and runner coverage.

   models = icemodel.namelists.smbmodel("test");

   testCase.verifyEqual(models(:), ["icemodel"; "skinmodel"]);
end

function test_testsmbmodel_completion_group_matches_formal_models(testCase)
   % Runner completions should expose only the formal smbmodel group plus
   % the aggregate selector used by suite entrypoints.

   models = icemodel.namelists.testsmbmodel();

   testCase.verifyEqual(models(:), ["all"; "icemodel"; "skinmodel"]);
end

function test_resolveRequestedSmbmodels_expands_virtual_all_selector(testCase)
   % Formal-suite entrypoints should expand the virtual aggregate selector
   % once, then iterate through the canonical single-model workflow.

   models = icemodel.test.helpers.resolveRequestedSmbmodels("all");

   testCase.verifyEqual(models(:), ["icemodel"; "skinmodel"]);
end

function test_resolveRequestedSmbmodels_preserves_single_model_selector(testCase)
   % Concrete formal model selectors should flow through unchanged.

   models = icemodel.test.helpers.resolveRequestedSmbmodels("skinmodel");

   testCase.verifyEqual(models(:), "skinmodel");
end

function test_setModelOptsForCase_defaults_to_two_year_contract(testCase)
   % Formal single-year case definitions should expand to one spinup year
   % plus the retained comparison year.

   c = struct('smbmodel', "icemodel", 'sitename', "kanm", ...
      'forcings', "kanm", 'userdata', "", 'uservars', "", ...
      'simyear', 2016, 'solver', 2);

   opts = icemodel.test.helpers.setModelOptsForCase(c);

   testCase.verifyEqual(opts.simyears, [2015 2016]);
   testCase.verifyEqual(opts.n_spinup_years, 1);
   testCase.verifyEqual(opts.output_years, 2016);
end

function test_bootstrapTestEnvironment_restores_caller_config(testCase)
   % The suite bootstrap should install the canonical demo config for the
   % run, then restore the caller's previous config on cleanup.

   previous_output = getenv('ICEMODEL_OUTPUT_PATH');
   restore_output = onCleanup(@() setenv('ICEMODEL_OUTPUT_PATH', previous_output));

   custom_output = fullfile(tempdir, 'codex_custom_output');
   setenv('ICEMODEL_OUTPUT_PATH', custom_output);
   test_cfg = icemodel.config('casename', 'test', 'setenv', false);

   [~, ~, output_path, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.verifyClass(suite_cleanup, 'onCleanup');
   testCase.verifyEqual(string(output_path), string(test_cfg.ICEMODEL_OUTPUT_PATH));

   clear suite_cleanup
   testCase.verifyEqual(string(getenv('ICEMODEL_OUTPUT_PATH')), ...
      string(custom_output));

   clear restore_output
end

function test_loadBaseline_perf_returns_saved_metadata(testCase)
   % Perf baseline loads should expose the saved build metadata so runners can
   % decide whether whole-model wall-time comparison is fair.

   filepath = [tempname '.mat'];
   cleanup = onCleanup(@() deleteIfExists(filepath));

   PerfBaseline = table( ...
      "icemodel_kanm_2016_solver2", ...
      12.3, ...
      'VariableNames', {'case_id', 'median_wall_s'});
   meta = struct('matlab_version', "24.2.0 (R2024b)", 'host', "MACA64");
   save(filepath, 'PerfBaseline', 'meta');

   [loaded, loaded_meta] = icemodel.test.helpers.loadBaseline("perf", ...
      smbmodel="icemodel", filename=filepath);

   testCase.verifyEqual(loaded.case_id, "icemodel_kanm_2016_solver2");
   testCase.verifyEqual(loaded_meta.matlab_version, "24.2.0 (R2024b)");
   testCase.verifyEqual(loaded_meta.host, "MACA64");

   clear cleanup
end

function test_baselineFilePath_returns_rolling_perf_by_default(testCase)
   % The default call should return the rolling perf baseline for icemodel.

   pathname = icemodel.test.helpers.baselineFilePath("perf");
   testCase.verifyTrue(contains(pathname, ...
      "perf_baseline_2016_rolling_icemodel.mat"));
end

function test_baselineFilePath_accepts_smbmodel_name_value(testCase)
   % smbmodel should be a name-value argument.

   pathname = icemodel.test.helpers.baselineFilePath("perf", ...
      smbmodel="skinmodel");
   testCase.verifyTrue(contains(pathname, ...
      "perf_baseline_2016_rolling_skinmodel.mat"));
end

function test_baselineFilePath_tag_implies_release(testCase)
   % Providing a baseline_tag without baseline_type should infer release.

   pathname = icemodel.test.helpers.baselineFilePath("perf", ...
      baseline_tag="v1.1");
   testCase.verifyTrue(contains(pathname, ...
      "perf_baseline_2016_v1_1_icemodel.mat"));
end

function test_baselineFilePath_returns_rolling_regression(testCase)
   % Regression baselines do not include simyear in the filename.

   pathname = icemodel.test.helpers.baselineFilePath("regression");
   testCase.verifyTrue(contains(pathname, ...
      "regression_baseline_rolling_icemodel.mat"));
end

function test_baselineFilePath_resolves_latest_release(testCase)
   % baseline_type="release" without a tag should resolve the latest version.

   pathname = icemodel.test.helpers.baselineFilePath("perf", ...
      baseline_type="release");
   testCase.verifyTrue(contains(pathname, "v1_1"));
end

function test_loadBaseline_perf_returns_nonempty_table(testCase)
   % Loading the rolling perf baseline should return a populated table.

   baseline = icemodel.test.helpers.loadBaseline("perf");
   testCase.verifyFalse(isempty(baseline));
   testCase.verifyTrue(istable(baseline));
end

function test_loadBaseline_regression_returns_nonempty_table(testCase)
   % Loading the rolling regression baseline should return a populated table.

   baseline = icemodel.test.helpers.loadBaseline("regression");
   testCase.verifyFalse(isempty(baseline));
   testCase.verifyTrue(istable(baseline));
end

function test_artifactFilePath_returns_existing_file(testCase)
   % The default call should find the most recent perf artifact.

   pathname = icemodel.test.helpers.artifactFilePath("perf");
   testCase.verifyTrue(exist(pathname, 'file') == 2);
end

function test_referenceFilePath_returns_runoff_path(testCase)
   % The runoff reference path should point to the references directory.

   pathname = icemodel.test.helpers.referenceFilePath("runoff");
   testCase.verifyTrue(contains(pathname, "references"));
   testCase.verifyTrue(contains(pathname, "runoff_reference"));
end

function test_loadReference_returns_nonempty_table(testCase)
   % Loading the runoff reference should return a populated table.

   ref = icemodel.test.helpers.loadReference("runoff");
   testCase.verifyFalse(isempty(ref));
   testCase.verifyTrue(istable(ref));
end

function deleteIfExists(filepath)
   %DELETEIFEXISTS Remove one file if it exists.

   if exist(filepath, 'file') == 2
      delete(filepath);
   end
end
