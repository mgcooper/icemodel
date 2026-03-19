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

function test_loadRegressionBaseline_returns_empty_for_blank_selector(testCase)
   % A blank selector is the explicit "no baseline" case.

   baseline = icemodel.test.helpers.loadRegressionBaseline("", "all");
   testCase.verifyTrue(istable(baseline));
   testCase.verifyEmpty(baseline);
end

function test_loadRegressionBaseline_normalizes_legacy_schema(testCase)
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

   loaded = icemodel.test.helpers.loadRegressionBaseline( ...
      "rolling", "icemodel", filepath);

   testCase.verifyEqual(loaded.case_id, "icemodel_kanm_2016_solver2");
   testCase.verifyEqual(loaded.baseline_tag, "");
   testCase.verifyEqual(loaded.baseline_type, "rolling");
   testCase.verifyEqual(loaded.smbmodel_filter, "icemodel");
   testCase.verifyTrue(isdatetime(loaded.last_updated_utc));
   clear cleanup
end

function test_loadRegressionBaseline_prefers_explicit_path(testCase)
   % Explicit baseline paths should still load rolling files even when the
   % rolling selector has already been normalized to a blank tag.

   filepath = [tempname '.mat'];
   cleanup = onCleanup(@() deleteIfExists(filepath));

   RegressionBaseline = table( ...
      "icemodel_kanm_2016_solver2", ...
      1.5, ...
      'VariableNames', {'case_id', 'runoff_final'});
   save(filepath, 'RegressionBaseline');

   loaded = icemodel.test.helpers.loadRegressionBaseline("", ...
      "icemodel", filepath);

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

function deleteIfExists(filepath)
   %DELETEIFEXISTS Remove one file if it exists.

   if exist(filepath, 'file') == 2
      delete(filepath);
   end
end
