function tests = test_baseline_contracts
%TEST_BASELINE_CONTRACTS Verify baseline-selection and legacy-load helpers.
   tests = functiontests(localfunctions);
end

function test_resolveBaselineSelector_handles_rolling_and_release(testCase)

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector("rolling");
   testCase.verifyEqual(baseline_type, "rolling");
   testCase.verifyEqual(baseline_tag, "");

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector("v1.1");
   testCase.verifyEqual(baseline_type, "release");
   testCase.verifyEqual(baseline_tag, "v1.1");
end

function test_loadRegressionBaseline_returns_empty_for_blank_selector(testCase)

   baseline = icemodel.test.helpers.loadRegressionBaseline("", "all");
   testCase.verifyTrue(istable(baseline));
   testCase.verifyEmpty(baseline);
end

function test_loadRegressionBaseline_normalizes_legacy_schema(testCase)

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

function deleteIfExists(filepath)

   if exist(filepath, 'file') == 2
      delete(filepath);
   end
end
