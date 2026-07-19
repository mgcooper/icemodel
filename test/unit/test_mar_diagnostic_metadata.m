function tests = test_mar_diagnostic_metadata
   %TEST_MAR_DIAGNOSTIC_METADATA Verify optional MAR diagnostic provenance.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Install the repository namespace for direct focused execution.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.bootstrap_cleanup = cleanup;
end

function test_full_contract_validates_hourly_melt(testCase)
   % SUH, SU, RZ, and daily ME together form the complete optional schema.
   times = hourlyTimes(48);
   melt_daily_rate = [repmat(0.001, 24, 1); repmat(0.002, 24, 1)];
   T = timetable(times, melt_daily_rate, ...
      linspace(-2e-4, 2e-4, 48)', ...
      [repmat(-1e-5, 24, 1); repmat(2e-5, 24, 1)], ...
      [repmat(1e-4, 24, 1); repmat(2e-4, 24, 1)], ...
      VariableNames=["melt", "subl", "subl_evap", ...
      "refreeze_deposition"]);

   metadata = icemodel.forcing.helpers.marDiagnosticMetadata( ...
      T, melt_daily_rate, sector=1);

   testCase.verifyEqual(string(metadata.mar_diagnostic_status), "applied");
   testCase.verifyEqual(metadata.mar_diagnostic_channels, ...
      ["subl"; "subl_evap"; "refreeze_deposition"]);
   testCase.verifyEqual(metadata.mar_diagnostic_native_variables, ...
      ["SUH"; "SU"; "RZ"; "ME"]);
   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_melt_validation_status), "validated");
   testCase.verifyEqual(metadata.mar_diagnostic_melt_day_status, ...
      uint8([1; 1]));
   testCase.verifyEqual(metadata.mar_diagnostic_melt_daily_reference_mwe, ...
      [0.024; 0.048], 'AbsTol', 1e-12);
   testCase.verifyEqual(metadata.mar_diagnostic_melt_residual_mwe_day, ...
      [0; 0], 'AbsTol', 1e-12);
   testCase.verifyEqual(string(metadata.mar_diagnostic_su_sector_name), ...
      "permanent_ice");
   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_refreeze_deposition_sign), ...
      "native_signed_combined_term_preserved_not_pure_refreeze");
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_negative_day_count, 0);
   testCase.verifyTrue(isnan( ...
      metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h));
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_material_negative_threshold_mwe_h, ...
      1e-8);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_material_negative_day_count, 0);
   testCase.verifyTrue(isnan( ...
      metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h));
end

function test_mismatch_and_partial_day_are_distinguished(testCase)
   % A complete mismatch fails validation while a one-row next day stays
   % explicitly unverified rather than being treated as a daily aggregate.
   times = hourlyTimes(25);
   daily_rate = repmat(0.001, 25, 1);
   melt = daily_rate;
   melt(1:24) = 0.002;
   T = timetable(times, melt, VariableNames="melt");

   metadata = icemodel.forcing.helpers.marDiagnosticMetadata( ...
      T, daily_rate, sector=2);

   testCase.verifyEqual(string(metadata.mar_diagnostic_status), "partial");
   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_melt_validation_status), "mismatch");
   testCase.verifyEqual(metadata.mar_diagnostic_melt_day_status, ...
      uint8([2; 3]));
   testCase.verifyEqual(metadata.mar_diagnostic_melt_daily_reference_mwe(1), ...
      0.024, 'AbsTol', 1e-12);
   testCase.verifyTrue(isnan( ...
      metadata.mar_diagnostic_melt_daily_reference_mwe(2)));
   testCase.verifyEqual(string(metadata.mar_diagnostic_su_sector_name), ...
      "tundra");
end

function test_reduced_and_partial_sources_remain_nonfatal(testCase)
   % No optional fields is a valid reduced schema; one SUH field is partial.
   times = hourlyTimes(24);
   T = timetable(times, zeros(24, 1), VariableNames="melt");
   reduced = icemodel.forcing.helpers.marDiagnosticMetadata(T, []);
   testCase.verifyEqual(string(reduced.mar_diagnostic_status), ...
      "not_available");
   testCase.verifyEqual(string( ...
      reduced.mar_diagnostic_melt_validation_status), "not_available");
   testCase.verifyEmpty(reduced.mar_diagnostic_channels);
   testCase.verifyEmpty(reduced.mar_diagnostic_melt_day_status);
   testCase.verifyEqual( ...
      reduced.mar_diagnostic_refreeze_negative_day_count, 0);
   testCase.verifyTrue(isnan( ...
      reduced.mar_diagnostic_refreeze_negative_minimum_mwe_h));
   testCase.verifyEqual( ...
      reduced.mar_diagnostic_refreeze_material_negative_day_count, 0);
   testCase.verifyTrue(isnan( ...
      reduced.mar_diagnostic_refreeze_material_negative_minimum_mwe_h));

   T.subl = linspace(-1e-4, 1e-4, 24)';
   partial = icemodel.forcing.helpers.marDiagnosticMetadata(T, []);
   testCase.verifyEqual(string(partial.mar_diagnostic_status), "partial");
   testCase.verifyEqual(partial.mar_diagnostic_channels, "subl");
   testCase.verifyEqual(partial.mar_diagnostic_native_variables, "SUH");
end

function test_signed_refreeze_metadata_counts_days_and_removes_legacy(testCase)
   % Signed native RZ provenance preserves strict roundoff-scale events while
   % separately classifying material UTC days and removing the old tolerance.
   times = hourlyTimes(72);
   values = 1e-4 * ones(72, 1);
   values(1:24) = -2e-4;
   values(25) = -1e-10;
   values(49:71) = -3e-4;
   values(72) = complex(-1, 1);
   T = timetable(times, values, VariableNames="refreeze_deposition");
   legacy = struct( ...
      'mar_diagnostic_refreeze_negative_tolerance_mwe_h', 1e-8);

   metadata = icemodel.forcing.helpers.marRefreezeMetadata(T, legacy);

   testCase.verifyFalse(isfield(metadata, ...
      'mar_diagnostic_refreeze_negative_tolerance_mwe_h'));
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_negative_day_count, 3);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h, -3e-4);
   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_refreeze_negative_statistics_basis), ...
      "distinct_utc_days_and_minimum_finite_negative_artifact_value");
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_material_negative_threshold_mwe_h, ...
      1e-8);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_material_negative_day_count, 2);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h, ...
      -3e-4);
   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_refreeze_material_negative_statistics_basis), ...
      "distinct_utc_days_and_minimum_below_declared_reporting_threshold");
end

function test_signed_refreeze_metadata_matches_three_source_cases(testCase)
   % Freeze the exact all-43 sweep expectations for KAN_U, SCO_U, and TAS_A.
   case_ids = ["kanu", "scou", "tasa"];
   negative_days = [3, 4, 1];
   minima = [-1.46520137786865e-5, ...
      -7.63545682032903e-6, -5.54515173037847e-6];
   for k = 1:numel(case_ids)
      times = hourlyTimes(24 * negative_days(k));
      values = -2e-8 * ones(numel(times), 1);
      values(1:24) = minima(k);
      T = timetable(times, values, VariableNames="refreeze_deposition");

      metadata = icemodel.forcing.helpers.marRefreezeMetadata(T);

      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_negative_day_count, ...
         negative_days(k), case_ids(k));
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h, ...
         minima(k), case_ids(k));
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_material_negative_day_count, ...
         negative_days(k), case_ids(k));
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h, ...
         minima(k), case_ids(k));
   end
end

function test_signed_refreeze_metadata_counts_utc_days(testCase)
   % UTC grouping must not depend on a timetable's current display time zone.
   times = datetime(2020, 1, 1, 18, 0, 0, ...
      'TimeZone', 'America/Los_Angeles') + hours([0; 7; 8]);
   T = timetable(times, -1e-4 * ones(3, 1), ...
      VariableNames="refreeze_deposition");

   metadata = icemodel.forcing.helpers.marRefreezeMetadata(T);

   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_negative_day_count, 1);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h, -1e-4);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_material_negative_day_count, 1);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h, ...
      -1e-4);
end

function test_daily_melt_without_hourly_melt_is_unverified(testCase)
   % Native ME availability alone cannot validate a missing MEH public field.
   times = hourlyTimes(24);
   T = timetable(times, zeros(24, 1), VariableNames="subl");
   metadata = icemodel.forcing.helpers.marDiagnosticMetadata( ...
      T, repmat(0.001, 24, 1));

   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_melt_validation_status), "unverified");
   testCase.verifyEqual(metadata.mar_diagnostic_melt_day_status, uint8(3));
   testCase.verifyEqual(metadata.mar_diagnostic_melt_daily_reference_mwe, ...
      0.024, 'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(metadata.mar_diagnostic_melt_residual_mwe_day));
end

function test_missing_daily_or_hourly_support_remains_unverified(testCase)
   % A missing native-ME day and a missing hourly-MEH day exercise distinct
   % unverified paths without manufacturing a comparison residual.
   times = hourlyTimes(48);
   daily_rate = [nan(24, 1); repmat(0.001, 24, 1)];
   melt = repmat(0.001, 48, 1);
   melt(25) = NaN;
   T = timetable(times, melt, VariableNames="melt");
   metadata = icemodel.forcing.helpers.marDiagnosticMetadata(T, daily_rate);

   testCase.verifyEqual(metadata.mar_diagnostic_melt_day_status, ...
      uint8([3; 3]));
   testCase.verifyTrue(isnan( ...
      metadata.mar_diagnostic_melt_daily_reference_mwe(1)));
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_melt_daily_reference_mwe(2), 0.024, ...
      'AbsTol', 1e-12);
   testCase.verifyTrue(all(isnan( ...
      metadata.mar_diagnostic_melt_residual_mwe_day)));
end

function test_rejects_daily_reference_with_wrong_height(testCase)
   % A daily-ME support vector must align one-for-one with the hourly table.
   times = hourlyTimes(24);
   T = timetable(times, zeros(24, 1), VariableNames="melt");
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.marDiagnosticMetadata(T, zeros(23, 1)), ...
      'icemodel:forcing:marDiagnosticMetadata:heightMismatch');
end

function times = hourlyTimes(count)
   %HOURLYTIMES Return a deterministic UTC interval-start fixture axis.
   start = datetime(2000, 1, 1, 0, 0, 0, TimeZone="UTC");
   times = start + hours((0:count - 1)');
end
