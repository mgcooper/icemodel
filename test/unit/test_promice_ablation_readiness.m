function tests = test_promice_ablation_readiness
   %TEST_PROMICE_ABLATION_READINESS Verify payload-derived cohort admission.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Install project dependencies and build one compact staged family tree.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   root = string(tempname);
   testCase.TestData.root = root;
   testCase.TestData.eval_root = fullfile(root, "eval");
   testCase.TestData.input_root = fullfile(root, "input");
   testCase.TestData.output_one = fullfile(root, "output-one");
   testCase.TestData.output_two = fullfile(root, "output-two");
   mkdir(testCase.TestData.eval_root)
   mkdir(testCase.TestData.input_root)
   writeFixtureTree(testCase.TestData.eval_root, ...
      testCase.TestData.input_root);
   testCase.TestData.report = ...
      icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=testCase.TestData.output_one);
end

function teardownOnce(~)
   % Remove synthetic payloads after every test has released loaded tables.
end

function test_payload_statuses_and_sensitivity_are_separate(testCase)
   % A gappy target remains evaluable while all four readiness lanes disagree.
   rows = testCase.TestData.report.rows;
   testCase.verifyEqual(height(rows), 5);

   ready = annualRecord(rows, "ready", 2019);
   testCase.verifyTrue(ready.ice_model_forcing_ready);
   testCase.verifyTrue(ready.evaluation_target_ready);
   testCase.verifyTrue(ready.snow_input_ready);
   testCase.verifyFalse(ready.production_snow_candidate_ready);
   testCase.verifyTrue(ready.admitted);
   testCase.verifyEqual(ready.requested_window_start, ...
      "2019-01-01 00:00:00");
   testCase.verifyNotEmpty(ready.forcing_readiness_artifact);
   testCase.verifyEqual(strlength(ready.forcing_readiness_sha256), 64);
   readiness_file = fullfile(testCase.TestData.root, ...
      ready.forcing_readiness_artifact);
   testCase.verifyEqual(ready.forcing_readiness_sha256, ...
      icemodel.verification.setup.fileSha256(readiness_file));
   testCase.verifyEqual(ready.gap_flag_count, 1);
   testCase.verifyEqual(ready.station_transition_flag_count, 1);
   testCase.verifyEqual(ready.step_detected_flag_count, 2);
   testCase.verifyEqual(ready.step_correctable_flag_count, 1);
   testCase.verifyEqual(ready.unresolved_step_flag_count, 2);
   testCase.verifyLessThan(ready.direct_target_count, ...
      ready.finite_target_count);
   testCase.verifyEqual(ready.snow_depth_min_m, 0.005, 'AbsTol', 1e-12);
   testCase.verifyEqual(ready.snow_depth_p50_m, 0.005, 'AbsTol', 1e-12);
   testCase.verifyGreaterThan(ready.snow_free_010mm_window_days, 30);
   testCase.verifyGreaterThan(ready.snow_free_050mm_window_days, 30);
   testCase.verifyEqual(ready.snow_free_100mm_window_days, ...
      ready.snow_free_050mm_window_days);
   window_start = icemodel.verification.setup.ensureUtc( ...
      ready.snow_free_window_start);
   window_end = icemodel.verification.setup.ensureUtc( ...
      ready.snow_free_window_end);
   gap_time = datetime(2019, 6, 1, 'TimeZone', 'UTC') ...
      + hours(floor(40 * 24 / 2) - 1);
   testCase.verifyLessThan(window_start, gap_time);
   testCase.verifyGreaterThan(window_end, gap_time);
   testCase.verifyGreaterThan(window_start, ...
      datetime(2019, 6, 2, 1, 0, 0, 'TimeZone', 'UTC'));
   testCase.verifyLessThan(window_end, ...
      datetime(2019, 7, 9, 23, 0, 0, 'TimeZone', 'UTC'));

   % ZACA proves a complete forcing product cannot stand in for an absent
   % ablation target, while the forcing-missing row proves the converse.
   zaca = annualRecord(rows, "zaca", 2019);
   testCase.verifyTrue(zaca.ice_model_forcing_ready);
   testCase.verifyFalse(zaca.target_available);
   testCase.verifyFalse(zaca.evaluation_target_ready);
   testCase.verifySubstring(zaca.evaluation_target_reason, ...
      "missing required variable");
   forcing_missing = annualRecord(rows, "forcing_missing", 2019);
   testCase.verifyFalse(forcing_missing.ice_model_forcing_ready);
   testCase.verifyTrue(forcing_missing.evaluation_target_ready);
   testCase.verifyFalse(forcing_missing.snow_input_ready);
   testCase.verifySubstring(forcing_missing.ice_model_forcing_reason, ...
      "producer manifest is missing");
   testCase.verifyEqual(forcing_missing.forcing_readiness_artifact, "");
   testCase.verifyEqual(forcing_missing.forcing_readiness_sha256, "");

   % A manifest and forcing payload that begin in June cannot redefine the
   % annual initialization request or be treated as full-year-ready.
   late = annualRecord(rows, "late_start", 2019);
   testCase.verifyEqual(late.requested_window_start, ...
      "2019-01-01 00:00:00");
   testCase.verifyFalse(late.evaluation_target_ready);
   testCase.verifySubstring(late.evaluation_target_reason, ...
      "case period does not contain requested annual window");
   testCase.verifyFalse(late.ice_model_forcing_ready);
   testCase.verifySubstring(late.ice_model_forcing_reason, ...
      "does not cover the requested window");

   % A missing observation remains represented rather than disappearing from
   % the canonical case/year inventory.
   missing = annualRecord(rows, "missing_obs", 2019);
   testCase.verifyFalse(missing.evaluation_target_ready);
   testCase.verifySubstring(missing.evaluation_target_reason, ...
      "observation artifact is missing");
end

function test_trace_snow_does_not_split_primary_summer_window(testCase)
   % Trace snow is allowed inside continuity but not counted as exposed ice.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "trace-snow");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   trace = 300:360;
   targets.data.snow_depth(trace) = 0.04;
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "trace-snow-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyTrue(ready.evaluation_target_ready);
   testCase.verifyGreaterThan(ready.snow_free_window_days, 30);
   testCase.verifyLessThan(ready.snow_free_direct_count, ...
      ready.direct_target_count);
end

function test_negative_snow_is_unknown_not_exposed(testCase)
   % An invalid negative posting cannot count as exposed-ice support.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "negative-snow");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   negative_row = 300;
   targets.data.snow_depth(negative_row) = -0.01;
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "negative-snow-output"));
   ready = annualRecord(report.rows, "ready", 2019);
   reference = annualRecord(testCase.TestData.report.rows, "ready", 2019);

   testCase.verifyTrue(ready.evaluation_target_ready);
   testCase.verifyEqual(ready.snow_free_window_start, ...
      reference.snow_free_window_start);
   testCase.verifyEqual(ready.snow_free_window_end, ...
      reference.snow_free_window_end);
   testCase.verifyEqual(ready.snow_free_direct_count, ...
      reference.snow_free_direct_count - 1);
end

function test_inventory_hashes_and_no_write_path_are_machine_readable(testCase)
   % Summary evidence is derived from rows and detects untracked staged cases.
   report = testCase.TestData.report;
   testCase.verifyEqual(report.summary.row_count, height(report.rows));
   testCase.verifyEqual(report.summary.admitted_count, 1);
   testCase.verifyTrue(ismember("extra", ...
      report.summary.inventory.extra_observation_cases));
   testCase.verifyTrue(ismember("missing_obs", ...
      report.summary.inventory.missing_observation_cases));
   testCase.verifyTrue( ...
      report.summary.reconstruction_source_audit.observation_write_path_verified);
   testCase.verifyGreaterThan( ...
      report.summary.reconstruction_source_audit.reconstruction_source_file_count, ...
      0);
   testCase.verifyGreaterThan( ...
      report.summary.reconstruction_source_audit.reconstruction_write_call_count, ...
      0);
   testCase.verifyGreaterThanOrEqual( ...
      report.summary.reconstruction_source_audit.runtime_evaluation_guard_call_count, ...
      2);
   testCase.verifyEqual( ...
      report.summary.reconstruction_source_audit.observation_write_match_count, 0);
   testCase.verifySubstring( ...
      report.summary.reconstruction_source_audit.evidence, ...
      "canonical runtime guards");

   decoded = jsondecode(fileread(report.files.json));
   testCase.verifyEqual(decoded.summary.row_count, height(report.rows));
   testCase.verifyEqual(string(decoded.csv_sha256), ...
      icemodel.verification.setup.fileSha256(report.files.csv));
   testCase.verifySubstring(string(decoded.observation_hash_scope), ...
      "forward baseline");
   decoded_ready = decoded.rows(strcmp({decoded.rows.case_id}, 'ready'));
   ready = annualRecord(report.rows, "ready", 2019);
   testCase.verifyEqual(string(decoded_ready.forcing_readiness_artifact), ...
      ready.forcing_readiness_artifact);
   testCase.verifyEqual(strlength(string( ...
      decoded_ready.forcing_readiness_sha256)), 64);
end

function test_outputs_are_byte_deterministic(testCase)
   % Repeating an explicit-root audit must preserve both portable artifacts.
   repeated = ...
      icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=testCase.TestData.output_two);
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(repeated.files.csv), ...
      icemodel.verification.setup.fileSha256( ...
      testCase.TestData.report.files.csv));
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(repeated.files.json), ...
      icemodel.verification.setup.fileSha256( ...
      testCase.TestData.report.files.json));
end

function test_missing_manifest_fails_before_writing(testCase)
   % An explicit empty evaluation root must fail with the public error id.
   empty_eval = fullfile(testCase.TestData.root, "empty-eval");
   mkdir(empty_eval)
   action = @() ...
      icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=empty_eval, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "unused"));
   testCase.verifyError(action, ...
      'icemodel:verification:promiceAblationReadiness:manifestMissing');
end

function test_shifted_hourly_coordinate_is_rejected(testCase)
   % Uniform one-hour spacing at :30 is not the native UTC observation grid.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "shifted-coordinate");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   targets.data.Time = targets.data.Time + minutes(30);
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "shifted-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyFalse(ready.evaluation_target_ready);
   testCase.verifySubstring(ready.evaluation_target_reason, ...
      "native whole-hour UTC coordinate");
end

function test_unzoned_observation_times_are_treated_as_utc(testCase)
   % Staged wall times without a zone use the shared UTC convention before
   % any requested-window filtering or coordinate-readiness assessment.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "unzoned-coordinate");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   row_times = targets.data.Time;
   row_times.TimeZone = '';
   targets.data.Time = row_times;
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "unzoned-output"));
   ready = annualRecord(report.rows, "ready", 2019);
   baseline = annualRecord(testCase.TestData.report.rows, "ready", 2019);

   testCase.verifyTrue(ready.evaluation_target_ready);
   testCase.verifyTrue(ready.admitted);
   testCase.verifyEqual(ready.observation_sample_count, ...
      baseline.observation_sample_count);
   testCase.verifyEqual(ready.snow_free_window_start, ...
      baseline.snow_free_window_start);
   testCase.verifyEqual(ready.snow_free_window_end, ...
      baseline.snow_free_window_end);
end

function test_non_utc_observation_times_are_converted_to_utc(testCase)
   % A staged zoned coordinate retains its instants while the readiness
   % writer converts its display zone to the native UTC comparison frame.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "non-utc-coordinate");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   row_times = targets.data.Time;
   row_times.TimeZone = 'America/New_York';
   targets.data.Time = row_times;
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "non-utc-output"));
   ready = annualRecord(report.rows, "ready", 2019);
   baseline = annualRecord(testCase.TestData.report.rows, "ready", 2019);

   testCase.verifyTrue(ready.evaluation_target_ready);
   testCase.verifyTrue(ready.admitted);
   testCase.verifyEqual(ready.observation_sample_count, ...
      baseline.observation_sample_count);
   testCase.verifyEqual(ready.snow_free_window_start, ...
      baseline.snow_free_window_start);
   testCase.verifyEqual(ready.snow_free_window_end, ...
      baseline.snow_free_window_end);
end

function test_malformed_observation_times_exclude_only_their_row(testCase)
   % A non-datetime timetable coordinate is a scientific row exclusion; it
   % must not abort the complete cohort CSV and JSON readiness audit.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "malformed-coordinate");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   targets.data.Properties.RowTimes = hours( ...
      (0:height(targets.data) - 1)');
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "malformed-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyEqual(height(report.rows), 5);
   testCase.verifyFalse(ready.evaluation_target_ready);
   testCase.verifyFalse(ready.admitted);
   testCase.verifyEqual(ready.evaluation_target_reason, ...
      "observation timestamps cannot be converted to a complete UTC " ...
      + "datetime coordinate");
   testCase.verifyTrue(isfile(report.files.csv));
   testCase.verifyTrue(isfile(report.files.json));
end

function test_unknown_datum_flag_splits_readiness_window(testCase)
   % An unknown station-transition flag splits the cumulative datum just like
   % an explicit transition and leaves no qualifying 30-day segment here.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "unknown-datum");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   middle = floor(height(targets.data) / 2);
   targets.data.surface_height_flag(middle + 10) = NaN;
   targets.data.station_transition_flag(middle + 10) = Inf;
   targets.data.step_detected_flag(middle + 10) = NaN;
   targets.data.step_correctable_flag(middle + 10) = Inf;
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "unknown-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyFalse(ready.evaluation_target_ready);
   testCase.verifyLessThan(ready.snow_free_window_days, 30);
   testCase.verifyEqual(ready.gap_flag_count, 1);
   testCase.verifyEqual(ready.station_transition_flag_count, 1);
   testCase.verifyEqual(ready.step_detected_flag_count, 2);
   testCase.verifyEqual(ready.step_correctable_flag_count, 1);
   testCase.verifyEqual(ready.unresolved_step_flag_count, 2);
end

function test_missing_surface_posting_remains_admissible(testCase)
   % Omitted ordinary postings stay inside the cumulative interval; elapsed
   % duration and direct exposed-ice endpoints come from the retained times.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "missing-posting");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   baseline = annualRecord(testCase.TestData.report.rows, "ready", 2019);
   expected_start = icemodel.verification.setup.ensureUtc( ...
      baseline.snow_free_window_start);
   expected_end = icemodel.verification.setup.ensureUtc( ...
      baseline.snow_free_window_end);
   middle = floor(height(targets.data) / 2);
   omitted = middle - 1:middle + 1;
   omitted_times = targets.data.Time(omitted);
   targets.data(omitted, :) = [];
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "missing-posting-output"));
   ready = annualRecord(report.rows, "ready", 2019);
   window_start = icemodel.verification.setup.ensureUtc( ...
      ready.snow_free_window_start);
   window_end = icemodel.verification.setup.ensureUtc( ...
      ready.snow_free_window_end);

   testCase.verifyTrue(ready.evaluation_target_ready);
   testCase.verifyTrue(ready.admitted);
   testCase.verifyEqual(window_start, expected_start);
   testCase.verifyEqual(window_end, expected_end);
   testCase.verifyEqual(ready.snow_free_window_days, ...
      days(window_end - window_start), 'AbsTol', 1e-12);
   testCase.verifyEqual(ready.snow_free_window_days, ...
      baseline.snow_free_window_days, 'AbsTol', 1e-12);
   testCase.verifyEqual(ready.snow_free_direct_count, ...
      baseline.snow_free_direct_count - 2);
   testCase.verifyEqual(ready.possible_support_count, ...
      baseline.possible_support_count);
   testCase.verifyEqual(ready.observation_sample_count, ...
      baseline.observation_sample_count - numel(omitted));
   testCase.verifyLessThan(window_start, min(omitted_times));
   testCase.verifyGreaterThan(window_end, max(omitted_times));

   % The retained endpoints must still satisfy the same direct exposed-ice
   % contract used by the production comparison, not merely bound a long gap.
   policy = icemodel.verification.namelists.promiceAblationReadiness();
   endpoints = ismember(targets.data.Time, [window_start; window_end]);
   support = targets.data{endpoints, ...
      cellstr(policy.support_flag_fields)};
   direct_flags = targets.data{endpoints, ...
      cellstr(policy.direct_zero_flag_fields)};
   testCase.verifyEqual(nnz(endpoints), 2);
   testCase.verifyTrue(all(isfinite(targets.data.ablation(endpoints))));
   testCase.verifyTrue(all(isfinite(support), 'all'));
   testCase.verifyTrue(all(direct_flags == 0, 'all'));
   testCase.verifyTrue(all(targets.data.snow_depth(endpoints) ...
      <= policy.ice_exposure_threshold_m));
end

function test_absent_target_table_metadata_excludes_only_that_case(testCase)
   % Empty optional timetable metadata arrays must create a row-level semantic
   % exclusion while the remaining canonical cohort is still audited/written.
   properties = ["VariableUnits", "VariableDescriptions"];
   expected = ["target units are not m", ...
      "target description does not identify geometric lowering"];
   for k = 1:numel(properties)
      [eval_root, observation_file] = copyEvaluationFixture( ...
         testCase, "missing-" + lower(properties(k)));
      saved = load(observation_file, 'targets');
      targets = saved.targets;
      if properties(k) == "VariableUnits"
         targets.data.Properties.VariableUnits = {};
      else
         targets.data.Properties.VariableDescriptions = {};
      end
      save(observation_file, 'targets')

      output_dir = fullfile(testCase.TestData.root, ...
         "missing-" + lower(properties(k)) + "-output");
      report = ...
         icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=eval_root, ...
         input_data_root=testCase.TestData.input_root, ...
         output_dir=output_dir);
      ready = annualRecord(report.rows, "ready", 2019);

      testCase.verifyEqual(height(report.rows), 5);
      testCase.verifyFalse(ready.evaluation_target_ready);
      testCase.verifyFalse(ready.admitted);
      testCase.verifySubstring(ready.evaluation_target_reason, expected(k));
      testCase.verifySubstring(ready.exclusion_reason, expected(k));
      testCase.verifyTrue(isfile(report.files.csv));
      testCase.verifyTrue(isfile(report.files.json));
   end
end

function test_nonnumeric_required_observation_excludes_only_that_case(testCase)
   % A malformed target or flag must produce one durable row-level exclusion,
   % not abort the remaining cohort or suppress its CSV/JSON evidence.
   variables = ["ablation", "surface_height_flag"];
   for variable = variables
      [eval_root, observation_file] = copyEvaluationFixture( ...
         testCase, "nonnumeric-" + variable);
      saved = load(observation_file, 'targets');
      targets = saved.targets;
      if variable == "ablation"
         targets.data.(variable) = string(targets.data.(variable));
      else
         targets.data.(variable) = categorical(targets.data.(variable));
      end
      save(observation_file, 'targets')

      output_dir = fullfile(testCase.TestData.root, ...
         "nonnumeric-" + variable + "-output");
      report = ...
         icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=eval_root, ...
         input_data_root=testCase.TestData.input_root, ...
         output_dir=output_dir);
      ready = annualRecord(report.rows, "ready", 2019);

      testCase.verifyEqual(height(report.rows), 5);
      testCase.verifyFalse(ready.evaluation_target_ready);
      testCase.verifyFalse(ready.admitted);
      testCase.verifySubstring(ready.evaluation_target_reason, ...
         "nonnumeric required variable(s): " + variable);
      testCase.verifyTrue(isfile(report.files.csv));
      testCase.verifyTrue(isfile(report.files.json));
   end
end

function test_unknown_snow_stays_inside_explicitly_bounded_window(testCase)
   % An explicit finite snow-present row splits continuity, while a later NaN
   % snow posting may remain inside a near-minimum cumulative endpoint window.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "unknown-snow");
   saved = load(observation_file, 'targets');
   targets = saved.targets;
   snow_break = 210;
   unknown_snow = 500;
   targets.data.snow_depth(snow_break) = 0.20;
   targets.data.snow_depth(unknown_snow) = NaN;
   expected_start = targets.data.Time(snow_break + 1);
   expected_end = targets.data.Time(end - 25);
   unknown_time = targets.data.Time(unknown_snow);
   save(observation_file, 'targets')

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, "unknown-snow-output"));
   ready = annualRecord(report.rows, "ready", 2019);
   window_start = icemodel.verification.setup.ensureUtc( ...
      ready.snow_free_window_start);
   window_end = icemodel.verification.setup.ensureUtc( ...
      ready.snow_free_window_end);

   testCase.verifyTrue(ready.evaluation_target_ready);
   testCase.verifyTrue(ready.admitted);
   testCase.verifyGreaterThan(ready.snow_free_window_days, 30);
   testCase.verifyLessThan(ready.snow_free_window_days, 31);
   testCase.verifyEqual(window_start, expected_start);
   testCase.verifyEqual(window_end, expected_end);
   testCase.verifyLessThan(window_start, unknown_time);
   testCase.verifyGreaterThan(window_end, unknown_time);
   endpoint_snow = targets.data.snow_depth( ...
      ismember(targets.data.Time, [window_start; window_end]));
   testCase.verifyTrue(all(isfinite(endpoint_snow)));
   testCase.verifyTrue(all(endpoint_snow <= 0.01));
end

function test_observation_station_and_site_identity_must_match_case(testCase)
   % Either crossed metadata identity invalidates the otherwise usable target.
   for variant = 1:2
      [eval_root, observation_file] = copyEvaluationFixture( ...
         testCase, "observation-identity-" + variant);
      saved = load(observation_file, 'targets');
      targets = saved.targets;
      if variant == 1
         targets.metadata.station = 'CROSSED';
      else
         targets.metadata.site_id = 'CROSSED';
      end
      save(observation_file, 'targets')

      report = icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=eval_root, ...
         input_data_root=testCase.TestData.input_root, ...
         output_dir=fullfile(testCase.TestData.root, ...
         "observation-identity-output-" + variant));
      ready = annualRecord(report.rows, "ready", 2019);

      testCase.verifyFalse(ready.evaluation_target_ready);
      testCase.verifyFalse(ready.admitted);
      testCase.verifySubstring(ready.evaluation_target_reason, ...
         "station/site identity does not match case");
   end
end

function test_requested_window_requires_canonical_period_containment(testCase)
   % A partial canonical period excludes rather than shrinking initialization.
   [eval_root, ~] = copyEvaluationFixture(testCase, "partial-case-period");
   manifest_file = fullfile(eval_root, "promice", "manifest.json");
   manifest = jsondecode(fileread(manifest_file));
   ready_index = find(string({manifest.cases.case_id}) == "ready");
   period_start = datetime(2019, 6, 2, 0, 0, 0, 'TimeZone', 'UTC');
   period_end = datetime(2019, 7, 9, 23, 0, 0, 'TimeZone', 'UTC');
   manifest.cases(ready_index).period = ...
      periodStruct(period_start, period_end);
   writeJson(manifest_file, manifest)

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, ...
      "partial-case-period-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyEqual(ready.requested_window_start, ...
      "2019-01-01 00:00:00");
   testCase.verifyEqual(ready.requested_window_end, ...
      "2019-12-31 23:00:00");
   testCase.verifyTrue(ready.ice_model_forcing_ready);
   testCase.verifyFalse(ready.evaluation_target_ready);
   testCase.verifyFalse(ready.admitted);
   testCase.verifySubstring(ready.evaluation_target_reason, ...
      "case period does not contain requested annual window");
end

function test_observation_path_outside_selected_root_is_rejected(testCase)
   % A canonical manifest must not hash or load an observation reached through
   % parent traversal outside the selected evaluation root.
   [eval_root, observation_file] = copyEvaluationFixture( ...
      testCase, "observation-path-escape");
   outside_file = fullfile(testCase.TestData.root, ...
      "outside-observations.mat");
   copyfile(observation_file, outside_file)
   manifest_file = fullfile(eval_root, "promice", "manifest.json");
   manifest = jsondecode(fileread(manifest_file));
   ready_index = find(string({manifest.cases.case_id}) == "ready");
   manifest.cases(ready_index).evaluation_file = ...
      char(fullfile("..", "..", "outside-observations.mat"));
   writeJson(manifest_file, manifest)

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      output_dir=fullfile(testCase.TestData.root, ...
      "observation-path-escape-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyFalse(ready.evaluation_target_ready);
   testCase.verifyFalse(ready.admitted);
   testCase.verifyEqual(ready.observation_artifact, "");
   testCase.verifyEqual(ready.observation_sha256, "");
   testCase.verifySubstring(ready.evaluation_target_reason, ...
      "escapes selected evaluation root");
end

function test_forcing_artifact_path_outside_selected_root_is_rejected(testCase)
   % A valid digest cannot authorize a producer path outside its selected root.
   [data_root, input_root] = copyForcingFixture( ...
      testCase, "forcing-path-escape");
   manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
      "plans", "ready-report-inputs.json");
   producer = jsondecode(fileread(manifest_file));
   filled_index = find(string({producer.artifacts.role}) == "filled");
   filled_file = fullfile(data_root, ...
      string(producer.artifacts(filled_index).path));
   outside_file = fullfile(testCase.TestData.root, "outside-filled.mat");
   copyfile(filled_file, outside_file)
   producer.artifacts(filled_index).path = char( ...
      fullfile("..", "outside-filled.mat"));
   producer.artifacts(filled_index).sha256 = char( ...
      icemodel.verification.setup.fileSha256(outside_file));
   writeJson(manifest_file, producer)

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=input_root, ...
      output_dir=fullfile(testCase.TestData.root, ...
      "forcing-path-escape-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyFalse(ready.ice_model_forcing_ready);
   testCase.verifyFalse(ready.admitted);
   testCase.verifySubstring(ready.ice_model_forcing_reason, ...
      "relative to its selected root");
end

function test_filled_payload_site_identity_must_match_case(testCase)
   % Hash-pinned bytes for another station must not admit the requested case.
   [data_root, input_root] = copyForcingFixture( ...
      testCase, "forcing-site-mismatch");
   manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
      "plans", "ready-report-inputs.json");
   producer = jsondecode(fileread(manifest_file));
   filled_index = find(string({producer.artifacts.role}) == "filled");
   filled_file = fullfile(data_root, ...
      string(producer.artifacts(filled_index).path));
   saved = load(filled_file, 'met', 'artifact_metadata');
   met = saved.met;
   met.Properties.UserData.site = 'crossed';
   saveFilledArtifactAndRepin( ...
      filled_file, manifest_file, producer, filled_index, met)

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=input_root, ...
      output_dir=fullfile(testCase.TestData.root, ...
      "forcing-site-mismatch-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyFalse(ready.ice_model_forcing_ready);
   testCase.verifyFalse(ready.admitted);
   testCase.verifySubstring(ready.ice_model_forcing_reason, ...
      "site identity does not match case");
end

function test_filled_payload_cadence_must_match_policy(testCase)
   % A regular artifact at another cadence must not pass the payload gate.
   [data_root, input_root] = copyForcingFixture( ...
      testCase, "forcing-cadence-mismatch");
   manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
      "plans", "ready-report-inputs.json");
   producer = jsondecode(fileread(manifest_file));
   filled_index = find(string({producer.artifacts.role}) == "filled");
   filled_file = fullfile(data_root, ...
      string(producer.artifacts(filled_index).path));
   saved = load(filled_file, 'met', 'artifact_metadata');
   met = saved.met(1:2:end, :);
   artifact_metadata = saved.artifact_metadata;
   save(filled_file, 'met', 'artifact_metadata')
   producer.artifacts(filled_index).sha256 = char( ...
      icemodel.verification.setup.fileSha256(filled_file));
   writeJson(manifest_file, producer)

   % The generic structural check accepts a complete regular native cadence;
   % the PROMICE policy gate below is what must reject the 30-minute payload.
   [native_ready, native_reason, ~, payload_cadence_seconds] = ...
      icemodel.verification.setup.metArtifactReadiness(filled_file);
   testCase.verifyTrue(native_ready, native_reason)
   testCase.verifyEqual(payload_cadence_seconds, 1800)

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=input_root, ...
      output_dir=fullfile(testCase.TestData.root, ...
      "forcing-cadence-mismatch-output"));
   ready = annualRecord(report.rows, "ready", 2019);
   policy = icemodel.verification.namelists.promiceAblationReadiness();

   testCase.verifyFalse(ready.ice_model_forcing_ready);
   testCase.verifyFalse(ready.admitted);
   testCase.verifySubstring(ready.ice_model_forcing_reason, ...
      "payload cadence");
   testCase.verifySubstring(ready.ice_model_forcing_reason, ...
      string(policy.forcing_cadence_seconds));
   testCase.verifySubstring(ready.ice_model_forcing_reason, "1800");
end

function test_forcing_readiness_is_requested_window_scoped(testCase)
   % A gap outside the requested year must not veto it, while a gap inside the
   % year must fail even when the producer readiness ledger says ready.

   variants = ["outside-window"; "inside-window"];
   for variant = variants'
      [data_root, input_root] = copyForcingFixture(testCase, variant);
      manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
         "plans", "ready-report-inputs.json");
      producer = jsondecode(fileread(manifest_file));
      filled_index = find(string({producer.artifacts.role}) == "filled");
      filled_file = fullfile(data_root, ...
         string(producer.artifacts(filled_index).path));
      saved = load(filled_file, 'met', 'artifact_metadata');
      met = saved.met;
      artifact_metadata = saved.artifact_metadata;

      % Mutate only the location of one missing forcing value; preserve cadence
      % and update the producer-pinned payload hash before consumer admission.
      if variant == "outside-window"
         prior = met(1, :);
         prior.Time = met.Time(1) - minutes(15);
         prior.tair(1) = NaN;
         prior.tair_provenance(1) = ...
            icemodel.forcing.reconstruct.provenanceCodes().missing;
         met = [prior; saved.met];
      else
         gap = floor(height(met) / 2);
         met.tair(gap) = NaN;
         met.tair_provenance(gap) = ...
            icemodel.forcing.reconstruct.provenanceCodes().missing;
      end
      save(filled_file, 'met', 'artifact_metadata')
      producer.artifacts(filled_index).sha256 = char( ...
         icemodel.verification.setup.fileSha256(filled_file));
      writeJson(manifest_file, producer)

      [whole_ready, ~, windows] = ...
         icemodel.verification.setup.metArtifactReadiness(filled_file);
      testCase.verifyFalse(whole_ready)
      testCase.verifyNotEmpty(windows)

      report = icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=testCase.TestData.eval_root, ...
         input_data_root=input_root, ...
         output_dir=fullfile(testCase.TestData.root, variant + "-output"));
      ready = annualRecord(report.rows, "ready", 2019);
      if variant == "outside-window"
         testCase.verifyTrue(ready.ice_model_forcing_ready)
         testCase.verifyTrue(ready.admitted)
      else
         testCase.verifyFalse(ready.ice_model_forcing_ready)
         testCase.verifyFalse(ready.admitted)
         testCase.verifySubstring(ready.ice_model_forcing_reason, ...
            "does not cover the requested window")
      end
   end
end

function test_acceptance_window_must_match_payload_support(testCase)
   % Producer window metadata cannot redefine the hash-pinned payload support.
   [data_root, input_root] = copyForcingFixture( ...
      testCase, "forcing-window-mismatch");
   manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
      "plans", "ready-report-inputs.json");
   producer = jsondecode(fileread(manifest_file));
   producer.acceptance_window.start = '2019-01-01 00:15:00';
   writeJson(manifest_file, producer)

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=input_root, ...
      output_dir=fullfile(testCase.TestData.root, ...
      "forcing-window-mismatch-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyFalse(ready.ice_model_forcing_ready);
   testCase.verifyFalse(ready.admitted);
   testCase.verifySubstring(ready.ice_model_forcing_reason, ...
      "does not match filled payload support");
end

function test_readiness_ledger_identity_must_match_site_and_window(testCase)
   % A crossed site or year range cannot borrow a producer readiness verdict.
   for variant = 1:3
      [data_root, input_root] = copyForcingFixture( ...
         testCase, "readiness-ledger-identity-" + variant);
      manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
         "plans", "ready-report-inputs.json");
      producer = jsondecode(fileread(manifest_file));
      readiness_index = find( ...
         string({producer.artifacts.role}) == "readiness");
      readiness_file = fullfile(data_root, ...
         string(producer.artifacts(readiness_index).path));
      ledger = readtable(readiness_file, TextType='string');
      if variant == 1
         ledger.site(:) = "crossed";
         expected_reason = "site identity does not match case";
      elseif variant == 2
         ledger.year(:) = 2018;
         expected_reason = "years do not match acceptance window";
      else
         ledger.site = [];
         expected_reason = "identity schema is incomplete";
      end
      writetable(ledger, readiness_file)
      producer.artifacts(readiness_index).sha256 = char( ...
         icemodel.verification.setup.fileSha256(readiness_file));
      writeJson(manifest_file, producer)

      report = icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=testCase.TestData.eval_root, ...
         input_data_root=input_root, ...
         output_dir=fullfile(testCase.TestData.root, ...
         "readiness-ledger-identity-output-" + variant));
      ready = annualRecord(report.rows, "ready", 2019);

      testCase.verifyFalse(ready.ice_model_forcing_ready);
      testCase.verifyFalse(ready.admitted);
      testCase.verifySubstring(ready.ice_model_forcing_reason, ...
         expected_reason);
   end
end

function test_invalid_acceptance_window_identity_is_rejected(testCase)
   % Missing, unparsable, and reversed producer windows are never identities.
   for variant = 1:3
      [data_root, input_root] = copyForcingFixture( ...
         testCase, "invalid-window-" + variant);
      manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
         "plans", "ready-report-inputs.json");
      producer = jsondecode(fileread(manifest_file));
      if variant == 1
         producer.acceptance_window = ...
            rmfield(producer.acceptance_window, 'end');
      elseif variant == 2
         producer.acceptance_window.start = 'not-a-timestamp';
      else
         producer.acceptance_window.start = '2020-01-01 00:00:00';
      end
      writeJson(manifest_file, producer)

      report = icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=testCase.TestData.eval_root, ...
         input_data_root=input_root, ...
         output_dir=fullfile(testCase.TestData.root, ...
         "invalid-window-output-" + variant));
      ready = annualRecord(report.rows, "ready", 2019);

      testCase.verifyFalse(ready.ice_model_forcing_ready);
      testCase.verifyFalse(ready.admitted);
      testCase.verifySubstring(ready.ice_model_forcing_reason, ...
         "acceptance-window identity is invalid");
   end
end

function test_incomplete_producer_artifact_schema_is_reported(testCase)
   % Missing producer artifact identity fields leave exact readiness provenance
   % empty and fail forcing readiness without aborting the cohort audit.
   [data_root, input_root] = copyForcingFixture( ...
      testCase, "incomplete-producer");
   manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
      "plans", "ready-report-inputs.json");
   producer = jsondecode(fileread(manifest_file));
   producer.artifacts = rmfield(producer.artifacts, 'sha256');
   writeJson(manifest_file, producer)

   report = icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=input_root, ...
      output_dir=fullfile(testCase.TestData.root, "incomplete-output"));
   ready = annualRecord(report.rows, "ready", 2019);

   testCase.verifyFalse(ready.ice_model_forcing_ready);
   testCase.verifySubstring(ready.ice_model_forcing_reason, ...
      "artifact schema is incomplete");
   testCase.verifyEqual(ready.forcing_readiness_artifact, "");
   testCase.verifyEqual(ready.forcing_readiness_sha256, "");
end

function test_policy_identity_mismatch_is_audited_without_aborting(testCase)
   % Missing and well-formed-but-wrong policy fingerprints must all fail the
   % same current-policy identity gate while the cohort audit is still saved.
   variants = ["stale", "missing", "forged"];
   for variant = variants
      [data_root, input_root] = copyForcingFixture( ...
         testCase, "policy-" + variant);
      [met_file, manifest_file, producer, filled_index] = ...
         loadFilledArtifact(data_root, "ready");
      saved = load(met_file, 'met');
      met = saved.met;
      if variant == "missing"
         met.Properties.UserData = rmfield( ...
            met.Properties.UserData, 'gapfill_policy_sha256');
      elseif variant == "stale"
         met.Properties.UserData.gapfill_policy_sha256 = ...
            string(repmat('a', 1, 64));
      else
         forged = char(icemodel.forcing.reconstruct.policySha256());
         if forged(end) == '0'
            forged(end) = '1';
         else
            forged(end) = '0';
         end
         met.Properties.UserData.gapfill_policy_sha256 = string(forged);
      end
      saveFilledArtifactAndRepin( ...
         met_file, manifest_file, producer, filled_index, met)

      report = ...
         icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=testCase.TestData.eval_root, ...
         input_data_root=input_root, ...
         output_dir=fullfile(testCase.TestData.root, ...
         "policy-" + variant + "-output"));
      ready = annualRecord(report.rows, "ready", 2019);

      testCase.verifyEqual(height(report.rows), 5);
      testCase.verifyFalse(ready.ice_model_forcing_ready);
      testCase.verifyFalse(ready.admitted);
      testCase.verifySubstring(ready.ice_model_forcing_reason, ...
         "canonical current promice_filled product");
      testCase.verifyTrue(isfile(report.files.csv));
      testCase.verifyTrue(isfile(report.files.json));
   end
end

function test_incomplete_filled_provenance_is_audited(testCase)
   % A missing registry or required channel ledger cannot be rescued by the
   % producer manifest hash or the separate readiness CSV.
   variants = ["registry", "channel"];
   for variant = variants
      [data_root, input_root] = copyForcingFixture( ...
         testCase, "incomplete-" + variant);
      [met_file, manifest_file, producer, filled_index] = ...
         loadFilledArtifact(data_root, "ready");
      saved = load(met_file, 'met');
      met = saved.met;
      if variant == "registry"
         met.Properties.UserData = rmfield( ...
            met.Properties.UserData, 'gapfill_registry');
      else
         met = removevars(met, 'lwd_provenance');
      end
      saveFilledArtifactAndRepin( ...
         met_file, manifest_file, producer, filled_index, met)

      report = ...
         icemodel.verification.setup.writePromiceAblationReadiness( ...
         evaluation_data_root=testCase.TestData.eval_root, ...
         input_data_root=input_root, ...
         output_dir=fullfile(testCase.TestData.root, ...
         "incomplete-" + variant + "-output"));
      ready = annualRecord(report.rows, "ready", 2019);

      testCase.verifyEqual(height(report.rows), 5);
      testCase.verifyFalse(ready.ice_model_forcing_ready);
      testCase.verifyFalse(ready.admitted);
      testCase.verifySubstring(ready.ice_model_forcing_reason, ...
         "complete canonical reconstruction provenance");
      testCase.verifyTrue(isfile(report.files.csv));
      testCase.verifyTrue(isfile(report.files.json));
   end
end

%% Fixture helpers
function writeFixtureTree(eval_root, input_root)
   %WRITEFIXTURETREE Create four cases spanning payload and readiness branches.
   family_root = fullfile(eval_root, "promice");
   mkdir(family_root)
   start_time = datetime(2019, 6, 1, 0, 0, 0, 'TimeZone', 'UTC');
   end_time = datetime(2019, 7, 10, 23, 0, 0, 'TimeZone', 'UTC');
   year_start = datetime(2019, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   year_end = datetime(2019, 12, 31, 23, 0, 0, 'TimeZone', 'UTC');

   cases = repmat(caseTemplate(), 5, 1);
   cases(1) = writeObservationCase(family_root, "ready", "READY", ...
      start_time, end_time, "ablation");
   cases(2) = writeObservationCase(family_root, "zaca", "ZAC_A", ...
      start_time, end_time, "surface_height");
   cases(3) = writeObservationCase(family_root, "forcing_missing", ...
      "FORCING_MISSING", start_time, end_time, "ablation");
   cases(4) = writeObservationCase(family_root, "late_start", ...
      "LATE_START", start_time, end_time, "ablation");
   for k = 1:3
      cases(k).period = periodStruct(year_start, year_end);
   end
   cases(5) = caseTemplate();
   cases(5).case_id = 'missing_obs';
   cases(5).site_id = 'MISSING_OBS';
   cases(5).surface_zone = 'ablation';
   cases(5).period = periodStruct(year_start, year_end);
   cases(5).evaluation_file = 'missing_obs/observations.mat';

   % The extra observations directory must be reported, never promoted into a
   % canonical ledger row.
   writeObservationCase(family_root, "extra", "EXTRA", ...
      start_time, end_time, "ablation");
   manifest = struct('dataset_family', 'promice', 'source_doi', '', ...
      'source_url', '', 'source_version', 'synthetic', ...
      'retrieval_date', '2026-08-03', 'cases', cases, ...
      'skipped', struct([]));
   writeJson(fullfile(family_root, "manifest.json"), manifest)

   % Ready and ZACA carry valid, producer-pinned forcing so their observation
   % verdicts can differ without conflating product availability.
   writeForcing(input_root, "ready", year_start, year_end)
   writeForcing(input_root, "zaca", year_start, year_end)
   writeForcing(input_root, "late_start", start_time, end_time)
end

function c = writeObservationCase( ...
      family_root, case_id, site_id, start_time, end_time, kind)
   %WRITEOBSERVATIONCASE Save one actual targets.data payload and manifest row.
   times = (start_time:hours(1):end_time).';
   n = numel(times);
   if kind == "ablation"
      ablation = linspace(0, 1, n).';
      snow_depth = repmat(0.005, n, 1);
      surface_height_flag = zeros(n, 1);
      station_transition_flag = zeros(n, 1);
      step_detected_flag = zeros(n, 1);
      step_correctable_flag = zeros(n, 1);
      middle = floor(n / 2);
      ablation(middle) = NaN;
      surface_height_flag(middle) = 1;
      step_detected_flag(25) = 1;
      step_detected_flag(26) = 1;
      step_correctable_flag(26) = 1;
      station_transition_flag(n - 24) = 1;
      data = timetable('RowTimes', times);
      data.ablation = ablation;
      data.snow_depth = snow_depth;
      data.surface_height_flag = surface_height_flag;
      data.station_transition_flag = station_transition_flag;
      data.step_detected_flag = step_detected_flag;
      data.step_correctable_flag = step_correctable_flag;
      data.Properties.VariableUnits = ...
         {'m', 'm', '1', '1', '1', '1'};
      data.Properties.VariableDescriptions = { ...
         'surface ablation (lowering) height', 'snow depth', ...
         'gap flag', 'station transition flag', ...
         'step detected flag', 'step correctable flag'};
   else
      surface_height = linspace(0, -1, n).';
      data = timetable('RowTimes', times);
      data.surface_height = surface_height;
      data.Properties.VariableUnits = {'m'};
      data.Properties.VariableDescriptions = {'surface height change'};
   end
   targets = struct('format', 'timeseries', 'data', data, ...
      'metadata', struct('source', 'synthetic', ...
      'source_family', 'promice', 'station', char(site_id), ...
      'site_id', char(site_id)));
   case_root = fullfile(family_root, case_id);
   if ~isfolder(case_root)
      mkdir(case_root)
   end
   save(fullfile(case_root, "observations.mat"), 'targets')

   c = caseTemplate();
   c.case_id = char(case_id);
   c.site_id = char(site_id);
   c.surface_zone = 'ablation';
   c.period = periodStruct(start_time, end_time);
   c.evaluation_file = char(fullfile(case_id, "observations.mat"));
end

function writeForcing(input_root, case_id, start_time, end_time)
   %WRITEFORCING Save one filled payload plus its final producer evidence.
   data_root = string(fileparts(input_root));
   met_dir = fullfile(input_root, 'met', 'promice_filled');
   qa_root = fullfile(data_root, 'preview', 'qa', 'gapfill');
   ledger_dir = fullfile(qa_root, 'ledger');
   plans_dir = fullfile(qa_root, 'plans');
   if ~isfolder(met_dir)
      mkdir(met_dir)
   end
   if ~isfolder(ledger_dir)
      mkdir(ledger_dir)
   end
   if ~isfolder(plans_dir)
      mkdir(plans_dir)
   end
   times = (start_time:minutes(15):end_time).';
   n = numel(times);
   met = timetable('RowTimes', times);
   met.tair = repmat(260, n, 1);
   met.swd = repmat(100, n, 1);
   met.lwd = repmat(250, n, 1);
   met.albedo = repmat(0.5, n, 1);
   met.wspd = repmat(5, n, 1);
   met.rh = repmat(80, n, 1);
   met.psfc = repmat(80000, n, 1);
   met.ppt = zeros(n, 1);
   met.rainf = zeros(n, 1);
   met.snowf = zeros(n, 1);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   channels = unique([ ...
      icemodel.forcing.reconstruct.icemodelRequiredChannels(), ...
      icemodel.forcing.helpers.precipitationVariables()], 'stable');
   for channel = channels
      met.(channel + "_provenance") = ...
         repmat(codes.observed, n, 1);
   end
   metadata = struct('site', char(case_id), ...
      'gapfill_registry', codes, 'gapfill_seed', 1, ...
      'gapfill_product', 'promice_filled', ...
      'gapfill_channels', ...
         icemodel.forcing.reconstruct.icemodelRequiredChannels(), ...
      'gapfill_engine_version', string(icemodel.internal.version()), ...
      'gapfill_policy_sha256', ...
         icemodel.forcing.reconstruct.policySha256(), ...
      'gapfill_donors', string.empty(1, 0));
   met.Properties.UserData = metadata;
   artifact_metadata = icemodel.forcing.helpers.artifactMetadata(met);
   met.Properties.UserData = artifact_metadata;
   met_name = "met_" + case_id + "_promice_filled_" ...
      + string(start_time, 'yyyyMMdd') + "_" ...
      + string(end_time, 'yyyyMMdd') + "_15m.mat";
   met_file = fullfile(met_dir, met_name);
   save(met_file, 'met', 'artifact_metadata')

   % The final producer ledger keeps IceModel and snow-input verdicts distinct.
   ledger = table(string(case_id), 2019, "ready", "", "ready", "", ...
      'VariableNames', {'site', 'year', 'verdict_icemodel', 'reason_icemodel', ...
      'verdict_snowmodel', 'reason_snowmodel'});
   ledger_name = case_id + "-readiness.csv";
   ledger_file = fullfile(ledger_dir, ledger_name);
   writetable(ledger, ledger_file)
   met_rel = "input/met/promice_filled/" + met_name;
   ledger_rel = "preview/qa/gapfill/ledger/" + ledger_name;
   artifacts = [artifactRecord("filled", met_rel, met_file), ...
      artifactRecord("readiness", ledger_rel, ledger_file)];
   producer = struct('site', char(case_id), ...
      'path_base', 'selected_data_root', 'artifacts', artifacts, ...
      'acceptance_window', struct( ...
      'start', char(formatTime(start_time)), ...
      'end', char(formatTime(end_time))));
   writeJson(fullfile(plans_dir, ...
      case_id + "-report-inputs.json"), producer)
end

function artifact = artifactRecord(role, path, pathname)
   %ARTIFACTRECORD Build one producer-pinned path and digest record.
   info = dir(pathname);
   artifact = struct('role', char(role), 'path', char(path), ...
      'bytes', info.bytes, 'sha256', char( ...
      icemodel.verification.setup.fileSha256(pathname)));
end

function c = caseTemplate()
   %CASETEMPLATE Keep the synthetic manifest case array homogeneous.
   c = struct('case_id', '', 'site_id', '', 'surface_zone', '', ...
      'period', struct('start', '', 'end', ''), 'evaluation_file', '', ...
      'colocation', struct());
end

function [eval_root, observation_file] = copyEvaluationFixture(testCase, name)
   %COPYEVALUATIONFIXTURE Isolate one readiness mutation from shared fixtures.
   eval_root = fullfile(testCase.TestData.root, "eval-" + name);
   copyfile(testCase.TestData.eval_root, eval_root)
   observation_file = fullfile( ...
      eval_root, "promice", "ready", "observations.mat");
end

function [data_root, input_root] = copyForcingFixture(testCase, name)
   %COPYFORCINGFIXTURE Isolate one producer mutation from shared fixtures.
   data_root = fullfile(testCase.TestData.root, name + "-producer-root");
   input_root = fullfile(data_root, "input");
   copyfile(testCase.TestData.input_root, input_root)
   copyfile(fullfile(testCase.TestData.root, "preview"), ...
      fullfile(data_root, "preview"))
end

function [met_file, manifest_file, producer, filled_index] = ...
      loadFilledArtifact(data_root, case_id)
   %LOADFILLEDARTIFACT Resolve one producer-pinned filled fixture for mutation.
   manifest_file = fullfile(data_root, "preview", "qa", "gapfill", ...
      "plans", case_id + "-report-inputs.json");
   producer = jsondecode(fileread(manifest_file));
   filled_index = find(string({producer.artifacts.role}) == "filled");
   met_file = fullfile(data_root, ...
      string(producer.artifacts(filled_index).path));
end

function saveFilledArtifactAndRepin( ...
      met_file, manifest_file, producer, filled_index, met)
   %SAVEFILLEDARTIFACTANDREPIN Keep both metadata copies and manifest coherent.
   artifact_metadata = met.Properties.UserData;
   save(met_file, 'met', 'artifact_metadata')
   producer.artifacts(filled_index) = artifactRecord( ...
      "filled", string(producer.artifacts(filled_index).path), met_file);
   writeJson(manifest_file, producer)
end

function period = periodStruct(first, last)
   %PERIODSTRUCT Encode one manifest window at second resolution.
   period = struct('start', char(formatTime(first)), ...
      'end', char(formatTime(last)));
end

function row = annualRecord(rows, case_id, year_value)
   %ANNUALRECORD Select one unique case-year from the returned public table.
   keep = rows.case_id == case_id & rows.year == year_value;
   assert(nnz(keep) == 1)
   row = rows(keep, :);
end

function text = formatTime(value)
   %FORMATTIME Format one manifest timestamp explicitly in UTC.
   value.TimeZone = 'UTC';
   text = string(value, 'yyyy-MM-dd HH:mm:ss');
end

function writeJson(pathname, value)
   %WRITEJSON Write one synthetic portable JSON artifact.
   fid = fopen(pathname, 'w', 'n', 'UTF-8');
   assert(fid >= 0)
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(value, PrettyPrint=true));
end
