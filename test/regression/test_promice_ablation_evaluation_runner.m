function tests = test_promice_ablation_evaluation_runner
   %TEST_PROMICE_ABLATION_EVALUATION_RUNNER Verify the cohort runner.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   %SETUPONCE Build a paired synthetic PROMICE data tree for all focused tests.
   root = string(tempname);
   eval_root = fullfile(root, 'eval');
   input_root = fullfile(root, 'input');
   artifact_root = fullfile(root, 'artifacts');
   mkdir(eval_root)
   mkdir(input_root)
   writeFixtureTree(eval_root, input_root)
   testCase.TestData.root = root;
   testCase.TestData.eval_root = eval_root;
   testCase.TestData.input_root = input_root;
   testCase.TestData.artifact_root = artifact_root;
end

function teardownOnce(testCase)
   %TEARDOWNONCE Remove only the test-owned synthetic tree.
   if isfolder(testCase.TestData.root)
      rmdir(testCase.TestData.root, 's')
   end
end

function test_policy_predeclares_runtime_and_nested_windows(testCase)
   % Runtime cadence and cumulative windows must be fixed before outcomes exist.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   testCase.verifyEqual(policy.model_substep_seconds, 900)
   testCase.verifyEqual(policy.nested_window_days, [30, 60, 90])
   testCase.verifyEqual(policy.endpoint_perturbation_days, [1, 3, 7])
   testCase.verifyEqual(policy.evaluation_season_start_month_day, [6, 1])
   testCase.verifyEqual(policy.evaluation_season_end_month_day, [10, 1])
end

function test_default_is_readiness_only_and_transient(testCase)
   % Empty case selection must never launch the model or leave artifacts behind.
   results = run_promice_ablation_evaluation( ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      model_provider=@unexpectedProvider);

   % The admitted row remains visible but not selected; the missing-forcing row
   % retains its readiness exclusion rather than disappearing from the cohort.
   kanm = results.summary.case_id == "kanm";
   kanl = results.summary.case_id == "kanl";
   testCase.verifyEqual(results.summary.status(kanm), "not_selected")
   testCase.verifyEqual(results.summary.status(kanl), "excluded")
   testCase.verifyTrue(isnan( ...
      results.summary.window_valid_sample_count(kanm)))
   testCase.verifyTrue(isnan(results.summary.coverage_fraction(kanm)))
   testCase.verifyTrue(isnan( ...
      results.summary.window_valid_sample_count(kanl)))
   testCase.verifyGreaterThan( ...
      results.summary.window_possible_sample_count(kanl), 0)
   testCase.verifyTrue(isnan(results.summary.coverage_fraction(kanl)))
   testCase.verifyTrue(all(isnan( ...
      results.summary.eligible_sample_count(kanm | kanl))))
   testCase.verifyEqual(results.paths.run_dir, "")
   testCase.verifyEqual(results.readiness.files.csv, "")

   % A readiness-only run executes no model, so it stamps no physics
   % fingerprint. Resolving one runs icemodel.setopts, which asserts that the
   % workspace exists, and this run reads its data roots from the caller.
   testCase.verifyFalse(isfield(results, 'physics_fingerprint'))
end

function test_selected_run_pins_filled_forcing_and_half_open_boundary(testCase)
   % One selected KAN row must use promice_filled and exclude the t1 interval.
   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      artifact_root=testCase.TestData.artifact_root, ...
      run_name="20260804-010203", ...
      write_artifacts=true, model_provider=@selectedModelProvider);

   % A run that executes a case stamps the default physics it ran under, so a
   % report built later can compare it against the code of the day.
   testCase.verifyTrue(isfield(results, 'physics_fingerprint'))
   testCase.verifyTrue( ...
      icemodel.verification.helpers.isPhysicsFingerprint( ...
      results.physics_fingerprint))
   testCase.verifyEqual(results.physics_fingerprint.opts_sha256, ...
      icemodel.verification.helpers.physicsFingerprint().opts_sha256)

   % All readiness rows remain in the flat summary while only KAN_M executes.
   kanm = results.summary.case_id == "kanm";
   kanl = results.summary.case_id == "kanl";
   testCase.verifyEqual(results.summary.status(kanm), "completed", ...
      results.summary.reason(kanm))
   testCase.verifyEqual(results.summary.status(kanl), "excluded")
   testCase.verifyEqual(results.summary.classification(kanm), ...
      "within_materiality")

   % The final row is an October 1 state checkpoint, not a future interval.
   result = results.site_year_results(kanm);
   testCase.verifyEqual(string( ...
      result.readiness_row.forcing_producer_manifest), ...
      "preview/qa/gapfill/plans/kanm-report-inputs.json")
   testCase.verifyEqual(string( ...
      result.readiness_row.forcing_readiness_artifact), ...
      "preview/qa/gapfill/ledger/kanm-readiness.csv")
   testCase.verifyNotEmpty(string( ...
      result.readiness_row.forcing_producer_manifest_sha256))
   testCase.verifyNotEmpty(string( ...
      result.readiness_row.forcing_readiness_sha256))
   model = result.model;
   t1 = parseTime(result.readiness_row.snow_free_window_end);
   t0 = parseTime(result.readiness_row.snow_free_window_start);
   display_start = datetime(2019, 6, 1, 'TimeZone', 'UTC');
   display_end = datetime(2019, 10, 1, 'TimeZone', 'UTC');
   initialization_start = parseTime( ...
      result.readiness_row.requested_window_start);
   testCase.verifyEqual(model.Time(1), display_start)
   testCase.verifyEqual(model.Time(end), display_end)
   testCase.verifyEqual( ...
      result.model_metadata.initialization_start, initialization_start)
   testCase.verifyEqual(result.model_metadata.evaluation_start, t0)
   testCase.verifyEqual(result.model_metadata.display_start, display_start)
   testCase.verifyEqual(result.model_metadata.display_end, display_end)
   testCase.verifyEqual(result.model_metadata.initialization_policy, ...
      "readiness requested_window_start")
   testCase.verifyEqual(initialization_start, ...
      datetime(2019, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'))
   testCase.verifyLessThan(initialization_start, t0)
   testCase.verifyEqual( ...
      result.model_metadata.run_end_inclusive, display_end - minutes(15))
   testCase.verifyFalse(result.model_metadata.future_interval_executed)
   sum_fields = icemodel.namelists.budgetoutputs('sum');
   testCase.verifyEqual(model{end, sum_fields}, ...
      zeros(1, numel(sum_fields)), AbsTol=0)
   testCase.verifyEqual(model.mass_budget_solid_start_mwe(end), ...
      model.mass_budget_solid_end_mwe(end - 1), AbsTol=0)

   % Saved observations retain the compact display season and policy fields;
   % annual initialization rows and unused payload channels stay out.
   saved_observations = result.observations.data;
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   testCase.verifyEqual(saved_observations.Time(1), display_start)
   testCase.verifyEqual(saved_observations.Time(end), display_end)
   testCase.verifyEqual(string(saved_observations.Properties.VariableNames), ...
      policy.required_observation_fields)
   testCase.verifyFalse(ismember( ...
      "unused_payload", saved_observations.Properties.VariableNames))
   expected_seasonal = ["observation_lowering_m", ...
      "observation_lower_mwe", "observation_upper_mwe", ...
      "observation_reference_mwe", ...
      "model_melt_mwe", "model_runoff_mwe", "model_freeze_mwe", ...
      "model_net_solid_loss_mwe", "model_layer_change_mwe", ...
      "model_surface_mass_loss_mwe", "model_ablation_proxy_mwe", ...
      "snow_depth_m", "ice_exposed", ...
      "direct_observation", "evaluation_window"];
   testCase.verifyEqual(string(result.seasonal.Properties.VariableNames), ...
      expected_seasonal)
   testCase.verifyEqual(result.seasonal.Time(1), display_start)
   testCase.verifyEqual(result.seasonal.Time(end), display_end)
   testCase.verifyTrue(all(isnan(result.seasonal.model_melt_mwe( ...
      ~result.seasonal.ice_exposed))))
   visible = find(isfinite(result.seasonal.observation_lowering_m));
   testCase.verifyGreaterThan(numel(visible), 1)
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   density_values = result.seasonal.observation_lowering_m(visible) ...
      .* policy.effective_density_kg_m3 ./ ro_liq;
   testCase.verifyEqual(result.seasonal.observation_lower_mwe(visible), ...
      min(density_values, [], 2), AbsTol=1e-14)
   testCase.verifyEqual(result.seasonal.observation_upper_mwe(visible), ...
      max(density_values, [], 2), AbsTol=1e-14)
   expected_step = (ro_ice / ro_liq) ...
      / hours(display_end - display_start);
   testCase.verifyEqual( ...
      result.seasonal.model_melt_mwe(visible(2)), expected_step, ...
      AbsTol=1e-14)
   testCase.verifyEqual( ...
      result.seasonal.model_net_solid_loss_mwe(visible(2)), ...
      expected_step, AbsTol=1e-14)
   testCase.verifyEqual( ...
      result.seasonal.model_surface_mass_loss_mwe(visible(2)), 0.01, ...
      AbsTol=0)

   % The snow-free denominator is distinct from full-year readiness support.
   summary = results.summary(kanm, :);
   expected_possible = floor(hours(t1 - t0)) + 1;
   testCase.verifyEqual(summary.window_possible_sample_count, expected_possible)
   testCase.verifyEqual(summary.window_direct_sample_count, expected_possible)
   testCase.verifyEqual(summary.window_valid_sample_count, expected_possible)
   testCase.verifyEqual(summary.coverage_fraction, 1, AbsTol=0)
   testCase.verifySubstring(summary.coverage_denominator, ...
      "inclusive [t0,t1]")
   testCase.verifyLessThan( ...
      summary.window_possible_sample_count, summary.possible_support_count)

   % Fixed endpoint diagnostics never shorten to hide missing support.
   nested = result.nested_windows;
   testCase.verifyTrue(nested.available(nested.window_label == "30d"))
   testCase.verifyFalse(nested.available(nested.window_label == "60d"))
   testCase.verifyFalse(nested.available(nested.window_label == "90d"))
   testCase.verifyTrue(nested.available(nested.window_label == "longest"))

   % Every endpoint sensitivity reuses the same model run and retains both the
   % requested and exact actual endpoints in a stable machine-readable schema.
   perturbations = result.endpoint_perturbations;
   perturbation_days = policy.endpoint_perturbation_days(:);
   expected_labels = string([ ...
      compose('start_plus_%gd', perturbation_days); ...
      compose('end_minus_%gd', perturbation_days)]);
   expected_axes = [repmat("start", numel(perturbation_days), 1); ...
      repmat("end", numel(perturbation_days), 1)];
   expected_offsets = [perturbation_days; -perturbation_days];
   testCase.verifyEqual(perturbations.perturbation_label, expected_labels)
   testCase.verifyEqual(perturbations.perturbation_axis, expected_axes)
   testCase.verifyEqual( ...
      perturbations.perturbation_days_signed, expected_offsets)
   testCase.verifyTrue(all(perturbations.available))
   testCase.verifyEqual(perturbations.actual_window_start, ...
      perturbations.requested_window_start)
   testCase.verifyEqual(perturbations.actual_window_end, ...
      perturbations.requested_window_end)
   testCase.verifyTrue(all(isfinite( ...
      perturbations.observation_sensitivity_min_mwe)))
   testCase.verifyTrue(all( ...
      perturbations.observation_sensitivity_min_mwe ...
      <= perturbations.observation_sensitivity_max_mwe))
   testCase.verifyEqual(results.endpoint_perturbations, perturbations)
   testCase.verifyTrue(isfile(results.paths.summary_csv))
   testCase.verifyTrue(isfile(results.paths.nested_windows_csv))
   testCase.verifyTrue(isfile(results.paths.endpoint_perturbations_csv))
   testCase.verifyTrue(isfile(results.paths.results_mat))
   artifact = load(results.paths.results_mat, 'results');
   testCase.verifyEqual(numel(results.site_year_results), ...
      height(results.summary))
   testCase.verifyEqual(numel(artifact.results.site_year_results), 1)
   testCase.verifyEqual( ...
      string(artifact.results.site_year_results.status), "completed")
   testCase.verifyEqual(artifact.results.summary, results.summary)
   testCase.verifyEqual(artifact.results.readiness, results.readiness)
   testCase.verifyEqual(artifact.results.nested_windows, ...
      results.nested_windows)
   testCase.verifyEqual(artifact.results.endpoint_perturbations, ...
      results.endpoint_perturbations)

   % The report must consume compact site payloads without positional alignment
   % to the full readiness and summary ledgers.
   report_dir = fullfile(testCase.TestData.root, 'compact-report');
   report_file = ...
      icemodel.verification.report.buildAblationEvaluationReport( ...
      results.paths.results_mat, render=false, output_dir=report_dir);
   testCase.verifyEqual(report_file, string(fullfile(report_dir, ...
      'promice-ablation-evaluation-report.html')))
   testCase.verifyTrue(isfile(fullfile(report_dir, ...
      'promice-ablation-evaluation-report.qmd')))
end

function test_a_stale_global_path_does_not_stop_an_explicit_root_run(testCase)
   % physicsFingerprint resolves its reference through icemodel.setopts,
   % which asserts ICEMODEL_INPUT_PATH exists. A cohort that supplies
   % explicit roots reaches its forcing through the manifest root instead, so
   % an unset or stale global path says nothing about whether its cases can
   % run. The runner guards the stamp, and the run must still complete.

   restore = onCleanup(@() setenv('ICEMODEL_INPUT_PATH', ...
      getenv('ICEMODEL_INPUT_PATH')));
   saved_path = getenv('ICEMODEL_INPUT_PATH');
   cleaner = onCleanup(@() setenv('ICEMODEL_INPUT_PATH', saved_path));
   setenv('ICEMODEL_INPUT_PATH', ...
      fullfile(tempdir, 'icemodel-no-such-workspace'));

   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      artifact_root=testCase.TestData.artifact_root, ...
      run_name="20260804-010204", ...
      write_artifacts=true, model_provider=@selectedModelProvider);

   % The case still executes.
   kanm = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(kanm), "completed", ...
      'an explicit-root run must not depend on the global path');

   % The stamp is absent rather than the run failing. The report reads an
   % absent stamp without warning.
   testCase.verifyFalse(isfield(results, 'physics_fingerprint') ...
      && ~isempty(results.physics_fingerprint));

   clear cleaner restore
end

function test_reduced_model_is_not_completed(testCase)
   % A model missing a required channel must remain unavailable.
   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      artifact_root=testCase.TestData.artifact_root, ...
      run_name="20260804-010209", ...
      write_artifacts=true, model_provider=@reducedModelProvider);

   kanm = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(kanm), "unavailable")
   testCase.verifyEqual(results.summary.error_identifier(kanm), ...
      "icemodel:verification:compareAblation:missingModelField")

   artifact = load(results.paths.results_mat, 'results');
   testCase.verifyEqual( ...
      string(artifact.results.site_year_results.status), "unavailable")
end

function test_future_display_interval_sentinel_is_rejected(testCase)
   % A provider row stamped at October 1 represents the forbidden interval.
   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      artifact_root=testCase.TestData.artifact_root, ...
      run_name="20260804-010204", write_artifacts=true, ...
      model_provider=@futureIntervalProvider);
   kanm = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(kanm), "unavailable")
   testCase.verifyTrue(isnan( ...
      results.summary.window_valid_sample_count(kanm)))
   testCase.verifyTrue(isnan(results.summary.coverage_fraction(kanm)))
   testCase.verifyTrue(isnan(results.summary.eligible_sample_count(kanm)))
   testCase.verifyEqual(results.summary.error_identifier(kanm), ...
      "icemodel:verification:promiceAblationEvaluation:futureInterval")
   testCase.verifyEqual(height(results.nested_windows), 4)
   testCase.verifyFalse(any(results.nested_windows.available))
   testCase.verifyEqual(height(results.endpoint_perturbations), 6)
   testCase.verifyFalse(any(results.endpoint_perturbations.available))
   artifact = load(results.paths.results_mat, 'results');
   testCase.verifyEqual(numel(artifact.results.site_year_results), 1)
   testCase.verifyEqual( ...
      string(artifact.results.site_year_results.status), "unavailable")
   testCase.verifyEqual(artifact.results.summary, results.summary)
end

function test_internal_display_gap_is_rejected(testCase)
   % A synthetic endpoint must not hide a missing hourly row inside the season.
   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      artifact_root=testCase.TestData.artifact_root, ...
      run_name="20260804-010205", write_artifacts=true, ...
      model_provider=@internalGapProvider);
   kanm = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(kanm), "unavailable")
   testCase.verifyEqual(results.summary.error_identifier(kanm), ...
      "icemodel:verification:promiceAblationEvaluation:noncontiguousDisplaySupport")
end

function test_missing_observation_posting_remains_a_visible_gap(testCase)
   % A missing source posting stays on the seasonal grid as censored NaN.
   fixture = isolatedFixture();
   cleanup = onCleanup(@() removeOwnedTree(fixture.root));
   observation_file = fullfile(fixture.eval_root, ...
      'promice', 'kanm', 'observations.mat');
   saved = load(observation_file, 'targets');
   gap_time = datetime(2019, 6, 10, 12, 0, 0, 'TimeZone', 'UTC');
   saved.targets.data(saved.targets.data.Time == gap_time, :) = [];
   targets = saved.targets;
   save(observation_file, 'targets')

   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=fixture.eval_root, ...
      input_data_root=fixture.input_root, ...
      artifact_root=fullfile(fixture.root, 'artifacts'), ...
      run_name="20260804-010206", write_artifacts=true, ...
      model_provider=@selectedModelProvider);
   selected = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(selected), "completed", ...
      results.summary.reason(selected))
   seasonal = results.site_year_results(selected).seasonal;
   gap = seasonal.Time == gap_time;
   testCase.verifyEqual(nnz(gap), 1)
   testCase.verifyFalse(seasonal.direct_observation(gap))
   testCase.verifyTrue(isnan(seasonal.observation_lowering_m(gap)))
   testCase.verifyTrue(isnan(seasonal.model_melt_mwe(gap)))
   testCase.verifyTrue(isfinite(seasonal.model_melt_mwe( ...
      seasonal.Time == gap_time + hours(1))))
end

function test_negative_snow_is_unknown_in_seasonal_diagnostics(testCase)
   % A negative source depth remains visible but never marks exposed ice.
   fixture = isolatedFixture();
   cleanup = onCleanup(@() removeOwnedTree(fixture.root));
   observation_file = fullfile(fixture.eval_root, ...
      'promice', 'kanm', 'observations.mat');
   saved = load(observation_file, 'targets');
   negative_time = datetime(2019, 6, 10, 13, 0, 0, ...
      'TimeZone', 'UTC');
   negative = saved.targets.data.Time == negative_time;
   saved.targets.data.snow_depth(negative) = -0.005;
   targets = saved.targets;
   save(observation_file, 'targets')

   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=fixture.eval_root, ...
      input_data_root=fixture.input_root, ...
      artifact_root=fullfile(fixture.root, 'artifacts'), ...
      run_name="20260804-010207", write_artifacts=true, ...
      model_provider=@selectedModelProvider);
   selected = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(selected), "completed", ...
      results.summary.reason(selected))
   seasonal = results.site_year_results(selected).seasonal;
   row = seasonal.Time == negative_time;
   testCase.verifyEqual(nnz(row), 1)
   testCase.verifyEqual(seasonal.snow_depth_m(row), -0.005)
   testCase.verifyTrue(seasonal.direct_observation(row))
   testCase.verifyFalse(seasonal.ice_exposed(row))
   testCase.verifyTrue(isnan(seasonal.observation_lowering_m(row)))
   testCase.verifyTrue(isnan(seasonal.model_melt_mwe(row)))
end

function test_negative_lowering_keeps_seasonal_density_band_ordered(testCase)
   % A direct exposed-ice dip below the common reference must swap density
   % endpoints while keeping the lower and upper fields numeric.
   fixture = isolatedFixture();
   cleanup = onCleanup(@() removeOwnedTree(fixture.root));
   observation_file = fullfile(fixture.eval_root, ...
      'promice', 'kanm', 'observations.mat');
   saved = load(observation_file, 'targets');
   negative_time = datetime(2019, 6, 10, 12, 0, 0, ...
      'TimeZone', 'UTC');
   saved.targets.data.ablation = -linspace( ...
      0, 1, height(saved.targets.data))';
   targets = saved.targets;
   save(observation_file, 'targets')

   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=fixture.eval_root, ...
      input_data_root=fixture.input_root, ...
      artifact_root=fullfile(fixture.root, 'artifacts'), ...
      run_name="20260804-010208", write_artifacts=true, ...
      model_provider=@selectedModelProvider);
   selected = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(selected), "completed", ...
      results.summary.reason(selected))
   seasonal = results.site_year_results(selected).seasonal;
   row = seasonal.Time == negative_time;
   [~, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   % Derive the band from the policy rather than restating its endpoints, so
   % this stays correct when the density decision changes.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   density_values = seasonal.observation_lowering_m(row) ...
      .* policy.effective_density_kg_m3 ./ ro_liq;
   testCase.verifyLessThan(seasonal.observation_lowering_m(row), 0)
   testCase.verifyEqual(seasonal.observation_lower_mwe(row), ...
      min(density_values), AbsTol=1e-14)
   testCase.verifyEqual(seasonal.observation_upper_mwe(row), ...
      max(density_values), AbsTol=1e-14)
   finite_band = isfinite(seasonal.observation_lower_mwe) ...
      & isfinite(seasonal.observation_upper_mwe);
   testCase.verifyTrue(all(seasonal.observation_lower_mwe(finite_band) ...
      <= seasonal.observation_upper_mwe(finite_band)))

   % The machine-readable endpoint ledger must retain the same numeric ordering
   % even though its density-specific endpoint sensitivities are negative.
   perturbations = readtable(results.paths.endpoint_perturbations_csv);
   available = logical(perturbations.available);
   testCase.verifyTrue(any( ...
      perturbations.observation_sensitivity_max_mwe(available) < 0))
   testCase.verifyTrue(all( ...
      perturbations.observation_sensitivity_min_mwe(available) ...
      <= perturbations.observation_sensitivity_max_mwe(available)))
end

function test_provider_failure_retains_predeclared_nested_rows(testCase)
   % A selected execution failure must remain visible at every nested duration.
   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      model_provider=@failingProvider);
   kanm = results.summary.case_id == "kanm";
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   nested = results.nested_windows;

   % Labels come from the central policy, and the same stable failure is copied
   % to the flat, per-site-year, and nested machine-readable artifacts.
   expected_labels = [compose('%dd', policy.nested_window_days(:)); "longest"];
   testCase.verifyEqual(results.summary.status(kanm), "unavailable")
   testCase.verifyEqual(nested.window_label, expected_labels)
   testCase.verifyFalse(any(nested.available))
   testCase.verifyTrue(all(nested.reason == "synthetic provider failure"))
   testCase.verifyTrue(all(nested.error_identifier == "test:providerFailure"))
   site_result = results.site_year_results(kanm);
   testCase.verifyEqual(site_result.nested_windows, nested)
   testCase.verifyEqual(site_result.error_identifier, "test:providerFailure")
   testCase.verifyTrue(all(~isnat(nested.window_start)))
   testCase.verifyTrue(all(~isnat(nested.window_end)))

   % Endpoint rows are also predeclared, so a model failure cannot erase the
   % requested sensitivity plan or its stable failure identity.
   perturbations = results.endpoint_perturbations;
   testCase.verifyEqual(height(perturbations), ...
      2 * numel(policy.endpoint_perturbation_days))
   testCase.verifyFalse(any(perturbations.available))
   testCase.verifyTrue(all( ...
      perturbations.reason == "synthetic provider failure"))
   testCase.verifyTrue(all( ...
      perturbations.error_identifier == "test:providerFailure"))
   testCase.verifyTrue(all(~isnat(perturbations.requested_window_start)))
   testCase.verifyTrue(all(~isnat(perturbations.requested_window_end)))
   testCase.verifyTrue(all(isnat(perturbations.actual_window_start)))
   testCase.verifyEqual(site_result.endpoint_perturbations, perturbations)
end

function test_invalid_selection_and_provider_fail_before_execution(testCase)
   % Typos and non-callable seams must fail as runner configuration errors.
   base = { ...
      'evaluation_data_root', testCase.TestData.eval_root, ...
      'input_data_root', testCase.TestData.input_root};
   testCase.verifyError(@() run_promice_ablation_evaluation( ...
      base{:}, case_ids="missing"), ...
      'icemodel:verification:promiceAblationEvaluation:unknownCase')
   testCase.verifyError(@() run_promice_ablation_evaluation( ...
      base{:}, case_ids="kanm", years=2020), ...
      'icemodel:verification:promiceAblationEvaluation:emptySelection')
   testCase.verifyError(@() run_promice_ablation_evaluation( ...
      base{:}, model_provider=42), ...
      'icemodel:verification:promiceAblationEvaluation:modelProvider')

   % Readiness-only runs are legitimate, but asking them to write artifacts
   % would leave a populated run directory and a renderable report describing
   % zero site-years, misrepresenting the run as a completed cohort evaluation.
   testCase.verifyError(@() run_promice_ablation_evaluation( ...
      base{:}, write_artifacts=true, model_provider=@unexpectedProvider), ...
      'icemodel:verification:promiceAblationEvaluation:emptySelection')
end

function test_tampered_selected_root_readiness_blocks_execution(testCase)
   % Changed producer-ledger bytes must fail before any model callback executes.
   fixture = isolatedFixture();
   cleanup = onCleanup(@() removeOwnedTree(fixture.root));
   readiness_file = fullfile(fixture.root, 'preview', 'qa', 'gapfill', ...
      'ledger', 'kanm-readiness.csv');

   % Keep the producer manifest pinned to the original digest, then mutate its
   % selected-root readiness artifact exactly as a stale transaction would.
   fid = fopen(readiness_file, 'a');
   testCase.assertGreaterThanOrEqual(fid, 0)
   file_cleanup = onCleanup(@() fclose(fid));
   fwrite(fid, uint8('tampered'), 'uint8');
   clear file_cleanup

   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=fixture.eval_root, ...
      input_data_root=fixture.input_root, ...
      model_provider=@unexpectedProvider);
   kanm = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(kanm), "excluded")
   testCase.verifySubstring(results.summary.reason(kanm), ...
      "producer-pinned artifact identity is invalid")
end

function test_escaping_selected_root_report_inputs_blocks_execution(testCase)
   % A producer-manifest symlink outside the selected root must never execute.
   fixture = isolatedFixture();
   cleanup = onCleanup(@() removeOwnedTree(fixture.root));
   outside_root = string(tempname);
   mkdir(outside_root)
   outside_cleanup = onCleanup(@() removeOwnedTree(outside_root));
   report_inputs_file = fullfile(fixture.root, 'preview', 'qa', 'gapfill', ...
      'plans', 'kanm-report-inputs.json');
   outside_file = fullfile(outside_root, 'kanm-report-inputs.json');

   % Readiness can parse the aliased producer bytes, but execution-time root
   % validation must resolve the symlink and reject the escaped verifier path.
   movefile(report_inputs_file, outside_file)
   [status, message] = system(sprintf('ln -s %s %s', ...
      char(outside_file), char(report_inputs_file)));
   testCase.assertEqual(status, 0, message)
   results = run_promice_ablation_evaluation( ...
      case_ids="kanm", years=2019, ...
      evaluation_data_root=fixture.eval_root, ...
      input_data_root=fixture.input_root, ...
      model_provider=@unexpectedProvider);
   kanm = results.summary.case_id == "kanm";
   testCase.verifyEqual(results.summary.status(kanm), "unavailable")
   testCase.verifyEqual(results.summary.error_identifier(kanm), ...
      "icemodel:verification:artifactIdentity:relativePath")
end

function test_pinned_hash_rejects_changed_bytes(testCase)
   % Execution-time identity checks must not trust an earlier readiness hash.
   pathname = fullfile(testCase.TestData.root, 'mutable-artifact.bin');
   fid = fopen(pathname, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, uint8('original'), 'uint8');
   clear cleaner
   expected = icemodel.verification.setup.fileSha256(pathname);
   icemodel.verification.helpers.assertArtifactSha256(pathname, expected)

   fid = fopen(pathname, 'a');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, uint8('-changed'), 'uint8');
   clear cleaner
   testCase.verifyError(@() ...
      icemodel.verification.helpers.assertArtifactSha256( ...
      pathname, expected), ...
      'icemodel:verification:artifactIdentity:mismatch')
   testCase.verifyError(@() ...
      icemodel.verification.helpers.assertArtifactSha256( ...
      pathname, "not-a-sha256"), ...
      'icemodel:verification:artifactIdentity:sha256')
   testCase.verifyError(@() ...
      icemodel.verification.helpers.assertArtifactSha256( ...
      pathname + ".missing", expected), ...
      'icemodel:verification:artifactIdentity:missing')

   % Root-relative identity checks cover the producer manifest and readiness
   % ledger independently, including changed bytes and unsafe path spellings.
   provenance_root = fullfile(testCase.TestData.root, 'provenance-helper');
   mkdir(provenance_root)
   names = ["producer-manifest.json", "producer-readiness.csv"];
   for k = 1:numel(names)
      relative_path = names(k);
      artifact_path = fullfile(provenance_root, relative_path);
      fid = fopen(artifact_path, 'w');
      cleaner = onCleanup(@() fclose(fid));
      fwrite(fid, uint8(char("original-" + relative_path)), 'uint8');
      clear cleaner
      expected = icemodel.verification.setup.fileSha256(artifact_path);
      resolved = ...
         icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
         provenance_root, relative_path, expected);
      testCase.verifyEqual(resolved, string(artifact_path))

      % Each producer link must fail under its originally pinned identity after
      % mutation; one passing sibling cannot mask another changed artifact.
      fid = fopen(artifact_path, 'a');
      cleaner = onCleanup(@() fclose(fid));
      fwrite(fid, uint8('-changed'), 'uint8');
      clear cleaner
      testCase.verifyError(@() ...
         icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
         provenance_root, relative_path, expected), ...
         'icemodel:verification:artifactIdentity:mismatch')
   end
   testCase.verifyError(@() ...
      icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
      provenance_root, "../mutable-artifact.bin", expected), ...
      'icemodel:verification:artifactIdentity:relativePath')
   testCase.verifyError(@() ...
      icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
      provenance_root, pathname, expected), ...
      'icemodel:verification:artifactIdentity:relativePath')
end

function model = unexpectedProvider(varargin)
   %UNEXPECTEDPROVIDER Fail if the readiness-only default invokes a model.
   model = timetable();
   assert(~isempty(model), 'test:unexpectedModelExecution', ...
      'readiness-only execution called the model provider')
end

function model = failingProvider(varargin)
   %FAILINGPROVIDER Raise one stable selected-row execution failure.
   model = timetable();
   assert(~isempty(model), 'test:providerFailure', ...
      'synthetic provider failure')
end

function model = selectedModelProvider(manifest, row, run_start, run_end)
   %SELECTEDMODELPROVIDER Return a closed hourly half-open diagnostic ledger.
   evaluation_start = parseTime(row.snow_free_window_start);
   initialization_start = parseTime(row.requested_window_start);
   display_start = datetime(double(row.year), 6, 1, 'TimeZone', 'UTC');
   display_end = datetime(double(row.year), 10, 1, 'TimeZone', 'UTC');
   % This fixture has five snow-covered days before its observed snow-free
   % interval; the canonical manifest period independently owns initialization.
   observation_start = evaluation_start - days(5);
   assert(run_end == display_end - minutes(15))
   assert(run_start == initialization_start)
   assert(run_start < evaluation_start)
   assert(string(manifest.forcing_sources) == "promice_filled")
   assert(isfield(manifest.colocation, 'promice_filled'))
   assert(startsWith(string( ...
      manifest.colocation.promice_filled.met_files), "promice_filled/"))

   % Resolve production options without running the solver. This proves the
   % overlay reaches opts.forcings and every exact producer-pinned input path.
   opts = icemodel.test.helpers.setModelOptsForCase( ...
      manifest, startdate=run_start, enddate=run_end, ...
      output_profile="diagnostic");
   assert(string(opts.forcings) == "promice_filled")
   assert(isscalar(opts.metfname))
   assert(contains(string(opts.metfname{1}), ...
      fullfile('met', 'promice_filled')))
   assert(string(opts.readiness_file) == string(manifest.readiness_file))
   assert(string(opts.report_inputs_file) ...
      == string(manifest.report_inputs_file))
   selected_root = string(fileparts(manifest.input_data_root));
   assert(icemodel.isPathInside(opts.readiness_file, selected_root))
   assert(icemodel.isPathInside( ...
      opts.report_inputs_file, selected_root))

   % Exercise the standard loadmet verifier against this alternate selected
   % root. Repository-local readiness defaults cannot satisfy these identities.
   [~, verified_opts] = icemodel.loadmet(opts);
   assert(verified_opts.promice_filled_readiness_verified)
   assert(verified_opts.promice_filled_manifest_verified)
   assert(string(verified_opts.report_inputs_file) ...
      == string(manifest.report_inputs_file))

   % Match the linear synthetic lowering after intact-ice conversion while
   % keeping every storage, phase, and remesh identity exactly closed.
   time = (run_start:hours(1):display_end - hours(1)).';
   n = numel(time);
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   model = timetable('RowTimes', time);
   fields = icemodel.namelists.budgetoutputs('all');
   for k = 1:numel(fields)
      model.(fields{k}) = zeros(n, 1);
   end
   solid_loss = zeros(n, 1);
   observation_rows = time >= observation_start;
   solid_loss(observation_rows) = (ro_ice / ro_liq) ...
      / hours(display_end - display_start);
   model.mass_budget_solid_start_mwe = ...
      10 - [0; cumsum(solid_loss(1:end - 1))];
   model.mass_budget_solid_end_mwe = ...
      model.mass_budget_solid_start_mwe - solid_loss;
   model.mass_budget_liquid_start_mwe(:) = 1;
   model.mass_budget_liquid_end_mwe(:) = 1;
   model.mass_budget_phase_solid_mwe = -solid_loss;
   model.mass_budget_phase_liquid_mwe = solid_loss;
   % remesh_liquid_mwe alone closes the liquid_storage checkpoint; the
   % schema carries no B/O liquid split (cloned_bottom_liquid,
   % merge_export_liquid).
   model.mass_budget_remesh_liquid_mwe = -solid_loss;
   deletion = time == evaluation_start;
   model.mass_budget_top_deletion_count(deletion) = 1;
   model.mass_budget_top_deletion_height_m(deletion) = 0.01;
   % Exported mass is a separate channel from the quantized cell height, so the
   % fixture gives it a distinct value that the surface-loss series must carry.
   model.mass_budget_top_export_solid_mwe(deletion) = 0.006;
   model.mass_budget_top_export_liquid_mwe(deletion) = 0.004;
   model.melt = cumsum(solid_loss);
   model.runoff = 0.8 .* model.melt;
   model.freeze = 0.2 .* model.melt;
   model.dlayer = cumsum(solid_loss) .* ro_liq ./ ro_ice;
end

function model = reducedModelProvider(manifest, row, run_start, run_end)
   %REDUCEDMODELPROVIDER Remove one required grid diagnostic.
   model = selectedModelProvider(manifest, row, run_start, run_end);
   model = removevars(model, "mass_budget_interior_merge_count");
end

function model = futureIntervalProvider(manifest, row, run_start, run_end)
   %FUTUREINTERVALPROVIDER Add a sentinel at the forbidden October 1 boundary.
   model = selectedModelProvider(manifest, row, run_start, run_end);
   future = model(end, :);
   future.Properties.RowTimes = datetime( ...
      double(row.year), 10, 1, 'TimeZone', 'UTC');
   future.mass_budget_phase_solid_mwe = -100;
   model = [model; future];
end

function model = internalGapProvider(manifest, row, run_start, run_end)
   %INTERNALGAPPROVIDER Remove one hourly interval inside the fixed season.
   model = selectedModelProvider(manifest, row, run_start, run_end);
   gap_time = datetime(double(row.year), 7, 15, 'TimeZone', 'UTC');
   model(model.Time == gap_time, :) = [];
end

function writeFixtureTree(eval_root, input_root)
   %WRITEFIXTURETREE Create admitted KAN_M and excluded KAN_L annual rows.
   family_root = fullfile(eval_root, 'promice');
   mkdir(family_root)
   start_time = datetime(2019, 6, 1, 0, 0, 0, 'TimeZone', 'UTC');
   end_time = datetime(2019, 10, 1, 0, 0, 0, 'TimeZone', 'UTC');
   forcing_start = datetime(2019, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   forcing_end = datetime(2019, 12, 31, 23, 0, 0, 'TimeZone', 'UTC');
   cases = repmat(caseTemplate(), 2, 1);
   cases(1) = writeObservationCase( ...
      family_root, "kanm", "KAN_M", start_time, end_time);
   cases(2) = writeObservationCase( ...
      family_root, "kanl", "KAN_L", start_time, end_time);
   for k = 1:numel(cases)
      cases(k).period = periodStruct(forcing_start, forcing_end);
   end
   manifest = struct( ...
      'dataset_family', 'promice', 'source_doi', '', ...
      'source_url', 'https://promice.org', ...
      'source_version', 'synthetic', 'retrieval_date', '2026-08-04', ...
      'cases', cases, 'skipped', struct([]));
   icemodel.verification.setup.writeJson( ...
      fullfile(family_root, 'manifest.json'), manifest)

   % Only KAN_M receives producer-pinned forcing; KAN_L exercises retained
   % readiness exclusion without introducing another model execution.
   writeForcing(input_root, "kanm", forcing_start, forcing_end)
end

function fixture = isolatedFixture()
   %ISOLATEDFIXTURE Create one independently mutable selected-root fixture.
   root = string(tempname);
   eval_root = fullfile(root, 'eval');
   input_root = fullfile(root, 'input');

   % Keep every generated producer artifact under the test-owned root so
   % tamper and escape cases cannot alter the shared setupOnce fixture.
   mkdir(eval_root)
   mkdir(input_root)
   writeFixtureTree(eval_root, input_root)
   fixture = struct( ...
      'root', root, 'eval_root', string(eval_root), ...
      'input_root', string(input_root));
end

function removeOwnedTree(root)
   %REMOVEOWNEDTREE Remove one test-created directory when it still exists.
   if isfolder(root)
      rmdir(root, 's')
   end
end

function c = writeObservationCase( ...
      family_root, case_id, site_id, start_time, end_time)
   %WRITEOBSERVATIONCASE Save one source-faithful hourly lowering target.
   time = (start_time:hours(1):end_time).';
   n = numel(time);
   data = timetable('RowTimes', time);
   data.ablation = linspace(0, 1, n).';
   data.snow_depth = repmat(0.005, n, 1);
   data.snow_depth(time < start_time + days(5)) = 0.2;
   data.snow_depth(time > start_time + days(40)) = 0.2;
   data.surface_height_flag = zeros(n, 1);
   data.station_transition_flag = zeros(n, 1);
   data.step_detected_flag = zeros(n, 1);
   data.step_correctable_flag = zeros(n, 1);
   data.unused_payload = transpose(1:n);
   data.Properties.VariableUnits = {'m', 'm', '1', '1', '1', '1', '1'};
   data.Properties.VariableDescriptions = { ...
      'surface ablation (lowering) height', 'snow depth', 'gap flag', ...
      'station transition flag', 'step detected flag', ...
      'step correctable flag', 'unused payload sentinel'};
   targets = struct( ...
      'format', 'timeseries', 'data', data, ...
      'metadata', struct('source', 'synthetic', ...
      'source_family', 'promice', 'station', char(site_id), ...
      'site_id', char(site_id)));
   case_root = fullfile(family_root, case_id);
   mkdir(case_root)
   save(fullfile(case_root, 'observations.mat'), 'targets')

   % The runner replaces this raw manifest forcing leg with readiness identity.
   c = caseTemplate();
   c.case_id = char(case_id);
   c.case_type = 'firn_observational';
   c.site_id = char(site_id);
   c.surface_zone = 'ablation';
   c.period = periodStruct(start_time, end_time);
   c.evaluation_file = char(fullfile(case_id, 'observations.mat'));
   c.forcing_sources = {'promice'};
end

function writeForcing(input_root, case_id, start_time, end_time)
   %WRITEFORCING Save one filled met file and its producer evidence.
   data_root = string(fileparts(input_root));
   met_dir = fullfile(input_root, 'met', 'promice_filled');
   ledger_dir = fullfile(data_root, 'preview', 'qa', 'gapfill', 'ledger');
   plans_dir = fullfile(data_root, 'preview', 'qa', 'gapfill', 'plans');
   mkdir(met_dir)
   mkdir(ledger_dir)
   mkdir(plans_dir)

   % The real 15-minute met file lets option resolution prove the
   % cadence and exact met path without invoking the production solver.
   time = (start_time:minutes(15):end_time).';
   n = numel(time);
   met = timetable('RowTimes', time);
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
   met.swu = met.albedo .* met.swd;
   met.boom_height = repmat(2.7, n, 1);

   % Stamp the same product hash and per-channel provenance that
   % production loadmet validates after the root-scoped producer files pass.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   defaults = icemodel.forcing.reconstruct.setopts();
   provenance_channels = unique([defaults.plan_channels, ...
      icemodel.forcing.helpers.precipitationVariables(), ...
      "boom_height"], 'stable');
   for channel = provenance_channels
      provenance = repmat(codes.observed, n, 1);
      provenance(~isfinite(met.(channel))) = codes.missing;
      met.(channel + "_provenance") = provenance;
   end
   met.Properties.UserData = struct( ...
      'site', string(case_id), ...
      'gapfill_registry', codes, ...
      'gapfill_seed', 1, ...
      'gapfill_product', "promice_filled", ...
      'gapfill_engine_version', string(icemodel.internal.version()), ...
      'gapfill_policy_sha256', ...
      icemodel.forcing.reconstruct.policySha256(), ...
      'gapfill_donors', string.empty(1, 0), ...
      'gapfill_channels', defaults.plan_channels);
   artifact_metadata = icemodel.forcing.helpers.artifactMetadata(met);
   met.Properties.UserData = artifact_metadata;
   met_name = "met_" + case_id ...
      + "_promice_filled_" + string(start_time, 'yyyyMMdd') + "_" ...
      + string(end_time, 'yyyyMMdd') + "_15m.mat";
   met_file = fullfile(met_dir, met_name);
   save(met_file, 'met', 'artifact_metadata')

   % Consumer verdicts and hashes are the only forcing admission source.
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
   producer = struct( ...
      'site', char(case_id), 'path_base', 'selected_data_root', ...
      'artifacts', artifacts, ...
      'acceptance_window', struct( ...
      'start', char(formatTime(start_time)), ...
      'end', char(formatTime(end_time))));
   icemodel.verification.setup.writeJson( ...
      fullfile(plans_dir, case_id + "-report-inputs.json"), producer)
end

function artifact = artifactRecord(role, path, pathname)
   %ARTIFACTRECORD Build one producer-pinned file identity.
   info = dir(pathname);
   artifact = struct( ...
      'role', char(role), 'path', char(path), 'bytes', info.bytes, ...
      'sha256', char( ...
      icemodel.verification.setup.fileSha256(pathname)));
end

function c = caseTemplate()
   %CASETEMPLATE Keep the two synthetic manifest entries homogeneous.
   c = struct( ...
      'case_id', '', 'case_type', '', 'site_id', '', 'surface_zone', '', ...
      'period', struct('start', '', 'end', ''), 'evaluation_file', '', ...
      'forcing_sources', {{}}, 'colocation', struct());
end

function period = periodStruct(first, last)
   %PERIODSTRUCT Encode one UTC manifest window.
   period = struct( ...
      'start', char(formatTime(first)), 'end', char(formatTime(last)));
end

function time = parseTime(value)
   %PARSETIME Convert one readiness timestamp to UTC.
   time = icemodel.verification.setup.ensureUtc(string(value));
end

function text = formatTime(value)
   %FORMATTIME Format one UTC timestamp for portable JSON.
   value.TimeZone = 'UTC';
   text = string(value, 'yyyy-MM-dd HH:mm:ss');
end
