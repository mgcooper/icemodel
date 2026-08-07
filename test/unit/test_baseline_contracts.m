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
      2, ...
      legacy_time, ...
      'VariableNames', {'case_id', 'runoff_final', 'solver_mode', ...
      'last_updated_utc'});
   save(filepath, 'baseline');

   loaded = icemodel.test.helpers.loadBaseline("regression", ...
      smbmodel="icemodel", filename=filepath);

   testCase.verifyEqual(loaded.case_id, "icemodel_kanm_2016_solver2");
   testCase.verifyEqual(loaded.baseline_tag, "");
   testCase.verifyEqual(loaded.baseline_type, "rolling");
   testCase.verifyEqual(loaded.smbmodel_filter, "icemodel");
   testCase.verifyEqual(loaded.solver_mode, 2);
   testCase.verifyEqual(loaded.solver, 2);
   testCase.verifyTrue(isdatetime(loaded.last_updated_utc));
   clear cleanup
end

function test_frozen_v11_baselines_normalize_without_mutation(testCase)
   % All four immutable release files must acquire the canonical solver alias
   % in memory and pass their registered forcing/case identity checks.
   for kind = ["regression", "perf"]
      for model = ["icemodel", "skinmodel"]
         baseline = icemodel.test.helpers.loadBaseline(kind, ...
            smbmodel=model, baseline_tag="v1.1", simyear=2016);
         testCase.verifyNotEmpty(baseline);
         testCase.verifyTrue(ismember( ...
            'solver_mode', baseline.Properties.VariableNames));
         testCase.verifyEqual(baseline.solver, ...
            double(baseline.solver_mode));
         icemodel.test.helpers.assertFormalBaselineForcing( ...
            baseline, "v1.1");
      end
   end
end

function test_loadBaseline_rejects_solver_and_selector_conflicts(testCase)
   % Existing persisted identity must never be silently replaced by requested
   % selector metadata or a canonical alias derived from another column.
   filepath = [tempname '.mat'];
   cleanup = onCleanup(@() deleteIfExists(filepath));
   RegressionBaseline = table( ...
      "icemodel_kanm_2016_solver1", 1, 2, ...
      'VariableNames', {'case_id', 'solver', 'solver_mode'});
   save(filepath, 'RegressionBaseline');
   testCase.verifyError(@() icemodel.test.helpers.loadBaseline( ...
      "regression", filename=filepath), ...
      'icemodel:test:baselineSolverSchemaMismatch');

   RegressionBaseline = table( ...
      "icemodel_kanm_2016_solver1", "release", "v1.1", ...
      'VariableNames', {'case_id', 'baseline_type', 'baseline_tag'});
   save(filepath, 'RegressionBaseline');
   testCase.verifyError(@() icemodel.test.helpers.loadBaseline( ...
      "regression", filename=filepath), ...
      'icemodel:test:baselineSelectorIdentityMismatch');

   RegressionBaseline.baseline_type = "rolling";
   RegressionBaseline.baseline_tag = "";
   meta = struct('baseline_type', "release", 'baseline_tag', "v1.1");
   save(filepath, 'RegressionBaseline', 'meta');
   testCase.verifyError(@() icemodel.test.helpers.loadBaseline( ...
      "regression", filename=filepath), ...
      'icemodel:test:baselineSelectorIdentityMismatch');
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

function test_getRegressionCaseMatrix_uses_official_forcing(testCase)
   % Default/rolling regression rows must select the accepted PROMICE product.

   cases = icemodel.test.helpers.getRegressionCaseMatrix(tier="full");

   testCase.verifyEqual(unique(cases.forcings), "promice_filled");
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

function test_getPerfCaseMatrix_uses_official_forcing(testCase)
   % Rolling performance acceptance must use the same forcing as regression.

   cases = icemodel.test.helpers.getPerfCaseMatrix(tier="full");

   testCase.verifyEqual(unique(cases.forcings), "promice_filled");
end

function test_getFormalForcing_returns_official_product(testCase)
   % Keep the shared formal forcing identity explicit and independently tested.

   testCase.verifyEqual( ...
      icemodel.test.helpers.getFormalForcing(), "promice_filled");
end

function test_formal_baseline_policy_owns_default_data_case(testCase)
   % The same registration selects forcing identity and its default data tree.

   rolling = icemodel.test.helpers.formalBaselinePolicy("rolling");
   release = icemodel.test.helpers.formalBaselinePolicy("v1.1");

   testCase.verifyEqual(rolling.config_case, "verification");
   testCase.verifyEqual(rolling.forcing, "promice_filled");
   testCase.verifyEmpty(rolling.required_fixture_capabilities);
   testCase.verifyEqual(release.config_case, "test");
   testCase.verifyEqual(release.site_forcings, ["kanm"; "kanl"]);
   testCase.verifyEqual(release.required_fixture_capabilities, "formal-core");
   testCase.verifyFalse(release.snapshot_from_rolling);
end

function test_release_v11_case_matrices_retain_historical_forcing(testCase)
   % Frozen v1.1 comparisons must recreate each station-specific forcing
   % identity rather than pairing promice_filled runs with legacy rows that
   % happen to retain the same case ids.

   regression_cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
      tier="full", baseline="v1.1");
   perf_cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier="full", baseline="v1_1");

   testCase.verifyEqual( ...
      regression_cases.forcings, lower(regression_cases.sitename));
   testCase.verifyEqual(perf_cases.forcings, lower(perf_cases.sitename));
end

function test_release_v11_matrices_match_frozen_baseline_forcing(testCase)
   % Check the registered identity against the preserved release files, not
   % merely against a second declaration in the test.

   regression_cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
      tier="full", baseline="v1.1");
   perf_cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier="full", baseline="v1.1");
   regression_baseline = icemodel.test.helpers.loadBaseline( ...
      "regression", smbmodel="all", baseline_tag="v1.1");
   perf_baseline = icemodel.test.helpers.loadBaseline( ...
      "perf", smbmodel="all", baseline_tag="v1.1", simyear=2016);

   verifyCaseForcingMatches(testCase, regression_cases, regression_baseline);
   verifyCaseForcingMatches(testCase, perf_cases, perf_baseline);
end

function test_unregistered_release_forcing_is_rejected(testCase)
   % New release tags require an explicit forcing registration; they must not
   % silently inherit the rolling product or the v1.1 station aliases.

   testCase.verifyError(@() ...
      icemodel.test.helpers.getRegressionCaseMatrix( ...
      tier="smoke", baseline="v9.9"), ...
      'icemodel:test:unregisteredReleaseForcing');
   testCase.verifyError(@() ...
      icemodel.test.helpers.getPerfCaseMatrix( ...
      tier="smoke", baseline="v9.9"), ...
      'icemodel:test:unregisteredReleaseForcing');
end

function test_v11_release_rejects_unregistered_site_forcing(testCase)
   % v1.1 has only the two historical formal station aliases.

   testCase.verifyError(@() ...
      icemodel.test.helpers.getFormalForcing( ...
      sitename="zaca", baseline="v1.1"), ...
      'icemodel:test:unsupportedReleaseForcingSite');
end

function test_formal_baseline_forcing_identity_is_checked(testCase)
   % Case-id matching must not hide a forcing mismatch or absent provenance.

   compatible = table( ...
      ["icemodel_kanm_2016_solver1"; "icemodel_kanl_2016_solver1"], ...
      ["icemodel"; "icemodel"], ["kanm"; "kanl"], [2016; 2016], ...
      [1; 1], ["kanm"; "kanl"], ...
      'VariableNames', {'case_id', 'smbmodel', 'sitename', 'simyear', ...
      'solver', 'forcings'});
   icemodel.test.helpers.assertFormalBaselineForcing(compatible, "v1.1");

   incompatible = compatible;
   incompatible.forcings(1) = "promice_filled";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing( ...
      incompatible, "v1.1"), ...
      'icemodel:test:baselineForcingIdentityMismatch');

   missing = removevars(compatible, 'forcings');
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing(missing, "v1.1"), ...
      'icemodel:test:baselineForcingIdentityMissing');
end

function test_obsolete_rolling_forcing_requires_acceptance(testCase)
   % The pre-promice_filled rolling state must stop explicitly rather than
   % compare rows that share case ids but not forcing identity.

   legacy = table("icemodel_kanm_2016_solver1", "icemodel", ...
      "kanm", 2016, 1, "kanm", ...
      'VariableNames', {'case_id', 'smbmodel', 'sitename', 'simyear', ...
      'solver', 'forcings'});
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing(legacy, "rolling"), ...
      'icemodel:test:rollingBaselineForcingAcceptanceRequired');

   missing = removevars(legacy, 'forcings');
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing(missing, "rolling"), ...
      'icemodel:test:rollingBaselineForcingAcceptanceRequired');
end

function test_formal_baseline_case_identity_is_checked(testCase)
   % A case id must not allow one site's metrics to masquerade as another.

   baseline = table("icemodel_kanm_2016_solver1", "icemodel", ...
      "kanm", 2016, 1, "promice_filled", ...
      'VariableNames', {'case_id', 'smbmodel', 'sitename', 'simyear', ...
      'solver', 'forcings'});
   icemodel.test.helpers.assertFormalBaselineForcing(baseline, "rolling");

   mislabeled = baseline;
   mislabeled.sitename = "kanl";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing( ...
      mislabeled, "rolling"), ...
      'icemodel:test:baselineCaseIdentityMismatch');

   release_mislabeled = mislabeled;
   release_mislabeled.forcings = "kanl";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing( ...
      release_mislabeled, "v1.1"), ...
      'icemodel:test:baselineCaseIdentityMismatch');

   legacy = baseline;
   legacy.case_id = "full_icemodel_kanm_2016_bc1";
   icemodel.test.helpers.assertFormalBaselineForcing(legacy, "rolling");

   missing = removevars(baseline, 'solver');
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing(missing, "rolling"), ...
      'icemodel:test:baselineCaseIdentityMissing');

   duplicate = [baseline; baseline];
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineForcing( ...
      duplicate, "rolling"), ...
      'icemodel:test:baselineCaseIdentityDuplicate');
end

function test_regression_metric_evidence_is_selector_aware(testCase)
   % Rolling state is complete and fail-closed; frozen releases may omit only
   % metrics that did not exist in their persisted schema.
   baseline = table(1, 'VariableNames', {'runoff_final'});
   [compare, passed] = ...
      icemodel.test.helpers.formalRegressionMetricEvidence( ...
      1, baseline, 1, "runoff_final", "rolling");
   testCase.verifyTrue(compare);
   testCase.verifyTrue(passed);

   [compare, passed] = ...
      icemodel.test.helpers.formalRegressionMetricEvidence( ...
      1, baseline, 1, "closure_seb_rmse", "rolling");
   testCase.verifyFalse(compare);
   testCase.verifyFalse(passed);
   [compare, passed] = ...
      icemodel.test.helpers.formalRegressionMetricEvidence( ...
      1, baseline, 1, "closure_seb_rmse", "v1.1");
   testCase.verifyFalse(compare);
   testCase.verifyTrue(passed);

   baseline.runoff_final = NaN;
   [compare, passed] = ...
      icemodel.test.helpers.formalRegressionMetricEvidence( ...
      1, baseline, 1, "runoff_final", "rolling");
   testCase.verifyTrue(compare);
   testCase.verifyFalse(passed);
   baseline.runoff_final = 1;
   [~, passed] = ...
      icemodel.test.helpers.formalRegressionMetricEvidence( ...
      NaN, baseline, 1, "runoff_final", "v1.1");
   testCase.verifyFalse(passed);
end

function test_performance_verdict_fails_closed_when_compatible(testCase)
   % Compatibility authorizes comparison only; one usable row must still exist.
   baseline = table("case", 10, 0.2, ...
      'VariableNames', {'case_id', 'median_wall_s', 'tol_perf'});
   [passed, ref] = icemodel.test.helpers.formalPerformanceVerdict( ...
      true, 10, baseline, 1, true, 0.2, "");
   testCase.verifyTrue(passed);
   testCase.verifyEqual(ref, 10);

   [passed, ~, ~, ~, reason] = ...
      icemodel.test.helpers.formalPerformanceVerdict( ...
      true, 10, baseline, [], true, 0.2, "");
   testCase.verifyFalse(passed);
   testCase.verifySubstring(reason, "missing or duplicated");
   baseline.median_wall_s = NaN;
   [passed, ~, ~, ~, reason] = ...
      icemodel.test.helpers.formalPerformanceVerdict( ...
      true, 10, baseline, 1, true, 0.2, "");
   testCase.verifyFalse(passed);
   testCase.verifySubstring(reason, "finite and positive");
   [passed, ~, ~, ~, reason] = ...
      icemodel.test.helpers.formalPerformanceVerdict( ...
      true, 10, table(), [], false, 0.2, "metadata incompatible");
   testCase.verifyTrue(passed);
   testCase.verifyEqual(reason, "metadata incompatible");
end

function test_baseline_candidate_validation_precedes_publication(testCase)
   % Pure candidate checks cover exact membership and invalid numerical state
   % without launching an expensive builder or touching a managed baseline.
   cases = table("icemodel_kanm_2016_solver1", ...
      'VariableNames', {'case_id'});
   identity = {"icemodel_kanm_2016_solver1", "full", "rolling", "", ...
      "icemodel", "kanm", "promice_filled", 2016, 1};
   regression = table(identity{:}, 1, ...
      'VariableNames', {'case_id', 'tier', 'baseline_type', ...
      'baseline_tag', 'smbmodel', 'sitename', 'forcings', 'simyear', ...
      'solver', 'runoff_final'});
   icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "regression", regression, cases, "rolling");
   invalid_regression = regression;
   invalid_regression.runoff_final = NaN;
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "regression", invalid_regression, cases, "rolling"), ...
      'icemodel:test:baselineCandidateMetricsInvalid');
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "regression", [regression; regression], cases, "rolling"), ...
      'icemodel:test:baselineCandidateCasesInvalid');

   performance = regression(:, 1:9);
   performance.n_runs = 3;
   performance.n_warmups = 1;
   performance.tol_perf = 0.2;
   performance.median_wall_s = 10;
   performance.mean_wall_s = 10;
   performance.min_wall_s = 9;
   performance.max_wall_s = 11;
   performance.valid = true;
   performance.passed_perf = true;
   icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "perf", performance, cases, "rolling");
   performance.valid = false;
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "perf", performance, cases, "rolling"), ...
      'icemodel:test:baselineCandidateMetricsInvalid');

   benchmark = table(true, 3, 1.5, ...
      'VariableNames', {'Valid', 'SampleSize', 'Mean'});
   icemodel.test.helpers.assertFormalBenchmarkCandidate(benchmark);
   benchmark.Valid = false;
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBenchmarkCandidate(benchmark), ...
      'icemodel:test:benchmarkBaselineCandidateInvalid');
end

function test_builder_shaped_regression_candidate_excludes_only_metadata(testCase)
   % The actual builder emits string baseline identity and a datetime beside
   % numeric metrics. Those registered fields must not be treated as metrics,
   % while an unknown nonnumeric column must still fail the closed contract.
   case_id = "icemodel_kanm_2016_solver1";
   cases = table(case_id, 'VariableNames', {'case_id'});
   row = struct( ...
      'case_id', case_id, 'tier', "full", ...
      'baseline_type', "rolling", 'baseline_tag', "", ...
      'smbmodel', "icemodel", 'sitename', "kanm", ...
      'forcings', "promice_filled", 'simyear', 2016, 'solver', 1, ...
      'melt_final', 1.25, 'runoff_final', 0.75, ...
      'last_updated_utc', datetime('now', 'TimeZone', 'UTC'));
   candidate = struct2table(row);

   icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "regression", candidate, cases, "rolling");

   candidate.unregistered_note = "not a metric";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "regression", candidate, cases, "rolling"), ...
      'icemodel:test:baselineCandidateMetricsInvalid');
end

function test_existing_release_target_is_immutable(testCase)
   % The shared target guard rejects overwrite=true without changing bytes.

   output_file = string(tempname) + ".mat";
   cleanup = onCleanup(@() deleteIfExists(output_file));
   marker = "unchanged";
   save(output_file, 'marker');

   testCase.verifyError(@() ...
      icemodel.test.helpers.assertNewReleaseBaselineTarget( ...
      output_file, true), ...
      'icemodel:test:releaseBaselineImmutable');
   saved = load(output_file, 'marker');
   testCase.verifyEqual(saved.marker, marker);
   clear cleanup
end

function test_new_snapshot_path_requires_release_registration(testCase)
   % A nonexistent target is not rejected as immutable, but a new release tag
   % must register its forcing identity before any snapshot can be written.

   output_file = string(tempname) + ".mat";
   testCase.verifyError(@() ...
       icemodel.test.helpers.snapshotBaseline( ...
       "regression", "vNext", "icemodel", false, output_file, 2016), ...
      'icemodel:test:unregisteredReleaseForcing');
   testCase.verifyFalse(isfile(output_file));
end

function test_prepareBaselineBuild_forwards_release_identity(testCase)
   % Both acceptance builders use the same baseline-aware matrices as their
   % compare runners. This call resolves setup only and writes no baseline.

   output_file = string(tempname) + ".mat";
   [regression_type, regression_tag, ~, ~, ~, regression_cases] = ...
      icemodel.test.helpers.prepareBaselineBuild( ...
      "regression", "v1.1", "", "full", "icemodel", output_file, ...
      2016, [], "kanm", ["kanm"; "kanl"]);
   [perf_type, perf_tag, ~, ~, ~, perf_cases] = ...
      icemodel.test.helpers.prepareBaselineBuild( ...
      "perf", "v1.1", "", "full", "icemodel", output_file, ...
      2016, [], "kanm", ["kanm"; "kanl"]);

   testCase.verifyEqual([regression_type; perf_type], ...
      ["release"; "release"]);
   testCase.verifyEqual([regression_tag; perf_tag], ["v1.1"; "v1.1"]);
   testCase.verifyEqual( ...
      regression_cases.forcings, lower(regression_cases.sitename));
   testCase.verifyEqual(perf_cases.forcings, lower(perf_cases.sitename));
   testCase.verifyFalse(isfile(output_file));
end

function test_direct_release_build_rejects_existing_target(testCase)
   % Direct versioned builders share snapshot immutability and must stop at
   % setup before model execution or replacement of existing bytes.

   kinds = ["regression"; "perf"];
   for i = 1:numel(kinds)
      output_file = string(tempname) + ".mat";
      cleanup = onCleanup(@() deleteIfExists(output_file));
      marker = "unchanged";
      save(output_file, 'marker');

      testCase.verifyError(@() ...
         icemodel.test.helpers.prepareBaselineBuild( ...
         kinds(i), "v1.1", "", "full", "icemodel", output_file, ...
         2016, [], "kanm", ["kanm"; "kanl"]), ...
         'icemodel:test:releaseBaselineImmutable');
      saved = load(output_file, 'marker');
      testCase.verifyEqual(saved.marker, marker);
      clear cleanup
   end
end

function test_baseline_builders_forward_data_root_without_writing(testCase)
   % Intercept setup before model execution to prove both builders forward the
   % selected root without creating or accepting a baseline artifact.

   fixture_root = string(tempname);
   fixture_source = fullfile(fixture_root, "source");
   helper_dir = fullfile(fixture_source, "+icemodel", "+test", "+helpers");
   mkdir(helper_dir);
   output_file = fullfile(fixture_root, "must-not-exist.mat");

   env_names = ["ICEMODEL_EXPECTED_BUILDER_ARGUMENT_ROOT"; ...
      "ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT"; ...
      "ICEMODEL_EXPECTED_BUILDER_KIND"; ...
      "ICEMODEL_EXPECTED_BUILDER_CASENAME"; ...
      "ICEMODEL_TEST_DATA_ROOT"];
   env_values = arrayfun(@(name) string(getenv(name)), env_names);
   original_path = path;
   cleanup = onCleanup(@() restoreBuilderFixture( ...
      original_path, fixture_root, env_names, env_values));

   % The bootstrap stub checks its received DATA_ROOT. The resolver stub also
   % checks that the perf builder exported the root for nested TestCase setup.
   writeBuilderBootstrapStub(fullfile(helper_dir, ...
      "bootstrapTestEnvironment.m"));
   writeBuilderResolverStub(fullfile(helper_dir, ...
      "resolveRequestedSmbmodels.m"));
   addpath(fixture_source, '-begin');
   clear('icemodel.test.helpers.bootstrapTestEnvironment', ...
      'icemodel.test.helpers.resolveRequestedSmbmodels')

   selected_root = fullfile(fixture_root, "selected-data");
   setenv('ICEMODEL_EXPECTED_BUILDER_ARGUMENT_ROOT', selected_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT', selected_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_CASENAME', 'verification');

   % Regression dispatch stops inside the bootstrap seam after verifying the
   % forwarded root, before case construction or output-file handling.
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'regression');
   testCase.verifyError(@() build_regression_baseline( ...
      data_root=selected_root, output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');

   % Performance dispatch continues through bootstrap, then stops at model
   % selection after proving the nested TestCase environment carries the root.
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'perf');
   testCase.verifyError(@() build_perf_baseline( ...
      data_root=selected_root, output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');

   % A blank selector resolves through the verification case before the perf
   % TestCase environment is exported; neither builder may fall back to test/data.
   verification_root = fullfile(icemodel.internal.fullpath(), 'data');
   setenv('ICEMODEL_EXPECTED_BUILDER_ARGUMENT_ROOT', '');
   setenv('ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT', verification_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'regression');
   testCase.verifyError(@() build_regression_baseline( ...
      data_root="", output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'perf');
   testCase.verifyError(@() build_perf_baseline( ...
      data_root="", output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');

   % Frozen v1.1 defaults to the historical test tree for both builders.
   historical_root = fullfile(icemodel.internal.fullpath(), 'test', 'data');
   setenv('ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT', historical_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_CASENAME', 'test');
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'regression');
   testCase.verifyError(@() build_regression_baseline( ...
      baseline_tag="v1.1", data_root="", output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'perf');
   testCase.verifyError(@() build_perf_baseline( ...
      baseline_tag="v1.1", data_root="", output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');

   % An explicit root remains authoritative even for a release registration.
   setenv('ICEMODEL_EXPECTED_BUILDER_ARGUMENT_ROOT', selected_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT', selected_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'regression');
   testCase.verifyError(@() build_regression_baseline( ...
      baseline_tag="v1.1", data_root=selected_root, ...
      output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');
   testCase.verifyFalse(isfile(output_file));

   clear cleanup
end

function test_baseline_runners_select_registered_data_case_without_running(testCase)
   % Stop both public runners at bootstrap and inspect their selected config
   % case plus explicit-root forwarding before any model or artifact work.

   fixture_root = string(tempname);
   helper_dir = fullfile(fixture_root, "source", ...
      "+icemodel", "+test", "+helpers");
   mkdir(helper_dir);

   env_names = ["ICEMODEL_EXPECTED_RUNNER_ROOT"; ...
      "ICEMODEL_EXPECTED_RUNNER_CASENAME"; ...
      "ICEMODEL_REGRESSION_BASELINE"; "ICEMODEL_TEST_DATA_ROOT"];
   env_values = arrayfun(@(name) string(getenv(name)), env_names);
   original_path = path;
   cleanup = onCleanup(@() restoreRunnerFixture( ...
      original_path, fixture_root, env_names, env_values));

   writeRunnerBootstrapStub(fullfile(helper_dir, ...
      "bootstrapTestEnvironment.m"));
   addpath(fullfile(fixture_root, "source"), '-begin');
   clear('icemodel.test.helpers.bootstrapTestEnvironment')

   % Blank rolling calls retain the modern verification case.
   setenv('ICEMODEL_EXPECTED_RUNNER_ROOT', '');
   setenv('ICEMODEL_EXPECTED_RUNNER_CASENAME', 'verification');
   verifyRunnerBootstrap(testCase, "rolling", "");

   % Blank v1.1 calls select the historical formal-data case.
   setenv('ICEMODEL_EXPECTED_RUNNER_CASENAME', 'test');
   verifyRunnerBootstrap(testCase, "v1.1", "");

   % Direct unittest-class setup follows the same release registration.
   setenv('ICEMODEL_REGRESSION_BASELINE', 'v1.1');
   setenv('ICEMODEL_TEST_DATA_ROOT', '');
   regression_test = IcemodelRegressionTest();
   testCase.verifyError(@() regression_test.configureCases(), ...
      'icemodel:test:baselineRunnerRootObserved');

   % An explicit root remains authoritative while the registered case stays
   % visible to the central bootstrap contract.
   selected_root = fullfile(fixture_root, "selected-data");
   setenv('ICEMODEL_EXPECTED_RUNNER_ROOT', selected_root);
   verifyRunnerBootstrap(testCase, "v1.1", selected_root);
   clear cleanup
end

function test_release_runners_verify_registered_fixture_capability(testCase)
   % Both public v1.1 runners must stop at the same network-free capability
   % and hash check before they can dispatch a formal model case.

   data_root = string(tempname);
   mkdir(data_root)
   cleanup = onCleanup(@() rmdir(data_root, 's'));

   verifyFixtureCapabilityError(testCase, @() run_regression_suite( ...
      tier="smoke", smbmodel="icemodel", baseline="v1.1", ...
      data_root=data_root, build_report=false), data_root);
   verifyFixtureCapabilityError(testCase, @() run_perf_suite( ...
      tier="smoke", smbmodel="icemodel", baseline="v1.1", ...
      data_root=data_root, include_benchmarks=false, build_report=false), ...
      data_root);

   clear cleanup
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

function test_summarizeIce1Metrics_prefers_full_surface_residual(testCase)
   % Regression closure must use the diagnosed full SEB residual rather than
   % the algebraically closed shortwave-partition identity from postprocess.

   Time = datetime(2026, 1, 1) + hours(0:1)';
   balance = [3; -4];
   Qbal = zeros(2, 1);
   ice1 = timetable(Time, balance, Qbal);

   metrics = icemodel.test.helpers.summarizeIce1Metrics(ice1);
   testCase.verifyEqual(metrics.closure_seb_mae, 3.5);
   testCase.verifyEqual(metrics.closure_seb_rmse, sqrt(12.5), ...
      AbsTol=1e-12);
   testCase.verifyEqual(metrics.closure_seb_max_abs, 4);

   % Keep compatibility with older outputs that expose only Qbal.
   legacy = removevars(ice1, "balance");
   legacy.Qbal = [1; -2];
   legacy_metrics = icemodel.test.helpers.summarizeIce1Metrics(legacy);
   testCase.verifyEqual(legacy_metrics.closure_seb_max_abs, 2);
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

function verifyCaseForcingMatches(testCase, cases, baseline)
   %VERIFYCASEFORCINGMATCHES Match each case id and compare its forcing identity.

   testCase.assertTrue(ismember('forcings', ...
      baseline.Properties.VariableNames));
   for i = 1:height(cases)
      idx = icemodel.test.helpers.findCaseRow( ...
         baseline, cases.case_id(i));
      testCase.assertNotEmpty(idx, ...
         sprintf('frozen baseline missing case %s', cases.case_id(i)));
      testCase.verifyEqual( ...
         string(baseline.forcings(idx)), cases.forcings(i));
   end
end

function deleteIfExists(filepath)
   %DELETEIFEXISTS Remove one file if it exists.

   if exist(filepath, 'file') == 2
      delete(filepath);
   end
end

function writeBuilderBootstrapStub(filename)
   %WRITEBUILDERBOOTSTRAPSTUB Intercept baseline setup before model execution.

   lines = [ ...
      "function [rootdir, input_path, output_path, eval_path, cleanup] = bootstrapTestEnvironment(varargin)"
      "names = string(varargin(1:2:end));"
      "values = string(varargin(2:2:end));"
      "data_root = values(names == 'data_root');"
      "casename = values(names == 'icemodel_config_casename');"
      "expected = string(getenv('ICEMODEL_EXPECTED_BUILDER_ARGUMENT_ROOT'));"
      "resolved = string(getenv('ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT'));"
      "expected_case = string(getenv('ICEMODEL_EXPECTED_BUILDER_CASENAME'));"
      "if data_root ~= expected || casename ~= expected_case"
      "   error('icemodel:test:baselineDataRootMismatch', 'unexpected data root')"
      "end"
      "if strcmp(getenv('ICEMODEL_EXPECTED_BUILDER_KIND'), 'regression')"
      "   error('icemodel:test:baselineDataRootObserved', 'observed data root')"
      "end"
      "rootdir = ''; input_path = fullfile(resolved, 'input');"
      "output_path = ''; eval_path = '';"
      "cleanup = onCleanup(@() false);"
      "end"];
   writeTextFixture(filename, join(lines, newline));
end

function writeBuilderResolverStub(filename)
   %WRITEBUILDERRESOLVERSTUB Verify perf TestCase root propagation.

   lines = [ ...
      "function models = resolveRequestedSmbmodels(~)"
      "models = string.empty();"
      "expected = string(getenv('ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT'));"
      "if string(getenv('ICEMODEL_TEST_DATA_ROOT')) ~= expected"
      "   error('icemodel:test:baselineDataRootMismatch', 'nested root mismatch')"
      "end"
      "error('icemodel:test:baselineDataRootObserved', 'observed data root')"
      "end"];
   writeTextFixture(filename, join(lines, newline));
end

function writeRunnerBootstrapStub(filename)
   %WRITERUNNERBOOTSTRAPSTUB Stop public runners at their setup boundary.

   lines = [ ...
      "function [rootdir, input_path, output_path, eval_path, cleanup] = bootstrapTestEnvironment(varargin)"
      "names = string(varargin(1:2:end));"
      "values = string(varargin(2:2:end));"
      "data_root = values(names == 'data_root');"
      "casename = values(names == 'icemodel_config_casename');"
      "expected_root = string(getenv('ICEMODEL_EXPECTED_RUNNER_ROOT'));"
      "expected_case = string(getenv('ICEMODEL_EXPECTED_RUNNER_CASENAME'));"
      "if data_root ~= expected_root || casename ~= expected_case"
      "   error('icemodel:test:baselineRunnerRootMismatch', 'unexpected root policy')"
      "end"
      "rootdir = ''; input_path = ''; output_path = ''; eval_path = '';"
      "cleanup = onCleanup(@() false);"
      "error('icemodel:test:baselineRunnerRootObserved', 'observed root policy')"
      "end"];
   writeTextFixture(filename, join(lines, newline));
end

function verifyRunnerBootstrap(testCase, baseline, data_root)
   %VERIFYRUNNERBOOTSTRAP Check both public entrypoints at the shared seam.

   testCase.verifyError(@() run_regression_suite( ...
      baseline=baseline, data_root=data_root, build_report=false), ...
      'icemodel:test:baselineRunnerRootObserved');
   testCase.verifyError(@() run_perf_suite( ...
      baseline=baseline, data_root=data_root, include_benchmarks=false, ...
      build_report=false), ...
      'icemodel:test:baselineRunnerRootObserved');
end

function verifyFixtureCapabilityError(testCase, operation, data_root)
   %VERIFYFIXTURECAPABILITYERROR Check the stable network-free repair contract.
   try
      operation();
      testCase.verifyFail('expected incomplete frozen fixture capability');
   catch err
      testCase.verifyEqual(string(err.identifier), ...
         "icemodel:verification:fetchFixtures:incompleteFixtures");
      testCase.verifySubstring(string(err.message), "formal-core");
      testCase.verifySubstring(string(err.message), data_root);
   end
end

function writeTextFixture(filename, text)
   %WRITETEXTFIXTURE Write one temporary MATLAB source fixture.

   file_id = fopen(filename, 'w');
   assert(file_id ~= -1, 'could not create fixture: %s', filename)
   cleanup = onCleanup(@() fclose(file_id));
   fprintf(file_id, '%s', text);
   clear cleanup
end

function restoreBuilderFixture(original_path, fixture_root, names, values)
   %RESTOREBUILDERFIXTURE Restore path/env state and remove the temporary seam.

   path(original_path);
   clear('icemodel.test.helpers.bootstrapTestEnvironment', ...
      'icemodel.test.helpers.resolveRequestedSmbmodels')
   for n = 1:numel(names)
      setenv(names(n), values(n));
   end
   if isfolder(fixture_root)
      rmdir(fixture_root, 's');
   end
end

function restoreRunnerFixture(original_path, fixture_root, names, values)
   %RESTORERUNNERFIXTURE Restore path/env state after runner interception.

   path(original_path);
   clear('icemodel.test.helpers.bootstrapTestEnvironment')
   for n = 1:numel(names)
      setenv(names(n), values(n));
   end
   if isfolder(fixture_root)
      rmdir(fixture_root, 's');
   end
end
