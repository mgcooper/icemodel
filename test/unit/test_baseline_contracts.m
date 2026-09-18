function tests = test_baseline_contracts
   %TEST_BASELINE_CONTRACTS Verify baseline-selection and legacy-load helpers.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   %SETUPONCE Put the baseline runners on the path.
   %
   % build_regression_baseline and build_perf_baseline live in test/tools,
   % IcemodelRegressionTest in test/regression, and run_regression_suite
   % plus run_perf_suite at the test root. None of these folders is on the
   % path by default, so running this file on its own would otherwise
   % report "Undefined function" instead of the expected error.

   root = icemodel.internal.fullpath();
   original_path = path;
   testCase.addTeardown(@() path(original_path));
   folders = {fullfile(root, 'test'), ...
      fullfile(root, 'test', 'tools'), ...
      fullfile(root, 'test', 'regression')};
   for k = 1:numel(folders)
      if isfolder(folders{k})
         addpath(folders{k});
      end
   end
end

function test_aa_acceptance_rejects_non_independent_or_foreign_runs(testCase)
   % The gate certifies two independent runs from one environment. The
   % same file, a reused run_id, another host, or a session-isolated
   % artifact must all be rejected before any ratio is computed. The
   % no-argument mode runs two full timed suites, so these tests cover
   % the artifact-comparison branches only.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   case_summary = aaCaseSummary(100.0, true);

   meta = aaArtifactMeta("run-a");
   file_a = fullfile(fixture.Folder, 'perf_results_a.mat');
   save(file_a, 'meta', 'case_summary');

   % The same file for A and B reproduces trivially: rejected.
   testCase.verifyError(@() run_aa_acceptance(file_a, file_a), ...
      'icemodel:test:aaAcceptance:reusedFile');

   % A copy of the same run under another file name carries the same
   % run_name: rejected.
   file_b = fullfile(fixture.Folder, 'perf_results_b.mat');
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:reusedRun');

   % A blank run name cannot identify an independent timed run.
   meta = aaArtifactMeta("   ");
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:unknownRunName');

   % The two sides must contain the same number of model artifacts.
   testCase.verifyError(@() run_aa_acceptance( ...
      [string(file_a), string(file_b)], file_b), ...
      'icemodel:test:aaAcceptance:artifactCountMismatch');

   % A side cannot repeat one artifact in place of another model artifact.
   testCase.verifyError(@() run_aa_acceptance( ...
      [string(file_a), string(file_a)], ...
      [string(file_b), string(file_b)]), ...
      'icemodel:test:aaAcceptance:duplicateArtifact');

   % An artifact without measured cases cannot certify the timing protocol.
   empty_summary = case_summary([], :);
   meta = aaArtifactMeta("empty-a");
   empty_a = fullfile(fixture.Folder, 'perf_results_empty_a.mat');
   case_summary = empty_summary;
   save(empty_a, 'meta', 'case_summary');
   meta = aaArtifactMeta("empty-b");
   empty_b = fullfile(fixture.Folder, 'perf_results_empty_b.mat');
   save(empty_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(empty_a, empty_b), ...
      'icemodel:test:aaAcceptance:emptyCaseSummary');
   case_summary = aaCaseSummary(100.0, true);

   % Another machine is not comparable: rejected.
   meta = aaArtifactMeta("run-b");
   meta.hostname = "other-host";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:environmentMismatch');

   % A session-isolated artifact certifies a different protocol: rejected.
   meta = aaArtifactMeta("run-b");
   meta.isolation = "session";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:notProcessIsolated');

   % An artifact without a recorded hostname certifies an unverified
   % machine: rejected.
   meta = aaArtifactMeta("run-b");
   meta.hostname = "   ";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:unknownHost');

   % A blank MATLAB version cannot identify one measurement environment.
   meta = aaArtifactMeta("run-b");
   meta.matlab_version = "   ";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:unknownMatlabVersion');

   % Two runs of different code certify an A/B change: rejected.
   meta = aaArtifactMeta("run-b");
   meta.git_revision = "other-rev";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:revisionMismatch');

   % A blank revision cannot certify that both runs used the same source.
   meta = aaArtifactMeta("run-b");
   meta.git_revision = "   ";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:unknownRevision');

   % Two measurement procedures are not comparable: rejected.
   meta = aaArtifactMeta("run-b");
   meta.n_runs = 30;
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:procedureMismatch');

   % Two input data trees are not comparable: rejected.
   meta = aaArtifactMeta("run-b");
   meta.data_root = "other-data";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:inputMismatch');

   % A blank data root cannot identify the measured input tree.
   meta = aaArtifactMeta("run-b");
   meta.data_root = "   ";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:unknownDataRoot');

   % A reused case ID can name a different forcing product: rejected.
   meta = aaArtifactMeta("run-b");
   case_summary.forcings = "legacy_forcing";
   save(file_b, 'meta', 'case_summary');
   testCase.verifyError(@() run_aa_acceptance(file_a, file_b), ...
      'icemodel:test:aaAcceptance:caseIdentityMismatch');
   case_summary.forcings = "promice_filled";

   % A side that mixes two runs could pair artifacts whose differences
   % cancel in the ratios: rejected.
   meta = aaArtifactMeta("run-c");
   file_a2 = fullfile(fixture.Folder, 'perf_results_a2.mat');
   save(file_a2, 'meta', 'case_summary');
   meta = aaArtifactMeta("run-b");
   save(file_b, 'meta', 'case_summary');
   file_b2 = fullfile(fixture.Folder, 'perf_results_b2.mat');
   save(file_b2, 'meta', 'case_summary');
   testCase.verifyError( ...
      @() run_aa_acceptance([string(file_a), string(file_a2)], ...
      [string(file_b), string(file_b2)]), ...
      'icemodel:test:aaAcceptance:mixedSide');
end

function test_aa_acceptance_no_argument_mode_compares_two_runs(testCase)
   % The no-argument form runs the formal suite twice and compares the
   % two runs' artifacts. A stub run_perf_suite shadows the real one and
   % writes one band-centered artifact per call, so the branch runs
   % without timing anything.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   stub_dir = fullfile(fixture.Folder, 'stub');
   icemodel.helpers.ensureDirExists(stub_dir);
   writeAaSuiteStub(fullfile(stub_dir, 'run_perf_suite.m'));
   testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stub_dir));

   [~] = evalc('report = run_aa_acceptance();');
   testCase.verifyTrue(report.passed);
   testCase.verifyEqual(height(report.cases), 1);
   testCase.verifyEqual(report.cases.ratio_b_over_a, 1.0, 'AbsTol', 0);
end

function test_aa_acceptance_pairs_artifacts_by_filename(testCase)
   % Pair model artifacts by filename when their parent directories differ.

   [status, scratch_root] = system('mktemp -d');
   testCase.assertEqual(status, 0);
   scratch_root = string(strtrim(scratch_root));
   cleanup = onCleanup(@() rmdir(scratch_root, 's'));
   testCase.addTeardown(@() delete(cleanup));

   dirs = fullfile(scratch_root, ["a_ice", "z_skin", "a_skin", "z_ice"]);
   arrayfun(@mkdir, dirs);

   meta = aaArtifactMeta("run-a");
   case_summary = aaCaseSummary(100.0, true);
   ice_a = fullfile(dirs(1), "perf_icemodel.mat");
   save(ice_a, 'meta', 'case_summary');
   case_summary.case_id = "skinmodel_kanm_2016_solver1";
   skin_a = fullfile(dirs(2), "perf_skinmodel.mat");
   save(skin_a, 'meta', 'case_summary');

   meta = aaArtifactMeta("run-b");
   skin_b = fullfile(dirs(3), "perf_skinmodel.mat");
   save(skin_b, 'meta', 'case_summary');
   case_summary.case_id = "icemodel_kanm_2016_solver1";
   ice_b = fullfile(dirs(4), "perf_icemodel.mat");
   save(ice_b, 'meta', 'case_summary');

   [~] = evalc(['report = run_aa_acceptance(' ...
      '[skin_a, ice_a], [ice_b, skin_b]);']);
   testCase.verifyTrue(report.passed);
   testCase.verifyEqual(height(report.cases), 2);
end

function test_aa_acceptance_fails_invalid_or_unstable_measurements(testCase)
   % An invalid measurement or an ambient-unstable run must fail the
   % verdict even when every ratio is inside the band.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   case_summary = aaCaseSummary(100.0, true);

   meta = aaArtifactMeta("run-a");
   file_a = fullfile(fixture.Folder, 'perf_results_a.mat');
   save(file_a, 'meta', 'case_summary');

   % Run B's measurement failed its validity gate: the verdict fails.
   meta = aaArtifactMeta("run-b");
   case_summary.valid = false;
   file_b = fullfile(fixture.Folder, 'perf_results_b.mat');
   save(file_b, 'meta', 'case_summary');
   [~] = evalc('report = run_aa_acceptance(file_a, file_b);');
   testCase.verifyFalse(report.passed);
   testCase.verifyFalse(all(report.cases.measurement_valid));

   % Run B is ambient-unstable: the verdict fails.
   meta.ambient_stable = false;
   case_summary.valid = true;
   save(file_b, 'meta', 'case_summary');
   [~] = evalc('report = run_aa_acceptance(file_a, file_b);');
   testCase.verifyFalse(report.passed);
   testCase.verifyFalse(report.ambient_stable);
end

function test_perf_builder_refuses_benchmarks_in_session_mode(testCase)
   % The component benchmark suite runs before the model cases, so an
   % in-session build with benchmarks would measure model medians in a
   % session the benchmarks already warmed. The builder must refuse. A
   % custom output file keeps the probe past the managed-build isolation
   % gate, which refuses every session build of the managed file first.

   testCase.verifyError(@() build_perf_baseline( ...
      isolation="session", include_benchmarks=true, ...
      output_file=string(tempname) + ".mat"), ...
      'icemodel:test:perf:benchmarksWarmSession');
end

function test_aa_acceptance_accepts_the_band_edge(testCase)
   % The A/A band is a closed interval: a ratio exactly at the band edge
   % certifies the measurement protocol, and one just outside fails it.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   band = 1 + icemodel.test.helpers.perfMeasurementPolicy().tol_perf;

   % Two one-case artifacts whose medians sit exactly at the band edge.
   meta = aaArtifactMeta("run-a");
   case_summary = aaCaseSummary(100.0, true);
   file_a = fullfile(fixture.Folder, 'perf_results_a.mat');
   save(file_a, 'meta', 'case_summary');
   meta = aaArtifactMeta("run-b");
   case_summary.median_wall_s = 100.0 * band;
   file_b = fullfile(fixture.Folder, 'perf_results_b.mat');
   save(file_b, 'meta', 'case_summary');

   % Discard the printed verdict so the suite log stays quiet.
   [~] = evalc('report = run_aa_acceptance(file_a, file_b);');
   testCase.verifyTrue(report.passed);
   testCase.verifyEqual(report.cases.ratio_b_over_a, band, 'AbsTol', 0);

   % A ratio just outside the closed interval fails.
   case_summary.median_wall_s = 100.0 * band * (1 + 1e-6);
   save(file_b, 'meta', 'case_summary');
   [~] = evalc('report = run_aa_acceptance(file_a, file_b);');
   testCase.verifyFalse(report.passed);
end

function test_resolveBaselineSelector_handles_rolling_and_release(testCase)
   % Test the public selector with rolling and release baselines.

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
   legacy_time = convertTo(datetime(2026, 3, 9, 4, 15, 39, ...
      'TimeZone', 'UTC'), 'datenum');

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
   % The four immutable release files must add solver from solver_mode in
   % memory. Each file must pass its forcing and case checks.
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
   % The loader must not replace baseline-file values with selector metadata
   % or an alias copied from another column.
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
   % The performance matrix must keep skinmodel at solver 1 when its solver
   % input is 2. This keeps build, run, and bootstrap paths consistent.

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

function test_solver_namelist_includes_solver_zero(testCase)
   % Solver 0, the Dirichlet single sweep, is a supported solver id. The
   % shared solver filter accepts it and rejects an unsupported id.

   returned = icemodel.namelists.solver();
   expected = [0 1 2 3];
   testCase.verifyEqual(returned, expected);
   testCase.verifyWarningFree( ...
      @() icemodel.validators.mustBeSolverFilter([0 3]));
   testCase.verifyError( ...
      @() icemodel.validators.mustBeSolverFilter(4), ...
      'icemodel:validators:mustBeSolverFilter');
end

function test_formal_case_matrices_register_solvers_per_baseline(testCase)
   % Rolling matrices run icemodel solvers 0 to 3. The frozen v1.1 and v1.2
   % matrices run solvers 1 to 3, the ids of their accepted rows. The frozen
   % v1.3 matrices run solvers 0 to 3, because the rolling rows it froze
   % carried solver 0.

   selectors = ["rolling", "v1.1", "v1.2", "v1.3"];
   expected = {[0; 1; 2; 3], [1; 2; 3], [1; 2; 3], [0; 1; 2; 3]};
   for k = 1:numel(selectors)
      regression_cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
         tier="smoke", baseline=selectors(k));
      perf_cases = icemodel.test.helpers.getPerfCaseMatrix( ...
         tier="smoke", baseline=selectors(k));
      verifyMatrixSolvers(testCase, regression_cases, expected{k});
      verifyMatrixSolvers(testCase, perf_cases, expected{k});

      % The policy registers the same ids that the matrices run.
      returned = icemodel.test.helpers.formalBaselinePolicy( ...
         selectors(k)).icemodel_solvers;
      testCase.verifyEqual(returned, expected{k}');
   end
end

function test_formal_case_ids_and_loader_accept_solver_zero(testCase)
   % A solver 0 case id and a saved solver 0 row are valid. A saved id outside
   % the supported solvers, or a fractional id, is a schema error.

   returned = icemodel.test.helpers.makeFormalCaseId( ...
      "icemodel", "kanm", 2016, 0);
   testCase.verifyEqual(returned, "icemodel_kanm_2016_solver0");
   testCase.verifyError(@() icemodel.test.helpers.makeFormalCaseId( ...
      "icemodel", "kanm", 2016, 4), ...
      'icemodel:validators:mustBeSolverFilter');

   filepath = [tempname '.mat'];
   cleanup = onCleanup(@() deleteIfExists(filepath));
   RegressionBaseline = table("icemodel_kanm_2016_solver0", 0, 1.5, ...
      'VariableNames', {'case_id', 'solver', 'runoff_final'});
   save(filepath, 'RegressionBaseline');
   loaded = icemodel.test.helpers.loadBaseline("regression", ...
      smbmodel="icemodel", filename=filepath);
   testCase.verifyEqual(loaded.solver, 0);

   % Solver 4 is not a supported id, and 1.5 is not an integer id.
   for solver_id = [4, 1.5]
      RegressionBaseline.solver = solver_id;
      save(filepath, 'RegressionBaseline');
      testCase.verifyError(@() icemodel.test.helpers.loadBaseline( ...
         "regression", filename=filepath), ...
         'icemodel:test:baselineSolverSchemaMismatch');
   end
   clear cleanup
end

function test_spectral_tools_share_the_solver_validator(testCase)
   % The spectral study tools validate their solver input with the shared
   % solver filter, so an unsupported id fails before any model run.

   tools = {@plot_spectral_variant_profiles, @run_spectral_study_bootstrap, ...
      @summarize_spectral_density_floor, @summarize_spectral_perf};
   for k = 1:numel(tools)
      testCase.verifyError(@() tools{k}('solver', 4), ...
         'icemodel:validators:mustBeSolverFilter', func2str(tools{k}));
   end
end

function test_getFormalForcing_returns_official_product(testCase)
   % Keep the shared formal forcing identity explicit and independently tested.

   testCase.verifyEqual( ...
      icemodel.test.helpers.getFormalForcing(), "promice_filled");
end

function test_formal_baseline_policy_owns_default_data_case(testCase)
   % The same registration selects forcing identity and its default data tree.

   rolling = icemodel.test.helpers.formalBaselinePolicy("rolling");
   release_v11 = icemodel.test.helpers.formalBaselinePolicy("v1.1");
   release_v12 = icemodel.test.helpers.formalBaselinePolicy("v1.2");
   release_v12_alias = icemodel.test.helpers.formalBaselinePolicy("V1_2");
   release_v13 = icemodel.test.helpers.formalBaselinePolicy("v1.3");

   testCase.verifyEqual(rolling.config_case, "verification");
   testCase.verifyEqual(rolling.forcing, "promice_filled");
   testCase.verifyEmpty(rolling.required_fixture_capabilities);
   testCase.verifyFalse(rolling.use_fixture_root_for_model);
   testCase.verifyEqual(rolling.promice_filled_policy_sha256, ...
      icemodel.forcing.reconstruct.policySha256());
   testCase.verifyEqual(release_v11.config_case, "test");
   testCase.verifyEqual(release_v11.site_forcings, ["kanm"; "kanl"]);
   testCase.verifyEqual( ...
      release_v11.required_fixture_capabilities, "formal-core");
   testCase.verifyFalse(release_v11.use_fixture_root_for_model);
   testCase.verifyFalse(release_v11.snapshot_from_rolling);
   testCase.verifyEqual(release_v12.config_case, "verification");
   testCase.verifyEqual(release_v12.forcing_mode, "fixed");
   testCase.verifyEqual(release_v12.forcing, "promice_filled");
   testCase.verifyEmpty(release_v12.sites);
   testCase.verifyEmpty(release_v12.site_forcings);
   testCase.verifyEqual( ...
      release_v12.required_fixture_capabilities, "formal-core");
   testCase.verifyTrue(release_v12.use_fixture_root_for_model);
   testCase.verifyTrue(release_v12.snapshot_from_rolling);
   testCase.verifyEqual(release_v12.promice_filled_policy_sha256, ...
      "bd336da0880474f1987facc2311c4f45a6c281877ae8b944a3fbdc7cfb68d513");
   testCase.verifyEqual(release_v12_alias.baseline_tag, "v1.2");
   testCase.verifyEqual(release_v13.baseline_tag, "v1.3");
   testCase.verifyEqual(release_v13.config_case, "verification");
   testCase.verifyEqual(release_v13.forcing, "promice_filled");
   testCase.verifyEqual( ...
      release_v13.required_fixture_capabilities, "formal-core");
   testCase.verifyTrue(release_v13.use_fixture_root_for_model);
   testCase.verifyTrue(release_v13.snapshot_from_rolling);
   testCase.verifyEqual(release_v13.icemodel_solvers, [0 1 2 3]);
   % The v1.3 pin is the policy digest at registration, so it must equal
   % the live value until POLICY.md changes after the release.
   testCase.verifyEqual(release_v13.promice_filled_policy_sha256, ...
      icemodel.forcing.reconstruct.policySha256());
   testCase.verifyFalse(isfield(release_v13, 'require_source_revision'));
end

function test_release_v12_case_matrices_use_rolling_forcing(testCase)
   % v1.2 snapshots use the accepted forcing from the rolling baselines.

   regression_cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
      tier="full", baseline="v1.2");
   perf_cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier="full", baseline="v1_2");

   testCase.verifyEqual(unique(regression_cases.forcings), ...
      "promice_filled");
   testCase.verifyEqual(unique(perf_cases.forcings), "promice_filled");

   % These rows check the plumbing that carries the registered digest into
   % each matrix and into opts, so read the digest from its registration
   % rather than declaring it a second time. The pin itself is checked in
   % test_formal_baseline_policy_owns_default_data_case.
   expected_sha256 = icemodel.test.helpers.formalBaselinePolicy( ...
      "v1.2").promice_filled_policy_sha256;
   testCase.verifyEqual( ...
      unique(regression_cases.promice_filled_expected_policy_sha256), ...
      expected_sha256);
   testCase.verifyEqual( ...
      unique(perf_cases.promice_filled_expected_policy_sha256), ...
      expected_sha256);
   opts = icemodel.test.helpers.setModelOptsForCase(regression_cases(1, :));
   testCase.verifyEqual( ...
      string(opts.promice_filled_expected_policy_sha256), expected_sha256);
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
   % inherit the rolling product or the v1.1 station aliases by default.

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
   % A case id must not let the metrics of one site pass as another site's.

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
end

function test_incompatible_performance_verdict_is_not_acceptance(testCase)
   % Valid samples do not pass when baseline timings are not comparable.
   [passed, ~, ~, ~, reason] = ...
      icemodel.test.helpers.formalPerformanceVerdict( ...
      true, 10, table(), [], false, 0.2, "metadata incompatible");
   testCase.verifyFalse(passed);
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
   % while an unknown nonnumeric column must still be rejected.
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

function test_bootstrap_release_rejects_partial_model_snapshots(testCase)
   % A partial aggregate release keeps the existing file and creates nothing.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   models = ["icemodel"; "skinmodel"];
   existing_file = fullfile(fixture.Folder, "icemodel.mat");
   missing_file = fullfile(fixture.Folder, "skinmodel.mat");
   writeTextFixture(existing_file, "frozen baseline bytes");
   original_bytes = fileread(existing_file);
   snapshot_calls = string.empty(0, 1);
   removed_models = string.empty(0, 1);
   fail_model = "";
   fail_load = false;
   empty_reload = false;
   force_empty = false;
   remove_fail_model = "";
   source_revisions = ["revision-a"; "revision-a"];

   testCase.verifyError(@() ...
      icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub), ...
      'icemodel:test:partialReleaseBaseline');

   testCase.verifyEmpty(snapshot_calls);
   testCase.verifyEqual(fileread(existing_file), original_bytes);
   testCase.verifyFalse(isfile(missing_file));

   % An existing file that filters to no rows must not be deleted on failure.
   writeTextFixture(missing_file, "second frozen baseline");
   second_bytes = fileread(missing_file);
   force_empty = true;
   testCase.verifyError(@() ...
      icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub, ...
      snapshot_remover=@removeStub), ...
      'test:immutableSnapshot');
   testCase.verifyEmpty(removed_models);
   testCase.verifyEqual(fileread(existing_file), original_bytes);
   testCase.verifyEqual(fileread(missing_file), second_bytes);

   % A failure while creating a fresh set removes completed snapshots.
   delete(existing_file);
   delete(missing_file);
   snapshot_calls = string.empty(0, 1);
   force_empty = false;
   fail_model = "skinmodel";
   testCase.verifyError(@() ...
      icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub, ...
      snapshot_remover=@removeStub), ...
      'test:snapshotFailure');
   testCase.verifyEqual(snapshot_calls, models);
   testCase.verifyEqual(removed_models, "icemodel");
   testCase.verifyFalse(isfile(existing_file));
   testCase.verifyFalse(isfile(missing_file));

   % Reload failures roll back both files even when one remover fails.
   snapshot_calls = string.empty(0, 1);
   removed_models = string.empty(0, 1);
   fail_model = "";
   fail_load = true;
   remove_fail_model = "icemodel";
   testCase.verifyError(@() ...
      icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub, ...
      snapshot_remover=@removeStub), ...
      'test:snapshotLoadFailure');
   testCase.verifyEqual(snapshot_calls, models);
   testCase.verifyEqual(removed_models, models);
   testCase.verifyFalse(isfile(existing_file));
   testCase.verifyFalse(isfile(missing_file));

   % An empty persisted reload is incomplete and rolls back the whole set.
   snapshot_calls = string.empty(0, 1);
   removed_models = string.empty(0, 1);
   fail_load = false;
   empty_reload = true;
   remove_fail_model = "";
   testCase.verifyError(@() ...
      icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub, ...
      snapshot_remover=@removeStub), ...
      'icemodel:test:emptyReleaseSnapshot');
   testCase.verifyEqual(snapshot_calls, models);
   testCase.verifyEqual(removed_models, models);
   testCase.verifyFalse(isfile(existing_file));
   testCase.verifyFalse(isfile(missing_file));

   % A fresh aggregate request creates both scalar model snapshots.
   snapshot_calls = string.empty(0, 1);
   removed_models = string.empty(0, 1);
   fail_load = false;
   empty_reload = false;
   remove_fail_model = "";
   baseline = icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub);
   testCase.verifyEqual(snapshot_calls, models);
   testCase.verifyEqual(sort(baseline.smbmodel), sort(models));

   % Existing aggregate snapshots load without calling a snapshotter.
   snapshot_calls = string.empty(0, 1);
   baseline = icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub);
   testCase.verifyEmpty(snapshot_calls);
   testCase.verifyEqual(sort(baseline.smbmodel), sort(models));

   % Git is the provenance of a tracked release set, so the loader does
   % not compare the recorded revisions of the two model files.
   source_revisions(2) = "revision-b";
   baseline = icemodel.test.helpers.resolveBootstrapRelease( ...
      "perf", "v1.2", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, perf_snapshotter=@snapshotStub);
   testCase.verifyEqual(sort(baseline.smbmodel), sort(models));
   source_revisions(2) = "revision-a";

   % A row-vector regression request uses the regression snapshotter.
   delete(existing_file);
   delete(missing_file);
   baseline = icemodel.test.helpers.resolveBootstrapRelease( ...
      "regression", "v1.2", models.', 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', true), ...
      loader=@loadStub, regression_snapshotter=@snapshotStub);
   testCase.verifyEqual(snapshot_calls, models);
   testCase.verifyEqual(sort(baseline.smbmodel), sort(models));

   % A preserved release reports a missing immutable model snapshot.
   delete(existing_file);
   delete(missing_file);
   testCase.verifyError(@() ...
      icemodel.test.helpers.resolveBootstrapRelease( ...
      "regression", "v1.1", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', false), ...
      loader=@loadStub, regression_snapshotter=@snapshotStub), ...
      'icemodel:test:preservedReleaseMissing');

   % Historical v1.1 files remain usable without source-revision metadata.
   writeTextFixture(existing_file, "historical icemodel baseline");
   writeTextFixture(missing_file, "historical skinmodel baseline");
   source_revisions(:) = "";
   baseline = icemodel.test.helpers.resolveBootstrapRelease( ...
      "regression", "v1.1", "all", 2016, ...
      policy_resolver=@(~) struct('snapshot_from_rolling', false), ...
      loader=@loadStub, regression_snapshotter=@snapshotStub);
   testCase.verifyEqual(sort(baseline.smbmodel), sort(models));

   function [baseline, meta] = loadStub(kind, varargin)
      %LOADSTUB Load temporary model markers as a baseline table.
      kwargs = struct(varargin{:});
      model_index = find(models == kwargs.smbmodel, 1);
      meta = struct('git_revision', source_revisions(model_index));
      assert(ismember(kind, ["perf", "regression"]) ...
         && ismember(kwargs.baseline_tag, ["rolling", "v1.1", "v1.2"]) ...
         && kwargs.simyear == 2016)
      if kwargs.baseline_tag == "rolling"
         baseline = table(kwargs.smbmodel, ...
            'VariableNames', {'smbmodel'});
         return
      end
      if force_empty
         baseline = table(strings(0, 1), ...
            'VariableNames', {'smbmodel'});
         return
      end
      all_present = all(arrayfun(@(model) isfile( ...
         fullfile(fixture.Folder, model + ".mat")), models));
      if fail_load && all_present
         error('test:snapshotLoadFailure', ...
            'synthetic snapshot reload failure')
      end
      if empty_reload && all_present
         baseline = table(strings(0, 1), ...
            'VariableNames', {'smbmodel'});
         return
      end
      present = arrayfun(@(model) isfile( ...
         fullfile(fixture.Folder, model + ".mat")), models);
      selected = models(present);
      if kwargs.smbmodel ~= "all"
         selected = selected(selected == kwargs.smbmodel);
      end
      baseline = table(selected, 'VariableNames', {'smbmodel'});
   end

   function snapshotStub(varargin)
      %SNAPSHOTSTUB Record and create one temporary model snapshot.
      kwargs = struct(varargin{:});
      assert(kwargs.baseline_tag == "v1.2")
      if isfield(kwargs, 'simyear')
         assert(kwargs.simyear == 2016)
      end
      snapshot_calls(end + 1, 1) = kwargs.smbmodel;
      pathname = fullfile(fixture.Folder, kwargs.smbmodel + ".mat");
      if isfile(pathname)
         error('test:immutableSnapshot', ...
            'synthetic immutable target rejection')
      end
      if kwargs.smbmodel == fail_model
         error('test:snapshotFailure', 'synthetic snapshot failure')
      end
      writeTextFixture(fullfile(fixture.Folder, ...
         kwargs.smbmodel + ".mat"), "new baseline bytes");
   end

   function removeStub(~, ~, model, ~)
      %REMOVESTUB Record and remove one attempted temporary snapshot.
      removed_models(end + 1, 1) = model;
      pathname = fullfile(fixture.Folder, model + ".mat");
      if isfile(pathname)
         delete(pathname)
      end
      if model == remove_fail_model
         error('test:removeFailure', 'synthetic snapshot removal failure')
      end
   end
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

function test_clean_snapshot_worktree_check_ignores_own_outputs(testCase)
   % A versioned snapshot refuses a dirty tree, but the release files the
   % snapshot sequence itself writes must not count, so running the perf
   % tool and then the regression tool succeeds.
   own_files = [ ...
      "?? test/baselines/perf_baseline_2016_v1_3_icemodel.mat"; ...
      "?? test/baselines/perf_baseline_2016_v1_3_skinmodel.mat"; ...
      "?? test/baselines/regression_baseline_v1_3_icemodel.mat"];
   % git status -z ends every record with NUL.
   porcelain = @(records) char(join(records, char(0)) + char(0));
   runner = @(~) deal(0, porcelain(own_files));
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=runner));

   % A clean tree passes.
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, '')));

   % Any other tracked change or untracked file refuses.
   dirty = [own_files; " M icemodel/icemodel.m"];
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, porcelain(dirty))), ...
      'icemodel:test:snapshot:dirtyWorktree');
   untracked = "?? test/unit/new_test.m";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, porcelain(untracked))), ...
      'icemodel:test:snapshot:dirtyWorktree');

   % Another release's files are not this snapshot's outputs.
   other_release = "?? test/baselines/regression_baseline_v1_4_icemodel.mat";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, porcelain(other_release))), ...
      'icemodel:test:snapshot:dirtyWorktree');

   % An untracked name that contains " -> " is a file, not a rename, and a
   % rename record carries its original path in a second NUL record.
   arrow_name = "?? scratch -> " + own_files(1).extractAfter(3);
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, porcelain(arrow_name))), ...
      'icemodel:test:snapshot:dirtyWorktree');
   % A rename is a tracked change, so it refuses even when its new name is
   % a release file, and the parser reports it once, without the original
   % path as a second entry.
   rename = ["R  " + own_files(1).extractAfter(3); "old/name.mat"];
   try
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
         command_runner=@(~) deal(0, porcelain(rename)));
      testCase.verifyFail("a rename record must refuse the snapshot");
   catch err
      testCase.verifyEqual(err.identifier, ...
         'icemodel:test:snapshot:dirtyWorktree');
      testCase.verifySubstring(err.message, "R  test/baselines/");
      testCase.verifyFalse(contains(err.message, "old/name.mat"));
   end

   % The perf tool may run for another benchmark year, and the regression
   % tool that follows names no year, so any four-digit year is an output.
   other_year = "?? test/baselines/perf_baseline_2017_v1_3_skinmodel.mat";
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, porcelain(other_year))));

   % A tracked release file that is modified or deleted is not an output of
   % the sequence; it is a change to an immutable file.
   deleted = " D test/baselines/perf_baseline_2016_v1_3_icemodel.mat";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, porcelain(deleted))), ...
      'icemodel:test:snapshot:dirtyWorktree');
   modified = " M test/baselines/regression_baseline_v1_3_icemodel.mat";
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(0, porcelain(modified))), ...
      'icemodel:test:snapshot:dirtyWorktree');

   % Another checkout root prints the same relative names, so its own
   % outputs are still recognized.
   other_root = string(tempname);
   testCase.verifyWarningFree(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      repo_root=other_root, command_runner=runner));

   % A failed git command is an unknown state, not a clean one.
   testCase.verifyError(@() ...
      icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3", ...
      command_runner=@(~) deal(128, 'fatal')), ...
      'icemodel:test:snapshot:worktreeStatusUnavailable');
end

function test_snapshot_refuses_a_dirty_rolling_source(testCase)
   % A rolling file built on a dirty tree records a revision no commit
   % names, so it cannot become a release file.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   source_file = fullfile(fixture.Folder, 'perf_rolling.mat');
   PerfBaseline = table();
   meta = struct('git_revision', "v1.0.0-500-gabcdef12-dirty");
   save(source_file, 'PerfBaseline', 'meta');
   testCase.verifyError(@() icemodel.test.helpers.snapshotBaseline( ...
      "perf", "vNext", "icemodel", false, string.empty(), 2016, ...
      source_file), 'icemodel:test:releaseBaselineSourceDirty');
end

function test_direct_regression_class_run_refuses_release_without_root(testCase)
   % A release selector whose policy runs the model from the provisioned
   % fixture root needs an explicit data root. A direct class run bypasses
   % run_regression_suite, which owns that resolution, so it must refuse.
   release_policy = icemodel.test.helpers.formalBaselinePolicy("v1.2");
   rolling_policy = icemodel.test.helpers.formalBaselinePolicy("rolling");
   testCase.verifyError(@() ...
      IcemodelRegressionTest.assertReleaseRootProvided(release_policy, ""), ...
      'icemodel:test:regression:releaseRootRequired');
   testCase.verifyWarningFree(@() ...
      IcemodelRegressionTest.assertReleaseRootProvided( ...
      release_policy, "/some/provisioned/root"));
   testCase.verifyWarningFree(@() ...
      IcemodelRegressionTest.assertReleaseRootProvided(rolling_policy, ""));
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
   setup_dir = fullfile(fixture_source, ...
      "+icemodel", "+verification", "+setup");
   mkdir(helper_dir);
   mkdir(setup_dir);
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
   writeBuilderFetchStub(fullfile(setup_dir, "fetchFixtures.m"));
   addpath(fixture_source, '-begin');
   clear('icemodel.test.helpers.bootstrapTestEnvironment', ...
      'icemodel.test.helpers.resolveRequestedSmbmodels', ...
      'icemodel.verification.setup.fetchFixtures')

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

   % Frozen v1.2 verifies and runs from the registered fixture root.
   setenv('ICEMODEL_EXPECTED_BUILDER_ARGUMENT_ROOT', historical_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_CASENAME', 'verification');
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'regression');
   testCase.verifyError(@() build_regression_baseline( ...
      baseline_tag="v1.2", data_root="", output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');
   setenv('ICEMODEL_EXPECTED_BUILDER_KIND', 'perf');
   testCase.verifyError(@() build_perf_baseline( ...
      baseline_tag="v1.2", data_root="", output_file=output_file), ...
      'icemodel:test:baselineDataRootObserved');

   % An explicit root still takes precedence, even for a release registration.
   setenv('ICEMODEL_EXPECTED_BUILDER_ARGUMENT_ROOT', selected_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_RESOLVED_ROOT', selected_root);
   setenv('ICEMODEL_EXPECTED_BUILDER_CASENAME', 'test');
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
      "ICEMODEL_REGRESSION_BASELINE"; "ICEMODEL_TEST_DATA_ROOT"; ...
      "ICEMODEL_TEST_SESSION_ACTIVITY"];
   env_values = arrayfun(@(name) string(getenv(name)), env_names);
   original_path = path;
   cleanup = onCleanup(@() restoreRunnerFixture( ...
      original_path, fixture_root, env_names, env_values));

   % This probe aborts the runners at bootstrap and times nothing, so the
   % contaminated-session refusal does not apply; present a clean session.
   setenv('ICEMODEL_TEST_SESSION_ACTIVITY', '');

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

   % An explicit root still takes precedence, while the registered case stays
   % visible to the central bootstrap.
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

   % This probe stops at the capability check and times nothing, so the
   % contaminated-session refusal does not apply; present a clean session.
   prior_activity = getenv('ICEMODEL_TEST_SESSION_ACTIVITY');
   activity_cleanup = onCleanup(@() ...
      setenv('ICEMODEL_TEST_SESSION_ACTIVITY', prior_activity));
   setenv('ICEMODEL_TEST_SESSION_ACTIVITY', '');

   verifyFixtureCapabilityError(testCase, @() run_regression_suite( ...
      tier="smoke", smbmodel="icemodel", baseline="v1.1", ...
      data_root=data_root, build_report=false), ...
      data_root);
   verifyFixtureCapabilityError(testCase, @() run_perf_suite( ...
      tier="smoke", smbmodel="icemodel", baseline="v1.1", ...
      data_root=data_root, ...
      include_benchmarks=false, build_report=false), ...
      data_root);

   clear cleanup
end

function test_release_data_roots_preserve_v1_1_and_unify_v1_2(testCase)
   % v1.1 keeps its explicit root; v1.2 verifies the root it executes.
   v1_1 = icemodel.test.helpers.formalBaselinePolicy("v1.1");
   [data_root, fixture_root] = ...
      icemodel.test.helpers.resolveReleaseDataRoots( ...
      v1_1, "/tmp/v1.1", "");
   testCase.verifyEqual([data_root; fixture_root], ...
      ["/tmp/v1.1"; "/tmp/v1.1"]);

   v1_2 = icemodel.test.helpers.formalBaselinePolicy("v1.2");
   canonical = icemodel.verification.setup.fixtureDataRoot("v1.2");
   [data_root, fixture_root] = ...
      icemodel.test.helpers.resolveReleaseDataRoots(v1_2, "", "");
   testCase.verifyEqual([data_root; fixture_root], ...
      [canonical; canonical]);

   [data_root, fixture_root] = ...
      icemodel.test.helpers.resolveReleaseDataRoots( ...
      v1_2, "", "/tmp/v1.2");
   testCase.verifyEqual([data_root; fixture_root], ...
      ["/tmp/v1.2"; "/tmp/v1.2"]);
   testCase.verifyError(@() ...
      icemodel.test.helpers.resolveReleaseDataRoots( ...
      v1_2, "/tmp/model", "/tmp/fixtures"), ...
      'icemodel:test:releaseDataRootMismatch');
   [data_root, fixture_root] = ...
      icemodel.test.helpers.resolveReleaseDataRoots( ...
      v1_2, "/tmp/v1.2/../v1.2", "/tmp/v1.2/");
   testCase.verifyEqual(data_root, "/tmp/v1.2/../v1.2");
   testCase.verifyEqual(fixture_root, "/tmp/v1.2/");
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
   % Formal-suite entrypoints expand "all" once, then run one model at a time.

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
   % The suite bootstrap must install the test demo configuration. Cleanup
   % must restore the caller's configuration.

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

function test_managedBaselineSiblings_returns_every_formal_model(testCase)
   % A build with no explicit output owns one file for each formal model, so
   % worktreeRevision must leave all of them out of the source identity.

   models = icemodel.test.helpers.resolveRequestedSmbmodels("all");
   returned = icemodel.test.helpers.managedBaselineSiblings( ...
      "regression", "rolling", string.empty());
   expected = arrayfun(@(model) string( ...
      icemodel.test.helpers.baselineFilePath("regression", ...
      smbmodel=model)), models);
   testCase.verifyEqual(sort(returned(:)), sort(expected(:)));
end

function test_managedBaselineSiblings_keys_perf_files_by_year(testCase)
   % Perf baseline paths carry the benchmark year; regression paths do not.

   models = icemodel.test.helpers.resolveRequestedSmbmodels("all");
   returned = icemodel.test.helpers.managedBaselineSiblings( ...
      "perf", "v1.2", string.empty(), simyear=2016);
   testCase.verifyNumElements(returned, numel(models));
   testCase.verifyTrue(all(contains(returned, "perf_baseline_2016_v1_2_")));
   testCase.verifyTrue(isstring(returned));
end

function test_managedBaselineSiblings_honors_an_explicit_output(testCase)
   % An explicit output file names the only managed file the build writes.

   returned = icemodel.test.helpers.managedBaselineSiblings( ...
      "perf", "rolling", "/tmp/one-off-baseline.mat");
   testCase.verifyEqual(returned, "/tmp/one-off-baseline.mat");
end

function test_baselineFilePath_resolves_latest_release(testCase)
   % baseline_type="release" without a tag should resolve the latest version.

   pathname = icemodel.test.helpers.baselineFilePath("perf", ...
      baseline_type="release");
   testCase.verifyTrue(contains(pathname, "v1_2"));
end

function test_scalar_snapshot_commands_reload_custom_outputs(testCase)
   % Scalar public snapshot commands validate their saved custom outputs.
   regression_file = string(tempname) + ".mat";
   perf_file = string(tempname) + ".mat";
   cleanup = onCleanup(@() removeSnapshotOutputs( ...
      [regression_file; perf_file]));

   regression = snapshot_regression_baseline( ...
      baseline_tag="v1.2", smbmodel="icemodel", ...
      output_file=regression_file);
   perf = snapshot_perf_baseline( ...
      baseline_tag="v1.2", smbmodel="icemodel", simyear=2016, ...
      output_file=perf_file);

   testCase.verifyTrue(isfile(regression_file));
   testCase.verifyTrue(isfile(perf_file));
   testCase.verifyFalse(isempty(regression));
   testCase.verifyFalse(isempty(perf));
   clear cleanup
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
   %VERIFYCASEFORCINGMATCHES Match each case id and compare its forcing hash.

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

function verifyMatrixSolvers(testCase, cases, expected_icemodel)
   %VERIFYMATRIXSOLVERS Compare the solver ids of each model in a case matrix.
   %
   % icemodel rows carry the expected ids in order. skinmodel rows carry
   % solver 1.

   returned = cases.solver(cases.smbmodel == "icemodel");
   testCase.verifyEqual(returned, expected_icemodel);
   returned = unique(cases.solver(cases.smbmodel == "skinmodel"));
   testCase.verifyEqual(returned, 1);
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
      "function [rootdir, input_path, output_path, eval_path, " ...
      + "cleanup] = bootstrapTestEnvironment(varargin)"
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

function writeBuilderFetchStub(filename)
   %WRITEBUILDERFETCHSTUB Keep ignored release data outside the unit seam.
   lines = [ ...
      "function result = fetchFixtures(varargin)"
      "result = struct();"
      "end"];
   writeTextFixture(filename, join(lines, newline));
end

function writeRunnerBootstrapStub(filename)
   %WRITERUNNERBOOTSTRAPSTUB Stop public runners after setup.

   lines = [ ...
      "function [rootdir, input_path, output_path, eval_path, " ...
      + "cleanup] = bootstrapTestEnvironment(varargin)"
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
   %VERIFYFIXTURECAPABILITYERROR Check the repair path works without network
   % access.

   % Mirror captureExpectedWarning's guard; this also anchors the
   % analyzer's view of the argument the evalc string consumes.
   assert(isa(operation, 'function_handle'))
   try
      % The empty data_root makes config resolution warn about
      % ICEMODEL_INPUT_PATH before the capability check errors. Run the
      % operation under evalc so that expected warning text stays out of
      % the suite log; evalc still propagates the error to the catch.
      evalc('operation();');
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
      'icemodel.test.helpers.resolveRequestedSmbmodels', ...
      'icemodel.verification.setup.fetchFixtures')
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

function meta = aaArtifactMeta(run_name)
   %AAARTIFACTMETA Build the metadata run_aa_acceptance requires.
   meta = struct('ambient_stable', true, 'isolation', "process", ...
      'run_name', string(run_name), 'hostname', "test-host", ...
      'matlab_version', string(version), 'git_revision', "test-rev", ...
      'tier', "smoke", 'simyear', 2016, 'n_runs', 3, 'n_warmups', 1, ...
      'tol_perf', 0.2, 'data_root', "test-data");
end

function case_summary = aaCaseSummary(median_wall_s, valid)
   %AACASESUMMARY Build one complete A/A workload row.
   case_summary = table("icemodel_kanm_2016_solver1", ...
      "promice_filled", median_wall_s, valid, ...
      'VariableNames', {'case_id', 'forcings', 'median_wall_s', 'valid'});
end

function writeAaSuiteStub(stub_file)
   %WRITEAASUITESTUB Write a run_perf_suite stub for the A/A branch test.
   %
   % Each call writes one valid process-isolated artifact beside the stub
   % and returns its path, with a fresh run_name per call.
   stub = [ ...
      "function results = run_perf_suite(varargin)" ...
      "   persistent n" ...
      "   if isempty(n); n = 0; end" ...
      "   n = n + 1;" ...
      "   meta = struct('ambient_stable', true, 'isolation', ""process"", ..." ...
      "      'run_name', ""20260101-00000"" + n, 'hostname', ""test-host"", ..." ...
      "      'matlab_version', string(version), 'git_revision', ""test-rev"", ..." ...
      "      'tier', ""smoke"", 'simyear', 2016, 'n_runs', 3, ..." ...
      "      'n_warmups', 1, 'tol_perf', 0.2, 'data_root', ""test-data"");" ...
      "   case_summary = table(""icemodel_kanm_2016_solver1"", ..." ...
      "      ""promice_filled"", 100.0, true, ..." ...
      "      'VariableNames', {'case_id', 'forcings', ..." ...
      "      'median_wall_s', 'valid'});" ...
      "   file = fullfile(fileparts(mfilename('fullpath')), ..." ...
      "      sprintf('perf_results_stub_%d.mat', n));" ...
      "   save(file, 'meta', 'case_summary');" ...
      "   results = struct('artifact_file', file);" ...
      "end"];
   fid = fopen(stub_file, 'w');
   assert(fid >= 0, 'cannot write the run_perf_suite stub')
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', stub{:});
end

function removeSnapshotOutputs(files)
   %REMOVESNAPSHOTOUTPUTS Remove temporary snapshot files and profiler copies.
   for pathname = files.'
      if isfile(pathname)
         delete(pathname)
      end
      profile_dir = icemodel.test.helpers.baselineProfilerDir(pathname);
      if isfolder(profile_dir)
         rmdir(profile_dir, 's')
      end
   end
end
