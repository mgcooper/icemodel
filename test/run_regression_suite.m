function results = run_regression_suite(kwargs)
   %RUN_REGRESSION_SUITE Run formal icemodel numerical regression test suite.
   %
   %  results = run_regression_suite()
   %  results = run_regression_suite(tier="smoke")
   %  results = run_regression_suite(tier="full")
   %  results = run_regression_suite(tier="smoke", smbmodel="skinmodel")
   %  results = run_regression_suite(tier="smoke", smbmodel="icemodel", solver=2)
   %  results = run_regression_suite(tier="smoke", smbmodel="icemodel", solver=[2 3])
   %  results = run_regression_suite(tier="smoke", smbmodel="icemodel", solver=[1 3])
   %  results = run_regression_suite(simyear=2017, smoke_sites="kanm", ...
   %     full_sites=["kanm"; "kanl"])
   %  results = run_regression_suite(tier="full", baseline="v1.1")
   %  results = run_regression_suite(data_root="/path/to/test/data")
   %  results = run_regression_suite(fixture_root="/path/to/provisioned/data")
   %
   % Use this for normal regression comparisons against an existing rolling or
   % release baseline.
   %
   % This function does not update baselines. It runs the formal cases,
   % compares core scalar outputs to the requested baseline, and writes one
   % artifact and one Quarto HTML report under test/artifacts/<run_name>/.
   % RESULTS.failed_cases lists the failed case IDs. RESULTS.failed_gates is a
   % table of each failed case and its failed gates. The saved report stores
   % the per-case gates in its failed_gates column. A unittest failure after
   % the report is saved appears only in RESULTS.failed_gates, as the gate
   % test_framework. A failure before the save raises
   % icemodel:test:regressionArtifactMissing.
   %
   % The optional solver filter accepts any subset of
   % icemodel.namelists.solver().
   % DATA_ROOT overrides the default test case for isolated fixture comparisons.
   % FIXTURE_ROOT selects where a release's required fixture capabilities are
   % verified. A release that runs the model from its own provisioned data
   % (v1.2 onward) rejects a DATA_ROOT that names a different tree with
   % icemodel:test:releaseDataRootMismatch.
   %
   % SMOKE_SITES and FULL_SITES are advanced overrides for the site lists
   % used by each formal tier. Most callers should leave them at the
   % defaults and select only TIER.
   %
   % CLI entrypoint:
   %  matlab -batch "run('/ABS/PATH/icemodel/test/run_regression_suite.m')"

   arguments (Input)

      kwargs.tier (1, :) string ...
         {icemodel.validators.mustBeTestTierName(kwargs.tier)} ...
         = "smoke"

      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} ...
         = "all"

      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} ...
         = []

      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} ...
         = 2016

      kwargs.smoke_sites string ...
         = "kanm"

      kwargs.full_sites string ...
         = ["kanm"; "kanl"]

      kwargs.baseline (1, :) string ...
         = "rolling"

      kwargs.run_name string ...
         = string.empty()

      kwargs.build_report (1, 1) logical ...
         = true

      kwargs.data_root (1, 1) string ...
         = ""

      % See the FIXTURE_ROOT paragraph above for the same-root rule that a
      % release with use_fixture_root_for_model enforces.
      kwargs.fixture_root (1, 1) string ...
         = ""
   end

   % Record this suite in the session activity so a later
   % in-session formal perf run can refuse the contaminated session.
   icemodel.test.helpers.markTestSessionDirty("run_regression_suite");

   % Deal out arguments.
   [tier, smbmodel, solver, simyear, smoke_sites, full_sites, baseline, ...
      run_name, build_report] = deal(kwargs.tier, kwargs.smbmodel, ...
      kwargs.solver, kwargs.simyear, reshape(kwargs.smoke_sites, [], 1), ...
      reshape(kwargs.full_sites, [], 1), kwargs.baseline, ...
      kwargs.run_name, kwargs.build_report);

   % Resolve full path to the test/ dir.
   testdir = icemodel.getpath('test');

   % The baseline registration defines the default data tree and the forcing.
   % An explicit DATA_ROOT still takes precedence inside the bootstrap helper.
   baseline_policy = ...
      icemodel.test.helpers.formalBaselinePolicy(baseline);

   [data_root, fixture_root] = ...
      icemodel.test.helpers.resolveReleaseDataRoots( ...
      baseline_policy, kwargs.data_root, kwargs.fixture_root);

   % Bootstrap the source/test trees once for CLI and interactive runs.
   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, input_path, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment( ...
      icemodel_config_casename=baseline_policy.config_case, ...
      data_root=data_root);

   % Propagate the configured verification root through the TestCase selector.
   % A caller-supplied root keeps precedence over the resolved default.
   if isblanktext(data_root)
      data_root = string(fileparts(input_path));
   end

   % Verify release assets separately from the formal model input tree.
   if ~isempty(baseline_policy.required_fixture_capabilities)
      icemodel.verification.setup.fetchFixtures( ...
         baseline_policy.baseline_tag, ...
         capabilities=baseline_policy.required_fixture_capabilities, ...
         root=fixture_root, download=false);
   end

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);
   if baseline_policy.require_source_revision
      icemodel.test.helpers.assertCommonBaselineRevision( ...
         "regression", baseline, models, simyear);
   end

   % Resolve one run directory before dispatching models so aggregate artifacts
   % and their report always share the same review surface.
   [~, ~, run_name] = icemodel.test.helpers.resolveRunStamp(run_name);

   % Build the unittest suite once, then run the canonical single-model flow
   % for each requested formal model.
   suite = testsuite(fullfile(testdir, 'regression', ...
      'IcemodelRegressionTest.m'));
   runner = matlab.unittest.TestRunner.withTextOutput;

   % Run the canonical single-model workflow for each requested model.
   per_model = arrayfun(@(mdl) runSingleModelRegression( ...
      runner, suite, tier, mdl, solver, simyear, ...
      smoke_sites, full_sites, baseline, run_name, data_root), ...
      models, 'UniformOutput', false);

   % Combine per-model results into a common struct and display.
   results = combineRegressionResults(per_model);

   % Display the results.
   icemodel.test.helpers.displayRegressionResults(results)

   % Render the saved comparison into the common report layer unless the
   % caller explicitly requested an artifact-only run.
   results.report_file = "";
   if build_report
      results.report_file = ...
         icemodel.verification.report.buildTestSuiteReport( ...
         "regression", results);
   end

   % Restore the caller's config last, after every path-dependent step has
   % run against the configured test environment. An early error still
   % restores it, because the object dies with the scope.
   delete(suite_cleanup)
end

function results = runSingleModelRegression(runner, suite, tier, smbmodel, ...
      solver, simyear, smoke_sites, full_sites, baseline, ...
      run_name, data_root)
   %RUNSINGLEMODELREGRESSION Configure one formal model regression run.

   % Reject an unaccepted forcing transition before unittest dispatch so the
   % top-level runner keeps the shared baseline error id.
   formal_baseline = icemodel.test.helpers.loadBaseline("regression", ...
      smbmodel=smbmodel, baseline_tag=baseline);
   icemodel.test.helpers.assertFormalBaselineForcing( ...
      formal_baseline, baseline);

   % Provide the requested regression selection to the unittest class.
   selector_cleanup = configureRegressionSelectorEnv( ...
      tier, smbmodel, solver, simyear, smoke_sites, full_sites, ...
      baseline, run_name, data_root); %#ok<NASGU>

   % Run the formal regression class for this concrete smbmodel. Clear the
   % artifact path first, so a class that fails before saving cannot return
   % the artifact of an earlier model run in this session.
   setenv('ICEMODEL_REGRESSION_ARTIFACT_FILE', '');
   test_result = runner.run(suite);

   % Load the artifact saved by IcemodelRegressionTest to build the
   % results struct that mirrors the perf suite contract.
   artifact_file = string(getenv('ICEMODEL_REGRESSION_ARTIFACT_FILE'));
   if artifact_file == "" || ~isfile(artifact_file)
      error('icemodel:test:regressionArtifactMissing', ...
         ['The regression test did not write its comparison artifact. ' ...
         'Inspect the unittest diagnostics above for the setup or test failure.'])
   end
   S = load(artifact_file, 'report', 'case_opts', 'meta');

   results = struct();
   results.report = S.report;
   results.case_opts = S.case_opts;
   results.meta = S.meta;
   results.artifact_file = artifact_file;
   results.test_result = test_result;
   results.passed = all([test_result.Passed]);

   % Name every failed case and its failed gates from the saved report rows.
   [results.failed_cases, results.failed_gates] = ...
      icemodel.test.helpers.regressionFailures(S.report, results.passed);
end

function results = combineRegressionResults(per_model)
   %COMBINEREGRESSIONRESULTS Merge one-or-more single-model regression results.

   if isscalar(per_model)
      results = per_model{1};
      return
   end

   % Extract each returned field once, then concatenate the per-model pieces.
   report = cellfun(@(s) s.report, per_model, 'UniformOutput', false);
   case_opts = cellfun(@(s) s.case_opts(:), per_model, 'UniformOutput', false);
   meta = cellfun(@(s) s.meta, per_model, 'UniformOutput', false);
   artifact_file = cellfun(@(s) string(s.artifact_file(:)), per_model, ...
      'UniformOutput', false);
   test_result = cellfun(@(s) s.test_result, per_model, 'UniformOutput', false);
   failed_cases = cellfun(@(s) string(s.failed_cases(:)), per_model, ...
      'UniformOutput', false);
   failed_gates = cellfun(@(s) s.failed_gates, per_model, ...
      'UniformOutput', false);
   pass_flags = cellfun(@(s) s.passed, per_model);

   results = struct();
   results.report = vertcat(report{:});
   results.case_opts = vertcat(case_opts{:});
   results.meta = vertcat(meta{:});
   results.artifact_file = vertcat(artifact_file{:});
   results.test_result = horzcat(test_result{:});
   results.failed_cases = vertcat(failed_cases{:});
   results.failed_gates = vertcat(failed_gates{:});
   results.passed = all(pass_flags);
end

function cleanup = configureRegressionSelectorEnv(tier, smbmodel, solver, ...
      simyear, smoke_sites, full_sites, baseline, run_name, data_root)
   %CONFIGUREREGRESSIONSELECTORENV Export one regression selection contract.
   %
   % This is file-local glue for the regression runner. The unittest class
   % reads one suite-selection contract from env, unlike the perf flow which
   % installs per-case env inside runPerfCase.

   % Snapshot the caller's current selector state before overwriting it.
   names = [ ...
      "ICEMODEL_TEST_TIER"
      "ICEMODEL_TEST_SMBMODEL_FILTER"
      "ICEMODEL_TEST_SOLVER_FILTER"
      "ICEMODEL_TEST_SIMYEAR_FILTER"
      "ICEMODEL_TEST_SMOKE_SITES"
      "ICEMODEL_TEST_FULL_SITES"
      "ICEMODEL_REGRESSION_BASELINE"
      "ICEMODEL_TEST_RUN_NAME"
      "ICEMODEL_TEST_DATA_ROOT"];
   oldvals = arrayfun(@(name) string(getenv(name)), names, ...
      'UniformOutput', false);

   % Export the requested regression selection for one runner invocation.
   setenv('ICEMODEL_TEST_TIER', char(tier));
   setenv('ICEMODEL_TEST_SMBMODEL_FILTER', char(smbmodel));
   setenv('ICEMODEL_TEST_SOLVER_FILTER', char(join(string(solver), ',')));
   setenv('ICEMODEL_TEST_SIMYEAR_FILTER', int2str(simyear));
   setenv('ICEMODEL_TEST_SMOKE_SITES', char(join(smoke_sites, ',')));
   setenv('ICEMODEL_TEST_FULL_SITES', char(join(full_sites, ',')));
   setenv('ICEMODEL_REGRESSION_BASELINE', char(baseline));
   setenv('ICEMODEL_TEST_RUN_NAME', char(run_name));
   setenv('ICEMODEL_TEST_DATA_ROOT', char(data_root));

   % Restore the prior selector state when the runner returns.
   cleanup = onCleanup(@() restoreSelectorEnv(names, oldvals));
end

function restoreSelectorEnv(names, values)
   %RESTORESELECTORENV Restore prior regression selector env values.

   for n = 1:numel(names)
      setenv(names(n), values{n});
   end
end
