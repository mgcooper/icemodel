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
   %
   % Use this for normal regression comparisons against an existing rolling or
   % release baseline.
   %
   % This function does not update baselines; it only runs the formal cases,
   % compares core scalar outputs to the requested baseline, and writes one
   % artifact under test/artifacts/<run_name>/.
   %
   % The optional solver filter accepts any subset of [1 2 3].
   %
   % SMOKE_SITES and FULL_SITES are advanced overrides for the site lists
   % used by each formal tier. Most callers should leave them at the
   % defaults and select only TIER.
   %
   % CLI entrypoint:
   %  matlab -batch "run('/ABS/PATH/icemodel/test/run_regression_suite.m')"

   arguments (Input)
      kwargs.tier (1, :) string ...
         {icemodel.validators.mustBeTestTierName(kwargs.tier)} = "smoke"
      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} = "all"
      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} = []
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.smoke_sites string = "kanm"
      kwargs.full_sites string = ["kanm"; "kanl"]
      kwargs.baseline (1, :) string = "rolling"
      kwargs.run_name string = string.empty()
   end

   % Deal out arguments.
   [tier, smbmodel, solver, simyear, smoke_sites, full_sites, baseline, ...
      run_name] = deal(kwargs.tier, kwargs.smbmodel, kwargs.solver, ...
      kwargs.simyear, reshape(kwargs.smoke_sites, [], 1), ...
      reshape(kwargs.full_sites, [], 1), kwargs.baseline, kwargs.run_name);

   % Resolve full path to the test/ dir.
   thisdir = icemodel.getpath('test');

   % Bootstrap the source/test trees once for CLI and interactive runs.
   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);

   % Build the unittest suite once, then run the canonical single-model flow
   % for each requested formal model.
   suite = testsuite(fullfile(thisdir, 'regression', ...
      'IcemodelRegressionTest.m'));
   runner = matlab.unittest.TestRunner.withTextOutput;

   % Run the canonical single-model workflow for each requested model.
   per_model = arrayfun(@(mdl) runSingleModelRegression( ...
      runner, suite, tier, mdl, solver, simyear, ...
      smoke_sites, full_sites, baseline, run_name), ...
      models, 'UniformOutput', false);

   % Combine results into a common struct.
   results = vertcat(per_model{:});
end

function results = runSingleModelRegression(runner, suite, tier, smbmodel, ...
      solver, simyear, smoke_sites, full_sites, baseline, run_name)
   %RUNSINGLEMODELREGRESSION Configure one formal model regression run.

   % Provide the requested regression selection to the unittest class.
   selector_cleanup = configureRegressionSelectorEnv( ...
      tier, smbmodel, solver, simyear, smoke_sites, full_sites, ...
      baseline, run_name); %#ok<NASGU>

   % Run the formal regression class for this concrete smbmodel.
   results = runner.run(suite);
end

function cleanup = configureRegressionSelectorEnv(tier, smbmodel, solver, ...
      simyear, smoke_sites, full_sites, baseline, run_name)
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
      "ICEMODEL_TEST_RUN_NAME"];
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

   % Restore the prior selector state when the runner returns.
   cleanup = onCleanup(@() restoreSelectorEnv(names, oldvals));
end

function restoreSelectorEnv(names, values)
   %RESTORESELECTORENV Restore prior regression selector env values.

   for n = 1:numel(names)
      setenv(names(n), values{n});
   end
end
