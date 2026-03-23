function tests = test_run_contracts
   %TEST_RUN_CONTRACTS Verify run-option and path-derivation contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build one synthetic workspace that exposes path, spinup, and output-
   % year contracts without depending on external project data.

   testCase.TestData.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      [2015; 2016; 2017], configure=true, nsteps=24, dt_seconds=3600);
end

function teardown(testCase)
   % Remove the shared synthetic workspace after the contract checks finish.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_setopts_normalizes_none_inputs(testCase)
   % Legacy "none" placeholders should be normalized into the current
   % userdata/uservars/testname contract.

   opts = icemodel.setopts('skinmodel', 'kanm', [2015 2016 2017], ...
      'kanm', "none", "", "none", false, false);

   testCase.verifyEqual(opts.userdata, 'kanm');
   testCase.verifyEqual(opts.uservars, 'albedo');
   testCase.verifyEqual(opts.testname, '');
   testCase.verifyEqual(opts.output_profile, 'standard');
   testCase.verifyEqual(opts.output_years, [2015 2016 2017]);
end

function test_resetopts_updates_output_years_and_coupler_defaults(testCase)
   % RESETOPTS should recompute retained output years and switch coupler
   % defaults when the solver mode changes.

   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', [2015 2016 2017], solver=3);

   opts = icemodel.resetopts(opts, 'n_spinup_years', 1, 'solver', 2);
   testCase.verifyEqual(opts.output_years, [2016 2017]);
   testCase.verifyEqual(opts.cpl_maxiter, 1);

   opts = icemodel.resetopts(opts, 'solver', 3);
   testCase.verifyEqual(opts.cpl_maxiter, 100);
end

function test_configureRun_preserves_explicit_overrides(testCase)
   % CONFIGURERUN should honor explicit caller overrides instead of
   % rebuilding those fields from the default path/case contract.

   workspace = testCase.TestData.workspace;
   custom_output = fullfile(workspace.rootdir, 'custom_output');
   custom_restart = fullfile(workspace.rootdir, 'custom_restart');
   custom_met = {fullfile(workspace.metdir, 'custom_met.mat')};

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   opts = icemodel.resetopts(opts, ...
      'pathoutput', custom_output, ...
      'pathrestart', custom_restart, ...
      'casename', 'custom_case', ...
      'metfname', custom_met, ...
      'vars1', {'Tsfc'}, ...
      'vars2', {'Tice'});
   opts = icemodel.configureRun(opts);

   testCase.verifyEqual(opts.pathoutput, custom_output);
   testCase.verifyEqual(opts.pathrestart, custom_restart);
   testCase.verifyEqual(opts.casename, 'custom_case');
   testCase.verifyEqual(opts.metfname, custom_met);
   testCase.verifyEqual(opts.vars1, {'Tsfc'});
   testCase.verifyEqual(opts.vars2, {'Tice'});
end

function test_getpath_builds_restart_path_without_blank_parts(testCase)
   % GETPATH should omit blank components instead of leaving empty folders
   % in the restart path.

   output_root = icemodel.getpath('output');
   path_plain = icemodel.getpath('restart', 'kanm', 'skinmodel');
   path_test = icemodel.getpath('restart', 'kanm', 'skinmodel', ...
      '', [], 'case01');

   testCase.verifyEqual(path_plain, ...
      fullfile(output_root, 'kanm', 'skinmodel', 'restart'));
   testCase.verifyEqual(path_test, ...
      fullfile(output_root, 'kanm', 'skinmodel', 'case01', 'restart'));
end

function test_setpath_remains_a_compatibility_alias(testCase)
   % SETPATH should continue to match GETPATH while older callers migrate.

   returned = icemodel.setpath('restart', 'kanm', 'skinmodel', '', [], ...
      'case01');
   expected = icemodel.getpath('restart', 'kanm', 'skinmodel', '', [], ...
      'case01');

   testCase.verifyEqual(returned, expected);
end

function test_configureRun_builds_default_restart_path(testCase)
   % When pathrestart is cleared, CONFIGURERUN should rebuild it from the
   % configured output root and case identity.

   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, testname='restartcase');
   opts = icemodel.resetopts(opts, 'pathrestart', []);
   opts = icemodel.configureRun(opts);

   testCase.verifyEqual(opts.pathrestart, fullfile( ...
      workspace.outputdir, 'kanm', 'skinmodel', 'restartcase', 'restart'));
end
