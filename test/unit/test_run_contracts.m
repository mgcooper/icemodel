function tests = test_run_contracts
%TEST_RUN_CONTRACTS Verify run-option and path-derivation contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)

   testCase.TestData.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      [2015; 2016; 2017], configure=true, nsteps=24, dt_seconds=3600);
end

function teardown(testCase)

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_setopts_normalizes_none_inputs(testCase)

   opts = icemodel.setopts('skinmodel', 'kanm', [2015 2016 2017], ...
      'kanm', "none", "", "none", false, false);

   testCase.verifyEqual(opts.userdata, 'kanm');
   testCase.verifyEqual(opts.uservars, 'albedo');
   testCase.verifyEqual(opts.testname, '');
   testCase.verifyEqual(opts.output_profile, 'standard');
   testCase.verifyEqual(opts.output_years, [2015 2016 2017]);
end

function test_resetopts_updates_output_years_and_coupler_defaults(testCase)

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

function test_setpath_builds_restart_path_without_blank_parts(testCase)

   output_root = getenv('ICEMODEL_OUTPUT_PATH');
   path_plain = icemodel.setpath('restart', 'kanm', 'skinmodel');
   path_test = icemodel.setpath('restart', 'kanm', 'skinmodel', ...
      '', [], 'case01');

   testCase.verifyEqual(path_plain, ...
      fullfile(output_root, 'kanm', 'skinmodel', 'restart'));
   testCase.verifyEqual(path_test, ...
      fullfile(output_root, 'kanm', 'skinmodel', 'case01', 'restart'));
end

function test_configureRun_builds_default_restart_path(testCase)

   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, testname='restartcase');
   opts = icemodel.resetopts(opts, 'pathrestart', []);
   opts = icemodel.configureRun(opts);

   testCase.verifyEqual(opts.pathrestart, fullfile( ...
      workspace.outputdir, 'kanm', 'skinmodel', 'restartcase', 'restart'));
end
