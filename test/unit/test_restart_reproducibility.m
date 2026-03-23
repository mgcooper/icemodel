function tests = test_restart_reproducibility
   %TEST_RESTART_REPRODUCIBILITY Verify year-boundary restarts reproduce output.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Keep restart workspaces under the canonical test-artifacts tmp root so
   % failed cases leave inspectable artifacts in the repo-local test tree.
   tmpdir = fullfile(icemodel.getpath('test'), 'artifacts', 'tmp');
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace([2015; 2016], ...
      configure=false, parentdir=tmpdir);

   testCase.TestData.rootdir = workspace.rootdir;
   testCase.TestData.inputdir = workspace.inputdir;
   testCase.TestData.evaldir = workspace.evaldir;
   testCase.TestData.outdir = workspace.outputdir;
end

function teardown(testCase)
   % Remove the restart workspace after the file-level checks complete.

   rootdir = testCase.TestData.rootdir;
   if exist(rootdir, 'dir') == 7
      rmdir(rootdir, 's');
   end
end

function test_skinmodel_year_boundary_restart(testCase)
   % Skinmodel should reproduce the retained output when resumed from the
   % year-boundary restart file.

   verifyRestartRun(testCase, "skinmodel", 1);
end

function test_icemodel_year_boundary_restart(testCase)
   % Icemodel should reproduce the retained output when resumed from the
   % year-boundary restart file.

   verifyRestartRun(testCase, "icemodel", 2);
end

function verifyRestartRun(testCase, smbmodel, solver)
   % Compare a continuous two-year run against a split seed/resume pair and
   % require exact agreement before and after postprocessing.

   simyears = [2015 2016];

   opts_full = buildOpts(testCase, smbmodel, simyears, ...
      run_tag="full", solver=solver, n_spinup_years=1);
   [ice1_full, ice2_full, opts_full] = icemodel.test.helpers.runSmbModel(opts_full);

   opts_seed = buildOpts(testCase, smbmodel, 2015, ...
      run_tag="seed", solver=solver, saverestart=true);
   [~, ~, opts_seed] = icemodel.test.helpers.runSmbModel(opts_seed);
   restart_file = icemodel.restartfile(opts_seed, 2015);
   testCase.verifyTrue(exist(restart_file, 'file') == 2, ...
      'expected a saved restart file after the seed run');
   restart = load(restart_file, 'restart');
   testCase.verifyTrue(isfield(restart.restart, 'opts'), ...
      'expected restart file to preserve resolved opts metadata');

   opts_resume = buildOpts(testCase, smbmodel, 2016, ...
      run_tag="resume", solver=solver, use_restart=true, restartfile=restart_file);
   [ice1_restart, ice2_restart, opts_restart] = icemodel.test.helpers.runSmbModel(opts_resume);

   [ice1_full_pp, ice2_full_pp] = icemodel.postprocess( ...
      ice1_full, ice2_full, opts_full, opts_full.output_years);
   [ice1_restart_pp, ice2_restart_pp] = icemodel.postprocess( ...
      ice1_restart, ice2_restart, opts_restart, opts_restart.output_years);

   icemodel.test.verify.verifyEqualNested(testCase, ice1_full, ice1_restart, 1e-12);
   icemodel.test.verify.verifyEqualNested(testCase, ice2_full, ice2_restart, 1e-12);
   icemodel.test.verify.verifyEqualNested(testCase, ice1_full_pp, ice1_restart_pp, 1e-12);
   icemodel.test.verify.verifyEqualNested(testCase, ice2_full_pp, ice2_restart_pp, 1e-12);
end

function opts = buildOpts(testCase, smbmodel, simyears, kwargs)
   %BUILDOPTS Build one restart-test OPTS struct inside the synthetic workspace.

   arguments
      testCase
      smbmodel (1, 1) string
      simyears
      kwargs.run_tag (1, 1) string
      kwargs.solver (1, 1) double
      kwargs.n_spinup_years (1, 1) double = 0
      kwargs.saverestart (1, 1) logical = false
      kwargs.use_restart (1, 1) logical = false
      kwargs.restartfile (1, :) char = ''
   end

   opts = icemodel.setopts(char(smbmodel), 'kanm', simyears, ...
      'kanm', 'kanm', 'albedo', char(kwargs.run_tag), false, false);

   % Override the workspace roots so the restart run stays isolated inside
   % the synthetic test workspace.
   pathoutput = fullfile(testCase.TestData.outdir, char(kwargs.run_tag));
   pathuserdata = fullfile(testCase.TestData.inputdir, 'userdata');
   if exist(pathuserdata, 'dir') ~= 7
      mkdir(pathuserdata);
   end

   opts = icemodel.resetopts(opts, ...
      'dt', 3600, ... % keep the synthetic restart test small and deterministic
      'solver', kwargs.solver, ...
      'n_spinup_years', kwargs.n_spinup_years, ...
      'pathinput', testCase.TestData.inputdir, ...
      'pathuserdata', pathuserdata, ...
      'patheval', testCase.TestData.evaldir, ...
      'pathoutput', pathoutput, ...
      'output_profile', 'standard', ...
      'saverestart', kwargs.saverestart, ...
      'use_restart', kwargs.use_restart, ...
      'restartfile', kwargs.restartfile);

   % Finalize the derived output/restart paths after overriding the roots.
   opts = icemodel.configureRun(opts);
end
