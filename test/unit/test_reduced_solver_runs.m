function tests = test_reduced_solver_runs
%TEST_REDUCED_SOLVER_RUNS Verify bounded output on small synthetic runs.
   tests = functiontests(localfunctions);
end

function setup(testCase)

   testCase.TestData.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      2016, configure=true, nsteps=24, dt_seconds=3600);
end

function teardown(testCase)

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_skinmodel_reduced_run_stays_bounded(testCase)

   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, solver=1, testname='skin_unit');

   [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
   [ice1_pp, ice2_pp] = icemodel.postprocess( ...
      ice1_raw, ice2_raw, opts, opts.output_years);

   icemodel.test.verify.verifyProcessedOutputBounds( ...
      testCase, ice1_pp, ice2_pp);
   testCase.verifyEqual(height(ice1_pp), workspace.nsteps);
end

function test_icemodel_reduced_runs_stay_bounded_across_solver_modes(testCase)

   workspace = testCase.TestData.workspace;
   for solver = 1:3
      opts = icemodel.test.helpers.buildSyntheticOpts( ...
         workspace, 'icemodel', 2016, ...
         solver=solver, testname=['ice_unit_s' int2str(solver)]);

      [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
      [ice1_pp, ice2_pp] = icemodel.postprocess( ...
         ice1_raw, ice2_raw, opts, opts.output_years);

      icemodel.test.verify.verifyProcessedOutputBounds( ...
         testCase, ice1_pp, ice2_pp);
      testCase.verifyTrue(all(isfinite(ice1_pp.Tice_numiter)));
      testCase.verifyGreaterThanOrEqual(min(ice1_pp.Tice_numiter), 0);
   end
end
