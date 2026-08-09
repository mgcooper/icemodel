function tests = test_reduced_solver_runs
   %TEST_REDUCED_SOLVER_RUNS Verify bounded output on small synthetic runs.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Reuse one reduced synthetic workspace so each solver-mode test runs on
   % the same controlled forcing and column geometry.

   testCase.TestData.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      2016, configure=true, nsteps=96, dt_seconds=900);
end

function teardown(testCase)
   % Remove the shared reduced workspace after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_skinmodel_reduced_run_stays_bounded(testCase)
   % A compact skinmodel run should finish with bounded, finite outputs on
   % the shared synthetic forcing file.

   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, solver=1, testname='skin_unit');

   [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
   [ice1_pp, ice2_pp] = icemodel.postprocess( ...
      ice1_raw, ice2_raw, opts, opts.output_years);

   icemodel.test.verify.verifyProcessedOutputBounds( ...
      testCase, ice1_pp, ice2_pp);
   testCase.verifyEqual(height(ice1_pp), workspace.nsteps / 4);
end

function test_shared_model_case_path_returns_postprocessed_output(testCase)
   % The shared formal-case path should resolve, run, and postprocess one case.

   workspace = testCase.TestData.workspace;
   c = struct( ...
      'smbmodel', "icemodel", ...
      'sitename', string(workspace.sitename), ...
      'forcings', string(workspace.forcings), ...
      'userdata', "", ...
      'uservars', "", ...
      'simyears', workspace.simyears, ...
      'n_spinup_years', 0, ...
      'solver', 0);

   [ice1, ice2, opts] = icemodel.test.helpers.runModelCase( ...
      c, output_profile="diagnostic");

   % Diagnostic profiles contain intentionally unavailable scheme-specific
   % fields, so verify the shared physical and mass-ledger outputs directly.
   testCase.verifyTrue(all(isfinite(ice1.tsfc)));
   testCase.verifyTrue(all(isfinite(ice2.Tice), 'all'));
   testCase.verifyTrue(all(isfinite(ice1.mass_budget_solid_start_mwe)));
   testCase.verifyTrue(all(isfinite(ice1.mass_budget_solid_end_mwe)));
   testCase.verifyEqual(height(ice1), workspace.nsteps / 4);
   testCase.verifyEqual(opts.output_years, workspace.simyears);
   testCase.verifyEqual(string(opts.output_profile), "diagnostic");
end

function test_icemodel_reduced_runs_stay_bounded_across_solver_modes(testCase)
   % Run the icemodel across all supported solver modes to make sure the
   % reduced synthetic case remains finite and bounded in each branch.
   % Modes: 0 = Dirichlet single-sweep, 1 = Dirichlet iterative,
   %        2 = Robin single-sweep, 3 = Robin iterative.

   workspace = testCase.TestData.workspace;
   for solver = [0, 1, 2, 3]
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
