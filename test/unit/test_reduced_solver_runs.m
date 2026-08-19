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

function test_skinmodel_forced_advance_restores_checkpoint(testCase)
   % A forced SkinModel advance may consume time only. The final state and
   % diagnosed conductive flux must therefore come from the last checkpoint.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=900);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'skinmodel', solver=1, testname='skin_forced_advance');
   f_liq_checkpoint = 0.01 * ones(size(state.f_liq));
   restart = struct('T', state.T, 'f_ice', state.f_ice, ...
      'f_liq', f_liq_checkpoint, 'Ts', state.Ts, 'r_eff', state.r_eff);
   restart_file = fullfile(workspace.rootdir, ...
      'skin-forced-advance-restart.mat');
   save(restart_file, 'restart');

   opts = state.opts;
   opts.output_profile = 'diagnostic';
   opts.vars1 = {};
   opts.vars2 = {};
   opts.maxiter = 2;
   opts.tol = 0;
   opts.dt = 1;
   opts.use_restart = true;
   opts.restartfile = restart_file;

   [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts);

   % Zero tolerance rejects the nonlinear solve. The one-second full step
   % reaches the bounded fallback immediately and must not accept trial state.
   testCase.verifyEqual(ice1.dt_sum, 1, 'AbsTol', 0);
   testCase.verifyEqual(ice1.n_subfail, 1, 'AbsTol', 0);
   testCase.verifyEqual(ice1.Tice_converged, 0, 'AbsTol', 0);
   testCase.verifyEqual(ice1.Tsfc, state.Ts, 'AbsTol', 0);
   testCase.verifyEqual(ice2.Tice, state.T, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_ice, state.f_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_liq, f_liq_checkpoint, 'AbsTol', 0);

   % Qc observes both restored temperature and restored conductivity, so it
   % catches a stale rejected-solve k_eff even when the saved state is correct.
   k_eff_checkpoint = icemodel.column.bulk_thermal_conductivity( ...
      state.T, state.f_ice, f_liq_checkpoint, 0);
   Qc_checkpoint = icemodel.surface.conductive_heat_flux( ...
      k_eff_checkpoint, state.T, state.dz, state.Ts);
   testCase.verifyEqual(ice1.Qc, Qc_checkpoint, 'AbsTol', 1e-12);

   clear cleanup
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
