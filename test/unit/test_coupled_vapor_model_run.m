function tests = test_coupled_vapor_model_run
   %TEST_COUPLED_VAPOR_MODEL_RUN Verify production vapor coupling end to end.
   %
   % The kernel tests isolate face transport and constrained phase transfer.
   % These tests run the model. They verify that the retired flag cannot
   % select alternate physics. They also verify surface-vapor closure,
   % per-phase storage closure, and nonzero interior redistribution.
   %
   % See also: test_coupled_vapor_transport,
   %  icemodel.column.couple_vapor_step
   tests = functiontests(localfunctions);
end

function test_option_contract_has_no_vapor_coupling_flag(testCase)
   % Coupled vapor transport is unconditional and has no selectable flag.

   opts = icemodel.setopts("icemodel", "kanm", 2016, "kanm");
   testCase.verifyFalse(isfield(opts, 'use_coupled_vapor'));
end

function test_stale_flag_field_cannot_select_alternate_physics(testCase)
   % Old saved option structs may carry the retired field. Ordinary struct
   % handling preserves it, but neither value may select different physics.

   base_opts = syntheticRunOpts(testCase);
   [ice1_default, ice2_default] = ...
      icemodel.test.helpers.runSmbModel(base_opts);

   stale_false = base_opts;
   stale_false.use_coupled_vapor = false;
   [ice1_false, ice2_false] = ...
      icemodel.test.helpers.runSmbModel(stale_false);

   stale_true = base_opts;
   stale_true.use_coupled_vapor = true;
   [ice1_true, ice2_true] = ...
      icemodel.test.helpers.runSmbModel(stale_true);

   testCase.verifyEqual(ice1_false, ice1_default);
   testCase.verifyEqual(ice2_false, ice2_default);
   testCase.verifyEqual(ice1_true, ice1_default);
   testCase.verifyEqual(ice2_true, ice2_default);
end

function test_production_run_closes_the_vapor_identity(testCase)
   % Potential surface energy equals realized phase storage, condensation
   % overflow, and unapplied energy. Interior redistribution stays separate.

   base_opts = syntheticRunOpts(testCase);
   [Ls, Lv, ro_liq] = icemodel.physicalConstant('Ls', 'Lv', 'ro_liq');
   [ice1, ~] = icemodel.test.helpers.runSmbModel(base_opts);

   accounted = ro_liq * ( ...
      Ls * ice1.mass_budget_vapor_solid_mwe ...
      + Lv * ice1.mass_budget_vapor_liquid_mwe ...
      + Lv * ice1.mass_budget_condensation_overflow_mwe) ...
      + ice1.mass_budget_unapplied_vapor_j_m2;
   potential = ice1.mass_budget_vapor_potential_j_m2;

   % Allow only column-integration subtraction roundoff.
   testCase.verifyEqual(potential, accounted, 'AbsTol', 1e-4);
   scale = max(abs(potential));
   testCase.verifyGreaterThan(scale, 0);
   testCase.verifyLessThan(max(abs(potential - accounted)) / scale, 1e-8);
end

function test_production_run_closes_the_per_phase_storage(testCase)
   % Endpoint phase storage must equal phase change, surface vapor exchange,
   % remeshing, and interior vapor redistribution together.

   base_opts = syntheticRunOpts(testCase);
   [ice1, ~] = icemodel.test.helpers.runSmbModel(base_opts);

   solid_delta = ice1.mass_budget_solid_end_mwe ...
      - ice1.mass_budget_solid_start_mwe;
   liquid_delta = ice1.mass_budget_liquid_end_mwe ...
      - ice1.mass_budget_liquid_start_mwe;

   testCase.verifyEqual(solid_delta, ...
      ice1.mass_budget_phase_solid_mwe ...
      + ice1.mass_budget_vapor_solid_mwe ...
      + ice1.mass_budget_remesh_solid_mwe ...
      + ice1.mass_budget_vapor_redistribution_solid_mwe, 'AbsTol', 1e-10);
   testCase.verifyEqual(liquid_delta, ...
      ice1.mass_budget_phase_liquid_mwe ...
      + ice1.mass_budget_vapor_liquid_mwe ...
      + ice1.mass_budget_remesh_liquid_mwe ...
      + ice1.mass_budget_vapor_redistribution_liquid_mwe, 'AbsTol', 1e-10);

   % Gross storage movement must bound the signed endpoint change.
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_solid_storage_gross_mwe + 1e-12, abs(solid_delta));
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_liquid_storage_gross_mwe + 1e-12, ...
      abs(liquid_delta));
end

function test_production_run_transports_interior_vapor(testCase)
   % A production run must exercise nonzero interior redistribution.

   base_opts = syntheticRunOpts(testCase);
   [ice1, ~] = icemodel.test.helpers.runSmbModel(base_opts);

   gross = ice1.mass_budget_vapor_redistribution_solid_gross_mwe ...
      + ice1.mass_budget_vapor_redistribution_liquid_gross_mwe;
   testCase.verifyGreaterThan(max(gross), 0);
end

function test_robin_outer_failure_recovers_from_checkpoint(testCase)
   % A healthy-inner Robin outer failure must restart from the prognostic
   % checkpoint. Its accepted result must equal a direct conservative run.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=2, dt_seconds=900);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='robin_recovery');

   % Alpha 2 without acceleration deterministically leaves the synthetic
   % inner solve healthy while exhausting the strong Robin outer loop.
   [ok_seb, ok_ieb, ok_cpl] = runRobinProbe(state, 2.0, false);

   testCase.verifyTrue(ok_seb);
   testCase.verifyTrue(ok_ieb);
   testCase.verifyFalse(ok_cpl);
   % The production retry must discard every failed-attempt output. Its result
   % is therefore exactly the same as starting with the conservative pair. The
   % first forcing step records the successful latch and the second proves that
   % the diagnostic resets while the latched policy remains active.
   primary_opts = icemodel.resetopts(state.opts, ...
      'cpl_alpha', 2.0, 'cpl_aitken', false, ...
      'output_profile', 'diagnostic');
   recovery_alpha = icemodel.parameterLookup('cpl_recovery_alpha');
   conservative_opts = icemodel.resetopts(state.opts, ...
      'cpl_alpha', recovery_alpha, 'cpl_aitken', false, ...
      'output_profile', 'diagnostic');
   [ice1_recovered, ice2_recovered] = ...
      icemodel.test.helpers.runSmbModel(primary_opts);
   [ice1_conservative, ice2_conservative] = ...
      icemodel.test.helpers.runSmbModel(conservative_opts);

   testCase.verifyEqual(ice1_recovered.cpl_recovery_count(:), [1; 0]);
   testCase.verifyEqual(ice1_conservative.cpl_recovery_count(:), [0; 0]);
   ice1_recovered = rmfield(ice1_recovered, 'cpl_recovery_count');
   ice1_conservative = rmfield(ice1_conservative, 'cpl_recovery_count');
   testCase.verifyEqual(ice1_recovered, ice1_conservative);
   testCase.verifyEqual(ice2_recovered, ice2_conservative);
   testCase.verifyEqual(ice1_recovered.n_subfail(:), zeros(2, 1));
   testCase.verifyEqual(ice1_recovered.Tsfc_converged(:), ones(2, 1));
   testCase.verifyEqual(ice1_recovered.Tice_converged(:), ones(2, 1));

end

function test_robin_recovery_preserves_forced_advance_boundary(testCase)
   % The last allowed failed attempt belongs to CHECKSUBSTEP. A conservative
   % retry must not replace its checkpoint-only forced advance.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=900);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='robin_recovery_boundary');

   % DT=1 makes MAXSUBSTEP=1. Alpha 2 leaves the inner solve healthy while the
   % primary Robin outer loop fails, so that failure must force immediately.
   opts = state.opts;
   opts.dt = 1;
   opts.cpl_alpha = 2.0;
   opts.cpl_aitken = false;
   opts.output_profile = 'diagnostic';
   opts.vars1 = {};
   opts.vars2 = {};
   [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts);

   testCase.verifyEqual(ice1.cpl_recovery_count, 0);
   testCase.verifyEqual(ice1.n_subfail, 1);
   testCase.verifyEqual(ice1.dt_sum, 1);
   testCase.verifyEqual(ice1.Tsfc, state.Ts, 'AbsTol', 0);
   testCase.verifyEqual(ice2.Tice, state.T, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_ice, state.f_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_liq, state.f_liq, 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_liq, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_evp, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_lyr, zeros(size(state.f_liq)), 'AbsTol', 0);

   % No physical ledger producer runs after the checkpoint-only advance.
   budget_fields = string(fieldnames(ice1));
   process_fields = budget_fields(startsWith(budget_fields, 'mass_budget_') ...
      & ~endsWith(budget_fields, {'_start_mwe', '_end_mwe'}));
   for k = 1:numel(process_fields)
      testCase.verifyEqual(ice1.(process_fields(k)), 0, 'AbsTol', 0, ...
         sprintf('%s changed during forced acceptance', process_fields(k)));
   end
end

function test_failed_robin_recovery_falls_through_substep_control(testCase)
   % A conservative retry can fail too. It must leave the checkpoint and latch
   % untouched, then use the ordinary reduced-step and forced-advance path.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=900);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='robin_recovery_failure');

   % Two strict outer sweeps keep both policies nonconverged while their inner
   % solves remain healthy. DT=2 permits one retry before the next failure
   % reaches the forced-advance boundary.
   probe = state;
   probe.opts.dt = 2;
   probe.cpl_maxiter = 2;
   probe.cpl_Ts_tol = 0;
   probe.cpl_seb_tol = 0;
   recovery_alpha = icemodel.parameterLookup('cpl_recovery_alpha');
   [primary_seb, primary_ieb, primary_cpl] = ...
      runRobinProbe(probe, 2.0, false);
   [retry_seb, retry_ieb, retry_cpl] = ...
      runRobinProbe(probe, recovery_alpha, false);

   testCase.verifyTrue(primary_seb && primary_ieb && ~primary_cpl);
   testCase.verifyTrue(retry_seb && retry_ieb && ~retry_cpl);

   opts = state.opts;
   opts.dt = 2;
   opts.cpl_maxiter = 2;
   opts.cpl_Ts_tol = 0;
   opts.cpl_seb_tol = 0;
   opts.cpl_alpha = 2.0;
   opts.cpl_aitken = false;
   opts.output_profile = 'diagnostic';
   opts.saverestart = true;
   opts.vars1 = {};
   opts.vars2 = {};

   % Four coupler calls distinguish primary+retry followed by two one-second
   % primaries from the three primary calls run when the retry is skipped.
   prior = profile('status');
   restore_profiler = onCleanup(@() restoreProfiler(prior));
   profile off
   profile clear
   profile on
   [ice1, ice2, returned_opts] = ...
      icemodel.test.helpers.runSmbModel(opts);
   profile off
   profile_data = profile('info');
   profile clear
   clear restore_profiler

   profile_names = string({profile_data.FunctionTable.FunctionName});
   coupler_rows = endsWith(profile_names, 'solve_surface_column_robin');
   coupler_calls = sum([profile_data.FunctionTable(coupler_rows).NumCalls]);
   testCase.verifyEqual(coupler_calls, 4);

   testCase.verifyEqual(ice1.cpl_recovery_count, 0);
   testCase.verifyEqual(ice1.n_subfail, 2);
   testCase.verifyEqual(ice1.dt_sum, 2);
   testCase.verifyEqual(ice1.Tsfc, state.Ts, 'AbsTol', 0);
   testCase.verifyEqual(ice2.Tice, state.T, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_ice, state.f_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_liq, state.f_liq, 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_liq, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_evp, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_lyr, zeros(size(state.f_liq)), 'AbsTol', 0);

   restart_file = icemodel.restartfile(returned_opts, 2016);
   saved = load(restart_file, 'restart');
   testCase.verifyFalse(saved.restart.use_conservative_cpl);

   % A failed retry and final forced advance must not run physical producers.
   budget_fields = string(fieldnames(ice1));
   process_fields = budget_fields(startsWith(budget_fields, 'mass_budget_') ...
      & ~endsWith(budget_fields, {'_start_mwe', '_end_mwe'}));
   for k = 1:numel(process_fields)
      testCase.verifyEqual(ice1.(process_fields(k)), 0, 'AbsTol', 0, ...
         sprintf('%s changed during failed recovery', process_fields(k)));
   end
end

function test_forced_advance_is_time_only(testCase)
   % A rejected solve may advance elapsed time, but it must not mutate the
   % restored checkpoint or record any physical process increment.
   % Bounded fixture probes produce only wholly accepted or wholly forced
   % steps, so no deterministic driver test covers a mixed accepted duration
   % without an artificial solver seam. The exact wet-duration dependence is
   % pinned separately by test_wet_growth_uses_both_liquid_branches.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=900);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=1, testname='forced_advance');
   f_liq_checkpoint = 0.01 * ones(size(state.f_liq));
   restart = struct('T', state.T, 'f_ice', state.f_ice, ...
      'f_liq', f_liq_checkpoint, 'Ts', state.Ts, 'r_eff', state.r_eff);
   restart_file = fullfile(workspace.rootdir, 'forced-advance-restart.mat');
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

   % Zero tolerance forces two nonlinear updates before rejection, leaving a
   % rejected-state conductivity to restore. The one-second full step makes
   % CHECKSUBSTEP take its bounded fallback immediately.
   testCase.verifyEqual(ice1.dt_sum, 1, 'AbsTol', 0);
   testCase.verifyEqual(ice1.n_subfail, 1, 'AbsTol', 0);
   testCase.verifyEqual(ice1.Tice_converged, 0, 'AbsTol', 0);
   testCase.verifyEqual(ice2.Tice, state.T, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_ice, state.f_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_liq, f_liq_checkpoint, 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_liq, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_evp, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_lyr, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.r_eff, state.r_eff, 'AbsTol', 0);

   % Final diagnostics must use the restored checkpoint's vapor-free solver
   % conductivity, not conductivity returned by the rejected solve.
   k_eff_checkpoint = icemodel.column.bulk_thermal_conductivity( ...
      state.T, state.f_ice, f_liq_checkpoint, 0);
   Qc_checkpoint = icemodel.surface.conductive_heat_flux( ...
      k_eff_checkpoint, state.T, state.dz, state.Ts);
   testCase.verifyEqual(ice1.Qc, Qc_checkpoint, 'AbsTol', 1e-12);

   % Every process ledger stays zero; only the unchanged storage endpoints are
   % nonzero because they describe the checkpoint itself.
   budget_fields = string(fieldnames(ice1));
   process_fields = budget_fields(startsWith(budget_fields, 'mass_budget_') ...
      & ~endsWith(budget_fields, {'_start_mwe', '_end_mwe'}));
   for k = 1:numel(process_fields)
      testCase.verifyEqual(ice1.(process_fields(k)), 0, 'AbsTol', 0, ...
         sprintf('%s changed during forced acceptance', process_fields(k)));
   end
   testCase.verifyEqual(ice1.mass_budget_solid_end_mwe, ...
      ice1.mass_budget_solid_start_mwe, 'AbsTol', 0);
   testCase.verifyEqual(ice1.mass_budget_liquid_end_mwe, ...
      ice1.mass_budget_liquid_start_mwe, 'AbsTol', 0);
end

function base_opts = syntheticRunOpts(testCase)
   %SYNTHETICRUNOPTS Build one short diagnostic production run.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=24, dt_seconds=3600);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   base_opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, output_profile='diagnostic', solver=1);
   assert(~isfield(base_opts, 'use_coupled_vapor'));
end

function [ok_seb, ok_ieb, ok_cpl] = runRobinProbe(s, cpl_alpha, cpl_aitken)
   %RUNROBINPROBE Exercise one explicit Robin coupling policy on a fixture.

   [~, ~, ~, ~, ~, ~, ok_seb, ok_ieb, ok_cpl] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.Ts, s.T, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, 3, s.tol, ...
      s.maxiter, s.alpha, s.use_aitken, s.jumpmax, s.cpl_Ts_tol, ...
      s.cpl_seb_tol, s.cpl_maxiter, cpl_alpha, cpl_aitken, ...
      s.cpl_jumpmax, s.ro_sfc, s.snow_depth, s.opts.f_res_pore_ice, ...
      s.opts);
end

function restoreProfiler(prior)
   %RESTOREPROFILER Put the profiler back in the caller's on/off state.

   profile off
   profile clear
   if strcmp(prior.ProfilerStatus, 'on')
      profile on
   end
end
