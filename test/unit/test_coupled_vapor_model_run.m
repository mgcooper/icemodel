function tests = test_coupled_vapor_model_run
   %TEST_COUPLED_VAPOR_MODEL_RUN Verify production vapor coupling end to end.
   %
   % The kernel tests isolate face transport and constrained phase transfer.
   % These tests run the model. They verify that the retired flag cannot
   % select alternate physics. They also verify surface-vapor closure and
   % per-phase storage closure. A kernel-level probe checks that
   % icemodel.column.couple_vapor_step's transport channels register a
   % clamped transfer, because a short production run gives no guarantee a
   % storage-limit clamp actually fires: absent a clamp the donor-receiver
   % scheme is exact and the net transport channels stay at zero.
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
   % Potential surface energy equals realized phase storage and
   % condensation overflow. Interior transport stays separate.

   base_opts = syntheticRunOpts(testCase);
   [Ls, Lv, ro_liq] = icemodel.physicalConstant('Ls', 'Lv', 'ro_liq');
   [ice1, ~] = icemodel.test.helpers.runSmbModel(base_opts);

   accounted = ro_liq * ( ...
      Ls * ice1.mass_budget_vapor_solid_mwe ...
      + Lv * ice1.mass_budget_vapor_liquid_mwe ...
      + Lv * ice1.mass_budget_condensation_overflow_mwe);
   potential = ice1.mass_budget_vapor_potential_j_m2;

   % Allow only column-integration subtraction roundoff.
   testCase.verifyEqual(potential, accounted, 'AbsTol', 1e-4);
   scale = max(abs(potential));
   testCase.verifyGreaterThan(scale, 0);
   testCase.verifyLessThan(max(abs(potential - accounted)) / scale, 1e-8);
end

function test_production_run_closes_the_per_phase_storage(testCase)
   % Endpoint phase storage must equal phase change, surface vapor exchange,
   % remeshing, and interior vapor transport together.

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
      + ice1.mass_budget_vapor_transport_solid_mwe, 'AbsTol', 1e-10);
   testCase.verifyEqual(liquid_delta, ...
      ice1.mass_budget_phase_liquid_mwe ...
      + ice1.mass_budget_vapor_liquid_mwe ...
      + ice1.mass_budget_remesh_liquid_mwe ...
      + ice1.mass_budget_vapor_transport_liquid_mwe, 'AbsTol', 1e-10);
end

function test_saved_conductivity_reproduces_the_conductive_flux(testCase)
   % Vapor exchange changes the phase state after the thermal solve. The saved
   % conductivity must use that state so it reproduces the model's heat flux.

   base_opts = syntheticRunOpts(testCase);
   [ice1, ice2, opts] = icemodel.test.helpers.runSmbModel(base_opts);
   k_eff = icemodel.column.bulk_thermal_conductivity( ...
      ice2.Tice, ice2.f_ice, ice2.f_liq, 0);
   Qc_expected = arrayfun(@(step) ...
      icemodel.surface.conductive_heat_flux(k_eff(:, step), ...
      ice2.Tice(:, step), opts.dz_thermal, ice1.Tsfc(step)), ...
      (1:numel(ice1.Tsfc))');

   testCase.verifyEqual(ice1.Qc, Qc_expected, 'AbsTol', 1e-12);
end

function test_couple_vapor_step_records_nonzero_transport_net(testCase)
   % icemodel.column.couple_vapor_step must show up in the surviving net
   % transport channels once a receiving node's storage-limit clamp blocks
   % part of a face transfer.
   %
   % Absent a clamp, the donor-receiver scheme is exact: a face's donor
   % loses what its receiver gains, in the same phase, so the column-total
   % net for that phase stays at zero. Nothing forces a clamp within a
   % short synthetic production run, so that run is not a reliable witness
   % for this channel. This kernel-level probe forces a clamp directly, so
   % the nonzero net is provable rather than assumed.

   [Lv, ro_liq] = icemodel.physicalConstant('Lv', 'ro_liq');

   % Two-node column. Node 1 is a generous liquid donor. Node 2 is fully
   % solid ice with no liquid pore space (f_ice = 1, f_liq = 0), so its
   % storage-limit clamp accepts none of an incoming liquid transfer.
   T_ice = [260; 260];
   f_ice = [0.3; 1.0];
   f_liq = [0.3; 0.0];
   dz = [0.05; 0.05];
   dt = 900;
   f_ice_min = 0.01;
   f_res_por = 0.02;

   % One interior face (JJ = 2 gives faces 1 and 3 as boundaries, face 2 as
   % the only interior face). L_vap(2) = Lv routes the transfer through the
   % liquid phase at both the donor (node 1) and the receiver (node 2), so
   % the ice channel is never touched and the solid net must stay exactly
   % zero.
   JJ = numel(f_ice);
   U_vap = zeros(JJ + 1, 1);
   U_vap(2) = 1e-4;
   L_vap = Lv * ones(JJ + 1, 1);
   d_vap_liq = zeros(JJ, 1);
   d_vap_ice = zeros(JJ, 1);

   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);

   [returned_f_ice, returned_f_liq, ~, ~, returned_budget] = ...
      icemodel.column.couple_vapor_step(f_ice, f_liq, U_vap, L_vap, ...
      d_vap_liq, d_vap_ice, dz, dt, f_ice_min, f_res_por, budget);

   % Node 2's liquid request is fully clamped, so its state is untouched
   % and node 1 gives up the full requested amount. The ice channel never
   % receives a request, so f_ice does not change at all.
   d_vap_face = U_vap(2) * dt / ro_liq;
   expected_f_liq = [f_liq(1) - d_vap_face / dz(1); f_liq(2)];
   testCase.verifyEqual(returned_f_liq, expected_f_liq, 'AbsTol', 1e-15);
   testCase.verifyEqual(returned_f_ice, f_ice, 'AbsTol', 0);

   % The clamp breaks the exact donor/receiver cancellation, so the liquid
   % net records the mass node 2 could not accept, and it must be nonzero.
   % The solid net is untouched and must stay exactly zero.
   expected_liquid_net = -d_vap_face;
   testCase.verifyEqual( ...
      returned_budget.mass_budget_vapor_transport_liquid_mwe, ...
      expected_liquid_net, 'AbsTol', 1e-15);
   testCase.verifyNotEqual( ...
      returned_budget.mass_budget_vapor_transport_liquid_mwe, 0);
   testCase.verifyEqual( ...
      returned_budget.mass_budget_vapor_transport_solid_mwe, 0, ...
      'AbsTol', 0);
end

function test_robin_recovery_converges_via_checksubstep_retry(testCase)
   % A healthy-inner Robin outer failure on the primary settings (ok_ieb
   % true, ok_cpl false) makes checksubstep restore the checkpoint and
   % return the recovery settings (cpl_alpha == cpl_recovery_alpha,
   % cpl_aitken false, cpl_recovery_active true) at the same dt with
   % diag.n_failed_substeps unchanged. Rerunning the coupler with those
   % settings converges.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=2, dt_seconds=900);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='robin_recovery');

   % Alpha 2 without acceleration deterministically leaves the synthetic
   % inner solve healthy while the primary outer loop exhausts its
   % iterations.
   primary_settings = state.settings;
   primary_settings.cpl_alpha = 2.0;
   primary_settings.cpl_aitken = false;
   [T_sfc1, T_ice1, f_ice1, f_liq1, k_eff1, ~, ~, diag1] = ...
      runRobinProbe(state, primary_settings);
   testCase.verifyTrue(diag1.ok_seb);
   testCase.verifyTrue(diag1.ok_ieb);
   testCase.verifyFalse(diag1.ok_cpl);

   % dt_full_step already matches this fixture's build-time opts.dt, so only
   % maxsubstep and debug need the small-probe override after the
   % initializer.
   dt = state.opts.dt;
   primary_settings.maxsubstep = 4;
   primary_settings.debug = false;
   step_diag1 = icemodel.couplers.initialize_solver_diag();
   step_diag1.substep = diag1;
   [T_sfc_r, T_ice_r, f_ice_r, f_liq_r, k_eff_r, ~, dt_out, ok, ...
      forced_advance, ~, retry_settings, diag_out] = ...
      icemodel.timestepping.checksubstep(T_sfc1, T_ice1, f_ice1, f_liq1, ...
      k_eff1, state.T_sfc, state.T_ice, state.f_ice, state.f_liq, ...
      state.k_eff, 0.0, dt, 1, 2, 1, 0.0, 'test', primary_settings, ...
      primary_settings, step_diag1);

   % checksubstep restores the checkpoint and offers the recovery-mode
   % retry without shortening dt or charging a failure.
   testCase.verifyFalse(ok);
   testCase.verifyFalse(forced_advance);
   testCase.verifyEqual(dt_out, dt, 'AbsTol', 0);
   testCase.verifyEqual(diag_out.n_failed_substeps, 0, 'AbsTol', 0);
   testCase.verifyTrue(retry_settings.cpl_recovery_active);
   testCase.verifyEqual(retry_settings.cpl_alpha, ...
      primary_settings.cpl_recovery_alpha, 'AbsTol', 0);
   testCase.verifyFalse(retry_settings.cpl_aitken);
   testCase.verifyEqual(T_sfc_r, state.T_sfc, 'AbsTol', 0);
   testCase.verifyEqual(T_ice_r, state.T_ice, 'AbsTol', 0);
   testCase.verifyEqual(f_ice_r, state.f_ice, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_r, state.f_liq, 'AbsTol', 0);
   testCase.verifyEqual(k_eff_r, state.k_eff, 'AbsTol', 0);

   % Rerunning the coupler with the recovery settings converges.
   [~, ~, ~, ~, ~, ~, ~, diag2] = runRobinProbe(state, retry_settings);
   testCase.verifyTrue(diag2.ok_seb);
   testCase.verifyTrue(diag2.ok_ieb);
   testCase.verifyTrue(diag2.ok_cpl);
end

function test_robin_recovery_rescues_the_last_allowed_attempt(testCase)
   % With recovery inside checksubstep's settings retry, a healthy-inner
   % outer failure on the only allowed substep attempt recovers instead of
   % forcing an advance: the retry runs at the same dt and charges no
   % failure, so it never spends the last allowed attempt.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=900);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='robin_recovery_boundary');

   % DT=1 makes MAXSUBSTEP=1: one attempt exists. Alpha 2 leaves the inner
   % solve healthy while the primary outer loop fails; the settings retry
   % rescues that only attempt.
   opts = state.opts;
   opts.dt = 1;
   opts.cpl_alpha = 2.0;
   opts.cpl_aitken = false;
   opts.output_profile = 'diagnostic';
   opts.vars1 = {};
   opts.vars2 = {};
   [ice1, ~] = icemodel.test.helpers.runSmbModel(opts);

   testCase.verifyEqual(ice1.cpl_recovery_count, 1);
   testCase.verifyEqual(ice1.n_failed_substeps, 0);
   testCase.verifyEqual(ice1.dt_sum, 1);
   testCase.verifyEqual(ice1.Tsfc_converged, 1);
   testCase.verifyEqual(ice1.Tice_converged, 1);
end

function test_failed_robin_recovery_falls_through_substep_control(testCase)
   % The recovery-mode retry can fail too. checksubstep then offers
   % the retry only once: a second consecutive healthy-inner failure under
   % recovery settings uses the ordinary dt reduction instead of another
   % retry. The model substep loop then proceeds through reduced-step and
   % forced-advance handling.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=900);
   testCase.addTeardown(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='robin_recovery_failure');

   % Two strict outer sweeps keep both policies nonconverged while their
   % inner solves remain healthy. DT=2 permits one retry before the next
   % failure reaches the forced-advance boundary.
   probe = state;
   probe.opts.dt = 2;
   probe.settings.cpl_maxiter = 2;
   probe.settings.cpl_Ts_tol = 0;
   probe.settings.cpl_seb_tol = 0;
   primary_settings = probe.settings;
   primary_settings.cpl_alpha = 2.0;
   primary_settings.cpl_aitken = false;

   [T_sfc1, T_ice1, f_ice1, f_liq1, k_eff1, ~, ~, diag1] = ...
      runRobinProbe(probe, primary_settings);
   testCase.verifyTrue(diag1.ok_seb && diag1.ok_ieb && ~diag1.ok_cpl);

   % probe.settings was resolved at the fixture's original build-time
   % opts.dt, before probe.opts.dt was overridden above, so dt_full_step
   % needs the same override as maxsubstep and debug.
   dt = probe.opts.dt;
   primary_settings.dt_full_step = dt;
   primary_settings.maxsubstep = 4;
   primary_settings.debug = false;
   step_diag1 = icemodel.couplers.initialize_solver_diag();
   step_diag1.substep = diag1;
   [~, ~, ~, ~, ~, ~, dt_out1, ok1, forced1, ~, retry_settings, ...
      diag_out1] = icemodel.timestepping.checksubstep( ...
      T_sfc1, T_ice1, f_ice1, f_liq1, k_eff1, probe.T_sfc, probe.T_ice, ...
      probe.f_ice, probe.f_liq, probe.k_eff, 0.0, dt, 1, 1, 1, 0.0, 'test', ...
      primary_settings, primary_settings, step_diag1);
   testCase.verifyFalse(ok1);
   testCase.verifyFalse(forced1);
   testCase.verifyEqual(dt_out1, dt, 'AbsTol', 0);
   testCase.verifyEqual(diag_out1.n_failed_substeps, 0, 'AbsTol', 0);
   testCase.verifyTrue(retry_settings.cpl_recovery_active);

   % The recovery-mode retry also fails on this deliberately oscillating
   % fixture.
   [T_sfc3, T_ice3, f_ice3, f_liq3, k_eff3, ~, ~, diag2] = ...
      runRobinProbe(probe, retry_settings);
   testCase.verifyTrue(diag2.ok_seb && diag2.ok_ieb && ~diag2.ok_cpl);

   % cpl_recovery_active is already true, so the guard blocks a second
   % retry: checksubstep shortens dt instead and charges the failure.
   step_diag2 = icemodel.couplers.initialize_solver_diag();
   step_diag2.substep = diag2;
   [~, ~, ~, ~, ~, ~, dt_out2, ok2, forced2, ~, settings_out, ...
      diag_out2] = icemodel.timestepping.checksubstep( ...
      T_sfc3, T_ice3, f_ice3, f_liq3, k_eff3, probe.T_sfc, probe.T_ice, ...
      probe.f_ice, probe.f_liq, probe.k_eff, 0.0, dt, 1, 1, 1, 0.0, 'test', ...
      retry_settings, primary_settings, step_diag2);
   testCase.verifyFalse(ok2);
   testCase.verifyFalse(forced2);

   % A dt reset returns the primary settings.
   testCase.verifyEqual(settings_out, primary_settings);
   testCase.verifyEqual(diag_out2.n_failed_substeps, 1, 'AbsTol', 0);
   testCase.verifyLessThan(dt_out2, dt);

   % Through the production driver, the same failure pattern falls through
   % to the ordinary reduced-step and forced-advance path: no substep in
   % this forcing step ever reaches ok_cpl, so nothing physical applies.
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

   [ice1, ice2, returned_opts] = icemodel.test.helpers.runSmbModel(opts);

   testCase.verifyEqual(ice1.cpl_recovery_count, 0);
   testCase.verifyEqual(ice1.n_failed_substeps, 2);
   testCase.verifyEqual(ice1.dt_sum, 2);
   testCase.verifyEqual(ice1.Tsfc, state.T_sfc, 'AbsTol', 0);
   testCase.verifyEqual(ice2.Tice, state.T_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_ice, state.f_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_liq, state.f_liq, 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_liq, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_evp, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_lyr, zeros(size(state.f_liq)), 'AbsTol', 0);

   % Solver policy is not prognostic state: no latch rides the restart.
   restart_file = icemodel.restartfile(returned_opts, 2016);
   saved = load(restart_file, 'restart');
   testCase.verifyFalse(isfield(saved.restart, 'use_conservative_cpl'));

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
   restart = struct('T_ice', state.T_ice, 'f_ice', state.f_ice, ...
      'f_liq', f_liq_checkpoint, 'T_sfc', state.T_sfc, 'r_eff', state.r_eff);
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
   testCase.verifyEqual(ice1.n_failed_substeps, 1, 'AbsTol', 0);

   % A step whose only content was a forced advance reports the fixed
   % defaults from initialize_solver_diag, never the rejected solve's
   % values: converged flags false, iteration counts and residuals NaN or
   % zero. update_solver_diag ignores a non-accepted record, so a forcing
   % step of only forced advances never overwrites the defaults.
   testCase.verifyEqual(ice1.Tsfc_converged, 0, 'AbsTol', 0);
   testCase.verifyEqual(ice1.Tice_converged, 0, 'AbsTol', 0);
   testCase.verifyTrue(isnan(ice1.Tice_numiter));
   testCase.verifyEqual(ice1.cpl_iters, 0, 'AbsTol', 0);
   testCase.verifyTrue(isnan(ice1.cpl_res));
   testCase.verifyTrue(isnan(ice1.seb_res));
   testCase.verifyEqual(ice1.cpl_recovery_count, 0, 'AbsTol', 0);
   testCase.verifyEqual(ice2.Tice, state.T_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_ice, state.f_ice, 'AbsTol', 0);
   testCase.verifyEqual(ice2.f_liq, f_liq_checkpoint, 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_liq, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_evp, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.df_lyr, zeros(size(state.f_liq)), 'AbsTol', 0);
   testCase.verifyEqual(ice2.r_eff, state.r_eff, 'AbsTol', 0);

   % Final diagnostics must use the restored checkpoint's vapor-free solver
   % conductivity, not conductivity returned by the rejected solve.
   k_eff_checkpoint = icemodel.column.bulk_thermal_conductivity( ...
      state.T_ice, state.f_ice, f_liq_checkpoint, 0);
   Qc_checkpoint = icemodel.surface.conductive_heat_flux( ...
      k_eff_checkpoint, state.T_ice, state.dz, state.T_sfc);
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

function [T_sfc, T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, diag] = ...
      runRobinProbe(s, settings)
   %RUNROBINPROBE Exercise the Robin coupler once under one settings policy.

   [T_sfc, T_ice, f_ice, f_liq, k_eff, U_vap, L_vap, diag] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings, s.opts);
end
