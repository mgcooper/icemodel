function tests = test_mesh_and_timestep_kernels
   %TEST_MESH_AND_TIMESTEP_KERNELS Verify mesh and timestep-update kernels.
   tests = functiontests(localfunctions);
end

function test_cvmesh_uniform_and_exponential_layout(testCase)
   % control_volume_mesh should preserve the requested depth while changing
   % only the spacing pattern between uniform and exponential layouts.

   [dz_u, ~, ~, z_edge_u, f_u] = ...
      icemodel.column.control_volume_mesh(1.0, 0.25);
   [dz_e, ~, ~, z_edge_e] = ...
      icemodel.column.control_volume_mesh(1.0, 0.10, 1.5);

   testCase.verifyEqual(sum(dz_u), 1.0, 'AbsTol', 1e-12);
   testCase.verifyEqual(z_edge_u(end), 1.0, 'AbsTol', 1e-12);
   testCase.verifyEqual(numel(f_u), numel(dz_u) + 1);
   testCase.verifyGreaterThanOrEqual(min(diff(z_edge_e)), 0);
   testCase.verifyGreaterThan(dz_e(end), dz_e(1));
end

function test_interp1_nearest_preserves_expected_spectral_remap(testCase)
   % The spectral density remap uses direct nearest-neighbor interp1; verify
   % the expected shape and values on a compact controlled example.

   ro_sno = interp1([0.2; 0.6; 1.0], [300; 400; 500], [0.1; 0.5; 0.9], ...
      'nearest', 'extrap');

   testCase.verifyEqual(numel(ro_sno), 3);
   testCase.verifyEqual(ro_sno, [300; 400; 500]);
end

function test_layerinds_selects_expected_merge_neighbors(testCase)
   % merge_layer_indices should pick the expected merge partner at the top,
   % over a zero-thickness layer, and for a nonzero interior layer.

   [j1_top, j2_top] = icemodel.column.merge_layer_indices(1, [0.0; 0.5; 0.4]);
   [j1_zero, j2_zero] = icemodel.column.merge_layer_indices(2, [0.4; 0.0; 0.2]);
   [j1_nonzero, j2_nonzero] = ...
      icemodel.column.merge_layer_indices(2, [0.0; 0.3; 0.5]);

   testCase.verifyEqual([j1_top j2_top], [1 2]);
   testCase.verifyEqual(j1_zero, 2);
   testCase.verifyEqual(j2_zero, 3);
   testCase.verifyEqual(j1_nonzero, 2);
   testCase.verifyEqual(j2_nonzero, 1);
end

function test_layerinds_merges_the_bottom_layer_upward(testCase)
   % The bottom layer has no layer below it, so it must merge with the one
   % above. Without this branch the interior path reads f_ice(j1 + 1) past the
   % end of the column.

   [j1, j2] = icemodel.column.merge_layer_indices(3, [0.4; 0.3; 0.0]);
   testCase.verifyEqual([j1 j2], [3 2]);

   % A nonzero bottom layer takes the same path, since the choice is forced by
   % the boundary rather than by the neighbour thicknesses.
   [j1_solid, j2_solid] = ...
      icemodel.column.merge_layer_indices(3, [0.4; 0.3; 0.2]);
   testCase.verifyEqual([j1_solid j2_solid], [3 2]);
end

function test_trisolve_matches_backslash(testCase)
   % The tridiagonal solver should reproduce MATLAB's dense solve on a
   % compact reference system.

   low = [0; -1; -1];
   mid = [4; 4; 4];
   upp = [-1; -1; 0];
   rhs = [2; 6; 2];
   A = [4 -1 0; -1 4 -1; 0 -1 4];

   x = icemodel.numerics.trisolve(low, mid, upp, rhs);

   testCase.verifyEqual(x, A \ rhs, 'AbsTol', 1e-12);
end

function test_conduct_matches_level_formulas(testCase)
   % conductive_heat_flux should return the top-boundary flux and its
   % T_sfc derivative matching the analytic finite-volume expressions.

   k_eff = [2; 4];
   T_ice = [270; 268];
   dz = [0.04; 0.04];
   T_sfc = 269;

   [Qc, dQc_dT_sfc] = icemodel.surface.conductive_heat_flux(k_eff, T_ice, ...
      dz, T_sfc);

   % Top-boundary flux: k_eff(1) * (T_ice(1) - T_sfc) / (dz(1)/2)
   testCase.verifyEqual(Qc, 2 * (270 - 269) / 0.02, 'AbsTol', 1e-12);

   % T_sfc derivative: -k_eff(1) / (dz(1)/2)
   testCase.verifyEqual(dQc_dT_sfc, -2 / 0.02, 'AbsTol', 1e-12);
end

function test_inittimesteps_and_newtimestep_follow_solver_contract(testCase)
   % Initialization and top-level timestep helpers must return the step
   % size, warmup count, and per-step reset state the drivers consume.

   opts = struct('dt', 900, 'numyears', 2, 'n_spinup_years', 1, ...
      'simyears', [2015 2016]);
   Time = transpose(datetime(2015, 1, 1) + minutes(15) * (0:7));

   [metstep, substep, numsteps, dt_new, numyears, numspinup, ...
      force_advance_streak_dt] = ...
      icemodel.timestepping.initialize_timesteps( ...
      opts, Time);

   [dt_sum, d_liq, d_evp, d_lyr, d_rof, d_vap_liq, ...
      d_vap_ice, diag] = icemodel.timestepping.newtimestep(zeros(3, 1));

   testCase.verifyEqual([metstep substep numsteps], [1 1 4]);
   testCase.verifyEqual([dt_new numyears numspinup], [900 2 1]);
   testCase.verifyEqual(force_advance_streak_dt, 0, 'AbsTol', 0);
   testCase.verifyEqual(dt_sum, 0);
   testCase.verifyEqual(diag.n_failed_substeps, 0);

   % No solver result is accepted at the start of a forcing step. ISEQUALN
   % treats the initial NaN values as equal.
   testCase.verifyTrue( ...
      isequaln(diag, icemodel.couplers.initialize_solver_diag()));

   % The couplers start each solve attempt from the raw attempt record, so
   % the second output must equal the substep field of the first.
   [returned_diag, returned_step_diag] = ...
      icemodel.couplers.initialize_solver_diag();
   testCase.verifyTrue(isequaln(returned_step_diag, returned_diag.substep));
   testCase.verifyEqual(d_liq, zeros(3, 1), 'AbsTol', 0);
   testCase.verifyEqual(d_evp, zeros(3, 1), 'AbsTol', 0);
   testCase.verifyEqual(d_lyr, zeros(3, 1), 'AbsTol', 0);
   testCase.verifyEqual(d_rof, 0.0, 'AbsTol', 0);
   testCase.verifyEqual(d_vap_liq, zeros(3, 1), 'AbsTol', 0);
   testCase.verifyEqual(d_vap_ice, zeros(3, 1), 'AbsTol', 0);
end

function test_nexttimestep_adapts_substep_divisor(testCase)
   % nexttimestep should shrink or grow the substep divisor based on the recent
   % convergence history and hard failures.

   settings = fixedNextTimestepSettings(900, 9);
   [~, substep_fast, dt_fast] = icemodel.timestepping.nexttimestep( ...
      1, 3, true, settings, fixedNextTimestepDiag(0, 1));
   [~, substep_slow, dt_slow] = icemodel.timestepping.nexttimestep( ...
      1, 3, true, settings, fixedNextTimestepDiag(2, 15));
   [~, substep_fail, dt_fail] = icemodel.timestepping.nexttimestep( ...
      1, 3, false, settings, fixedNextTimestepDiag(0, 0));

   testCase.verifyLessThan(substep_fast, 3);
   testCase.verifyGreaterThan(substep_slow, 3);
   testCase.verifyGreaterThan(substep_fail, 3);
   testCase.verifyGreaterThan(dt_fast, 300);
   testCase.verifyLessThan(dt_slow, 300);
   testCase.verifyLessThan(dt_fail, 300);
end

function test_resetsubstep_and_acceptsubstep_restore_and_advance(testCase)
   % The reset helper restores failed-substep state; the accept helper
   % checkpoints the accepted state, advances time consistently, records
   % the accepted solve, and restores the primary settings.

   Ls = icemodel.physicalConstant('Ls');
   ro_atm_val = 1.2;
   De_e_val = 1e-5;
   substep_opts = struct( ...
      'f_res_pore_snow', 0.07, 'f_res_pore_ice', 0.01);

   [T_sfc, T_ice, f_ice, f_liq, k_eff, n_subfail, substep, dt_new] = ...
      icemodel.timestepping.resetsubstep( ...
      270, [269; 268], [0.9; 0.9], [0.01; 0.01], [2.1; 2.2], ...
      900, 2, 9, 0, 450);

   % The accepted solve used the recovery settings, so acceptsubstep must
   % credit cpl_recovery_count and then hand back the primary settings.
   [settings, settings0] = defaultCheckSettings();
   settings.cpl_recovery_active = true;
   diag = fixedDiag(true, true, true);
   diag.substep.n_iters = 4;

   [T_sfc_up, T_ice_up, f_ice_up, f_liq_up, k_eff_up, dt_sum, dt_next, ...
      settings_up, diag_up] = icemodel.timestepping.acceptsubstep(T_sfc, ...
      T_ice, f_ice, f_liq, k_eff, 450, 300, 1e-12, settings, settings0, diag);

   testCase.verifyEqual([T_sfc_up; T_ice_up], [270; 269; 268], 'AbsTol', 0);
   testCase.verifyEqual([f_ice_up; f_liq_up], [0.9; 0.9; 0.01; 0.01], ...
      'AbsTol', 0);
   testCase.verifyEqual(k_eff_up, [2.1; 2.2], 'AbsTol', 0);
   testCase.verifyEqual(n_subfail, 1);
   testCase.verifyEqual(substep, 3);
   testCase.verifyEqual(dt_new, 300, 'AbsTol', 1e-12);
   testCase.verifyEqual(dt_sum, 750, 'AbsTol', 1e-12);
   testCase.verifyEqual(dt_next, 150, 'AbsTol', 1e-12);

   % The accepted solve replaces the default diagnostic values and credits the
   % recovery count because SETTINGS marked this substep cpl_recovery_active.
   testCase.verifyEqual(diag_up.ok_cpl, true);
   testCase.verifyEqual(diag_up.n_iters, 4);
   testCase.verifyEqual(diag_up.cpl_recovery_count, 1.0, 'AbsTol', 0);

   % acceptsubstep restores the primary settings after recording the solve.
   testCase.verifyEqual(settings_up, settings0);

   % update_surface_state derives the surface running state.
   [liqflag, ~, hv_atm_val, H_e, ~] = ...
      icemodel.surface.update_surface_state( ...
      f_ice(1), f_liq(1), ro_atm_val, De_e_val, 0, substep_opts);

   testCase.verifyFalse(liqflag);
   testCase.verifyEqual(hv_atm_val, ro_atm_val * Ls, 'AbsTol', 0);
   testCase.verifyEqual(H_e, hv_atm_val * De_e_val, 'AbsTol', 0);
end

function test_checksubstep_forces_advance_at_maxsubstep(testCase)
   % Once the max-substep limit is reached, checksubstep should force the
   % accepted state forward instead of stalling the timestep.

   settings = defaultCheckSettings();
   settings.maxsubstep = 2;
   diag = fixedDiag(false, false, false);
   diag.n_failed_substeps = 1;

   [T_sfc, T_ice, f_ice, f_liq, k_eff, substep, dt_new, ok, forced_advance, ...
      ~, ~, diag] = icemodel.timestepping.checksubstep(270, [269; 268], ...
      [0.9; 0.9], [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 150, 450, 1, 10, 1, 0.0, 'icemodel', ...
      settings, settings, diag);

   testCase.verifyTrue(ok);
   testCase.verifyEqual(T_sfc, 271);
   testCase.verifyEqual(T_ice, [270; 269]);
   testCase.verifyEqual(f_ice, [0.8; 0.8]);
   testCase.verifyEqual(f_liq, [0.02; 0.02]);
   testCase.verifyEqual(k_eff, [1; 2], 'AbsTol', 0);
   testCase.verifyEqual(diag.n_failed_substeps, 2);
   testCase.verifyEqual(substep, 2);
   testCase.verifyEqual(dt_new, 450, 'AbsTol', 1e-12);
   testCase.verifyTrue(forced_advance);
end

function test_checksubstep_retries_failed_solve_from_checkpoint(testCase)
   % Before maxsubstep, a rejected inner solve must restore the accepted
   % checkpoint, shorten dt, and remain rejected so the caller retries.

   settings = defaultCheckSettings();
   settings.maxsubstep = 4;
   diag = fixedDiag(false, false, false);
   diag.n_failed_substeps = 0;

   [T_sfc, T_ice, f_ice, f_liq, k_eff, substep, dt_new, ok, forced_advance, ...
      ~, ~, diag] = icemodel.timestepping.checksubstep(270, [269; 268], ...
      [0.9; 0.9], [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 0, 900, 1, 10, 1, 0.0, 'icemodel', ...
      settings, settings, diag);

   testCase.verifyFalse(ok);
   testCase.verifyEqual(T_sfc, 271);
   testCase.verifyEqual(T_ice, [270; 269]);
   testCase.verifyEqual(f_ice, [0.8; 0.8]);
   testCase.verifyEqual(f_liq, [0.02; 0.02]);
   testCase.verifyEqual(k_eff, [1; 2], 'AbsTol', 0);
   testCase.verifyEqual(diag.n_failed_substeps, 1);
   testCase.verifyEqual(substep, 2);
   testCase.verifyEqual(dt_new, 450, 'AbsTol', 1e-12);
   testCase.verifyFalse(forced_advance);
end

function test_checksubstep_seb_failure_skips_settings_retry(testCase)
   % A failed SEB solve goes straight to dt shortening: the settings retry
   % changes only the outer relaxation, which cannot repair the SEB.

   settings = defaultCheckSettings();
   settings.maxsubstep = 4;
   diag = fixedDiag(false, true, false);
   diag.n_failed_substeps = 0;

   [~, ~, ~, ~, ~, ~, dt_new, ok, forced_advance, ~, settings, ~] = ...
      icemodel.timestepping.checksubstep(270, [269; 268], [0.9; 0.9], ...
      [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 0, 900, 1, 10, 1, 0.0, 'icemodel', ...
      settings, settings, diag);

   testCase.verifyFalse(ok);
   testCase.verifyFalse(forced_advance);
   testCase.verifyEqual(dt_new, 450, 'AbsTol', 1e-12);
   testCase.verifyFalse(settings.cpl_recovery_active);
end

function test_checksubstep_recovery_primary_skips_settings_retry(testCase)
   % Skip recovery mode when it matches the primary settings. The next
   % attempt uses a shorter timestep instead of repeating the same solve.

   settings = defaultCheckSettings();
   settings.maxsubstep = 4;
   settings.cpl_aitken = false;
   settings.cpl_alpha = settings.cpl_recovery_alpha;
   diag = fixedDiag(true, true, false);
   diag.n_failed_substeps = 0;

   [~, ~, ~, ~, ~, ~, dt_new, ok, forced_advance, ~, settings, ~] = ...
      icemodel.timestepping.checksubstep(270, [269; 268], [0.9; 0.9], ...
      [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 0, 900, 1, 10, 1, 0.0, 'icemodel', ...
      settings, settings, diag);

   testCase.verifyFalse(ok);
   testCase.verifyFalse(forced_advance);
   testCase.verifyEqual(dt_new, 450, 'AbsTol', 1e-12);
   testCase.verifyFalse(settings.cpl_recovery_active);
end

function test_checksubstep_clamps_overshot_failure_count(testCase)
   % Even if a caller enters checksubstep with an already-overshot failure
   % count, the timestep controller should clamp to the accepted dt_min state
   % and force advance instead of stalling forever at dt_min.

   settings = defaultCheckSettings();
   settings.maxsubstep = 10;
   diag = fixedDiag(false, false, false);
   diag.n_failed_substeps = 10;

   [T_sfc, T_ice, f_ice, f_liq, k_eff, substep, dt_new, ok, forced_advance, ...
      ~, ~, diag] = icemodel.timestepping.checksubstep(270, [269; 268], ...
      [0.9; 0.9], [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 150, 450, 1, 10, 10, 0.0, 'icemodel', ...
      settings, settings, diag);

   testCase.verifyTrue(ok);
   testCase.verifyEqual(T_sfc, 271);
   testCase.verifyEqual(T_ice, [270; 269]);
   testCase.verifyEqual(f_ice, [0.8; 0.8]);
   testCase.verifyEqual(f_liq, [0.02; 0.02]);
   testCase.verifyEqual(k_eff, [1; 2], 'AbsTol', 0);
   testCase.verifyEqual(diag.n_failed_substeps, 10);
   testCase.verifyEqual(substep, 10);
   testCase.verifyEqual(dt_new, 90, 'AbsTol', 1e-12);
   testCase.verifyTrue(forced_advance);
end

function test_checksubstep_debug_dump_records_force_advance_context(testCase)
   % The maxsubstep debug dump should record whether the current failure
   % triggered force advance and the projected cross-timestep streak.

   debug_file = [tempname '.mat'];
   cleanup = onCleanup(@() cleanupDebugFile(debug_file));
   setenv('ICEMODEL_DEBUG_MAXSUBSTEP_FILE', debug_file);

   settings = defaultCheckSettings();
   settings.maxsubstep = 2;
   settings.debug = true;
   diag = fixedDiag(false, false, false);
   diag.n_failed_substeps = 1;

   icemodel.timestepping.checksubstep(270, [269; 268], [0.9; 0.9], ...
      [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 150, 450, 1, 10, 1, 450, 'icemodel', ...
      settings, settings, diag);

   loaded = load(debug_file, 'debug_state');
   debug_state = loaded.debug_state;

   testCase.verifyTrue(debug_state.forced_advance);
   testCase.verifyEqual(debug_state.projected_force_advance_dt, 900, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(debug_state.dt_full_step, 900, 'AbsTol', 1e-12);
   testCase.verifyEqual(debug_state.failed.k_eff, [3; 4], 'AbsTol', 0);
   testCase.verifyEqual(debug_state.checkpoint.k_eff, [1; 2], 'AbsTol', 0);
end

function test_force_advance_streak_resets_on_acceptance(testCase)
   % An accepted substep clears any earlier force-advance streak, so a
   % transient recovery does not affect later timesteps. CHECKSUBSTEP owns the
   % streak guard.

   settings = defaultCheckSettings();
   settings.maxsubstep = 2;
   diag = fixedDiag(true, true, true);
   diag.n_failed_substeps = 1;

   [~, ~, ~, ~, ~, ~, ~, ok, forced_advance, streak_dt] = ...
      icemodel.timestepping.checksubstep(270, [269; 268], [0.9; 0.9], ...
      [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 150, 450, 1, 10, 1, 300, 'icemodel', ...
      settings, settings, diag);

   testCase.verifyTrue(ok);
   testCase.verifyFalse(forced_advance);
   testCase.verifyEqual(streak_dt, 0, 'AbsTol', 0);
end

function test_force_advance_streak_errors_after_full_timestep(testCase)
   % Consecutive forced-advance time beyond one full forcing step fails
   % fast instead of letting a long broken run integrate checkpoints
   % forward. The failing call below forces an advance with the streak
   % already at the one-step limit.

   settings = defaultCheckSettings();
   settings.maxsubstep = 2;
   diag = fixedDiag(false, false, false);
   diag.n_failed_substeps = 1;

   testCase.verifyError(@() icemodel.timestepping.checksubstep( ...
      270, [269; 268], [0.9; 0.9], ...
      [0.01; 0.01], [3; 4], 271, [270; 269], [0.8; 0.8], ...
      [0.02; 0.02], [1; 2], 150, 450, 1, 10, 1, 900, 'icemodel', ...
      settings, settings, diag), 'icemodel:ForceAdvanceStreakExceeded');
end

function cleanupDebugFile(debug_file)
   %CLEANUPDEBUGFILE Restore the debug env var and remove the temp MAT file.

   setenv('ICEMODEL_DEBUG_MAXSUBSTEP_FILE', '');
   if exist(debug_file, 'file') == 2
      delete(debug_file);
   end
end

function diag = fixedDiag(ok_seb, ok_ieb, ok_cpl)
   %FIXEDDIAG Set the three attempt flags that checksubstep reads.

   diag = icemodel.couplers.initialize_solver_diag();
   diag.substep.ok_seb = ok_seb;
   diag.substep.ok_ieb = ok_ieb;
   diag.substep.ok_cpl = ok_cpl;
end

function [settings, settings0] = defaultCheckSettings()
   %DEFAULTCHECKSETTINGS Build a minimal settings struct for timestep probes.
   % dt=900 matches the probes. Each test can override maxsubstep or the
   % coupler recovery pair.

   opts = struct('solver', 3, 'maxiter', 50, 'tol', 1e-6, 'alpha', 1, ...
      'use_aitken', true, 'jumpmax', 10, 'cpl_maxiter', 20, ...
      'cpl_Ts_tol', 1e-3, 'cpl_seb_tol', 1e-2, 'cpl_alpha', 0.5, ...
      'cpl_aitken', true, 'cpl_jumpmax', 10, 'dt', 900, 'debug', false);
   [settings, settings0] = icemodel.couplers.initialize_solver_settings(opts);
end

function settings = fixedNextTimestepSettings(dt_full_step, maxsubstep)
   %FIXEDNEXTTIMESTEPSETTINGS Set the two fields nexttimestep reads.

   settings = struct('dt_full_step', dt_full_step, 'maxsubstep', maxsubstep);
end

function diag = fixedNextTimestepDiag(n_failed_substeps, n_iters)
   %FIXEDNEXTTIMESTEPDIAG Set the two diagnostic fields nexttimestep reads.

   diag = icemodel.couplers.initialize_solver_diag();
   diag.n_failed_substeps = n_failed_substeps;
   diag.substep.n_iters = n_iters;
end

function test_bottom_layer_merge_removes_it_and_conserves_mass(testCase)
   % A deepest layer below f_ice_min must actually be removed. Take the clone
   % that preserves column length AFTER the deletion. Taking it first would
   % copy the removed layer back into the column, and the layer above it would
   % then hold only half the pair's mass.

   [ro_ice, ro_liq, Tf] = icemodel.physicalConstant('ro_ice', 'ro_liq', 'Tf');
   dz = 0.04;
   f_ice_min = 0.1;
   T_ice = [Tf - 1; Tf - 2; Tf - 3];
   f_ice = [0.9; 0.8; 0.02];
   f_liq = [0.01; 0.01; 0.0];
   zeros_col = zeros(3, 1);

   water_equivalent = @(fi, fl) sum(ro_ice / ro_liq * fi + fl) * dz;
   expected = water_equivalent(f_ice, f_liq);

   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);
   [~, returned_f_ice, returned_f_liq] = icemodel.column.merge_thin_layers( ...
      T_ice, f_ice, f_liq, zeros_col, zeros_col, dz, 0.0, zeros_col, ...
      f_ice_min, budget);

   testCase.verifyFalse(any(returned_f_ice < f_ice_min))
   testCase.verifyEqual( ...
      water_equivalent(returned_f_ice, returned_f_liq), expected, ...
      AbsTol=1e-12)
   testCase.verifyNumElements(returned_f_ice, 3)
end
