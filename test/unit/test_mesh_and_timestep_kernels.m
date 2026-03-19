function tests = test_mesh_and_timestep_kernels
   %TEST_MESH_AND_TIMESTEP_KERNELS Verify mesh and timestep-update kernels.
   tests = functiontests(localfunctions);
end

function test_cvmesh_uniform_and_exponential_layout(testCase)
   % CVMESH should preserve the requested depth while changing only the
   % spacing pattern between uniform and exponential layouts.

   [dz_u, ~, ~, z_edge_u, f_u] = CVMESH(1.0, 0.25);
   [dz_e, ~, ~, z_edge_e] = CVMESH(1.0, 0.10, 1.5);

   testCase.verifyEqual(sum(dz_u), 1.0, 'AbsTol', 1e-12);
   testCase.verifyEqual(z_edge_u(end), 1.0, 'AbsTol', 1e-12);
   testCase.verifyEqual(numel(f_u), numel(dz_u) + 1);
   testCase.verifyGreaterThanOrEqual(min(diff(z_edge_e)), 0);
   testCase.verifyGreaterThan(dz_e(end), dz_e(1));
end

function test_gridforward_and_gridinverse_transform_shapes(testCase)
   % The spectral grid forward/inverse helpers should preserve the expected
   % vector sizes on simple controlled inputs.

   ro_sno = GRIDFORWARD([300; 400; 500], [0.2; 0.6; 1.0], [0.1; 0.5; 0.9]);
   dQp = GRIDINVERSE([10; 9; 8; 7], 0.5, 1.0, 4, 2);

   testCase.verifyEqual(numel(ro_sno), 3);
   testCase.verifyEqual(ro_sno, [300; 400; 500]);
   testCase.verifyEqual(dQp, [2; 2], 'AbsTol', 1e-12);
end

function test_layerinds_selects_expected_merge_neighbors(testCase)
   % LAYERINDS should pick the expected merge partner at the top, over a
   % zero-thickness layer, and for a nonzero interior layer.

   [j1_top, j2_top] = LAYERINDS(1, [0.0; 0.5; 0.4]);
   [j1_zero, j2_zero] = LAYERINDS(2, [0.4; 0.0; 0.2]);
   [j1_nonzero, j2_nonzero] = LAYERINDS(2, [0.0; 0.3; 0.5]);

   testCase.verifyEqual([j1_top j2_top], [1 2]);
   testCase.verifyEqual(j1_zero, 2);
   testCase.verifyEqual(j2_zero, 3);
   testCase.verifyEqual(j1_nonzero, 2);
   testCase.verifyEqual(j2_nonzero, 1);
end

function test_trisolve_matches_backslash(testCase)
   % The tridiagonal solver should reproduce MATLAB's dense solve on a
   % compact reference system.

   low = [0; -1; -1];
   mid = [4; 4; 4];
   upp = [-1; -1; 0];
   rhs = [2; 6; 2];
   A = [4 -1 0; -1 4 -1; 0 -1 4];

   x = TRISOLVE(low, mid, upp, rhs);

   testCase.verifyEqual(x, A \ rhs, 'AbsTol', 1e-12);
end

function test_conduct_matches_level_formulas(testCase)
   % CONDUCT should reproduce the top-boundary and interior finite-volume
   % forms used elsewhere in the column model.

   k_eff = [2; 4];
   T = [270; 268];
   dz = [0.04; 0.04];
   Ts = 269;

   Qc_top = CONDUCT(k_eff, T, dz, Ts, 1);
   Qc_int = CONDUCT(k_eff, T, dz, Ts, 2);

   testCase.verifyEqual(Qc_top, 2 * (270 - 269) / 0.02, 'AbsTol', 1e-12);
   testCase.verifyEqual(Qc_int, 3 * (268 - 270) / 0.08 / 2, ...
      'AbsTol', 1e-12);
end

function test_inittimesteps_and_newtimestep_follow_solver_contract(testCase)
   % Initialization and top-level timestep helpers should honor the public
   % solver contract for step size, warmup count, and flags.

   opts = struct('dt', 900, 'numyears', 2, 'n_spinup_years', 1, ...
      'simyears', [2015 2016]);
   Time = transpose(datetime(2015, 1, 1) + minutes(15) * (0:7));

   [metstep, substep, numsteps, maxsubstep, dt_new, dt_full, numyears, ...
      numspinup, simyears] = INITTIMESTEPS(opts, Time);
   [dt_sum, n_subfail, ok_seb_1, ok_ieb_1] = NEWTIMESTEP(zeros(3, 1), 1);
   [~, ~, ok_seb_2, ok_ieb_2] = NEWTIMESTEP(zeros(3, 1), 2);

   testCase.verifyEqual([metstep substep numsteps maxsubstep], [1 1 4 900]);
   testCase.verifyEqual([dt_new dt_full numyears numspinup], [900 900 2 1]);
   testCase.verifyEqual(simyears(:), [2015; 2016]);
   testCase.verifyEqual(dt_sum, 0);
   testCase.verifyEqual(n_subfail, 0);
   testCase.verifyFalse(ok_seb_1);
   testCase.verifyFalse(ok_ieb_1);
   testCase.verifyTrue(ok_seb_2);
   testCase.verifyFalse(ok_ieb_2);
end

function test_nextstep_adapts_substep_divisor(testCase)
   % NEXTSTEP should shrink or grow the substep divisor based on the recent
   % convergence history and hard failures.

   [~, substep_fast, dt_fast] = NEXTSTEP(1, 3, 300, 900, 9, true, 0, 1);
   [~, substep_slow, dt_slow] = NEXTSTEP(1, 3, 300, 900, 9, true, 2, 15);
   [~, substep_fail, dt_fail] = NEXTSTEP(1, 3, 300, 900, 9, false, 0, 0);

   testCase.verifyLessThan(substep_fast, 3);
   testCase.verifyGreaterThan(substep_slow, 3);
   testCase.verifyGreaterThan(substep_fail, 3);
   testCase.verifyGreaterThan(dt_fast, 300);
   testCase.verifyLessThan(dt_slow, 300);
   testCase.verifyLessThan(dt_fail, 300);
end

function test_resetsubstep_and_updatesubstep_restore_and_advance(testCase)
   % Reset and update helpers should restore failed-substep state, then
   % advance the accepted state and diagnostics consistently.

   [ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs] = ...
      icemodel.physicalConstant('ro_ice', 'ro_liq', 'ro_air', 'cv_ice', ...
      'cv_liq', 'roLv', 'roLs');
   [Ts, T, f_ice, f_liq, n_subfail, substep, dt_new] = RESETSUBSTEP( ...
      270, [269; 268], [0.9; 0.9], [0.01; 0.01], 900, 2, 9, 0, 450);

   [Ts_up, T_up, f_ice_up, f_liq_up, dt_sum, dt_next, liqflag, roL, ...
      ro_sno, cp_sno] = UPDATESUBSTEP(Ts, T, f_ice, f_liq, 900, 450, 300, ...
      1e-12, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, roLv, roLs);

   testCase.verifyEqual([Ts_up; T_up], [270; 269; 268], 'AbsTol', 0);
   testCase.verifyEqual([f_ice_up; f_liq_up], [0.9; 0.9; 0.01; 0.01], ...
      'AbsTol', 0);
   testCase.verifyEqual(n_subfail, 1);
   testCase.verifyEqual(substep, 3);
   testCase.verifyEqual(dt_new, 300, 'AbsTol', 1e-12);
   testCase.verifyEqual(dt_sum, 750, 'AbsTol', 1e-12);
   testCase.verifyEqual(dt_next, 150, 'AbsTol', 1e-12);
   testCase.verifyFalse(liqflag);
   testCase.verifyEqual(roL, roLs, 'AbsTol', 0);
   testCase.verifyTrue(all(isfinite(ro_sno)));
   testCase.verifyTrue(all(isfinite(cp_sno)));
end

function test_checksubstep_forces_advance_at_maxsubstep(testCase)
   % Once the max-substep limit is reached, CHECKSUBSTEP should force the
   % accepted state forward instead of stalling the timestep.

   [Ts, T, f_ice, f_liq, n_subfail, substep, dt_new, ok] = CHECKSUBSTEP( ...
      270, [269; 268], [0.9; 0.9], [0.01; 0.01], 271, [270; 269], ...
      [0.8; 0.8], [0.02; 0.02], 917, 1000, 150, 450, 900, 1, 10, 1, 2, 1, ...
      false, eps, false);

   testCase.verifyTrue(ok);
   testCase.verifyEqual(Ts, 271);
   testCase.verifyEqual(T, [270; 269]);
   testCase.verifyEqual(f_ice, [0.8; 0.8]);
   testCase.verifyEqual(f_liq, [0.02; 0.02]);
   testCase.verifyEqual(n_subfail, 2);
   testCase.verifyEqual(substep, 2);
   testCase.verifyEqual(dt_new, 450, 'AbsTol', 1e-12);
end
