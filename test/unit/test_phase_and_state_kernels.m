function tests = test_phase_and_state_kernels
   %TEST_PHASE_AND_STATE_KERNELS Verify phase-change and state-update kernels.
   tests = functiontests(localfunctions);
end

function test_liquid_fraction_function_and_derivative_are_self_consistent(testCase)
   % icemodel.column.liquid_fraction_function and
   % icemodel.column.liquid_fraction_derivative should agree when applied
   % to the same melt-zone state.

   [ro_ice, ro_liq, Tf] = icemodel.physicalConstant( ...
      'ro_ice', 'ro_liq', 'Tf');
   fcp = icemodel.parameterLookup('fcp');
   T_in = [Tf - 1.5; Tf - 0.8; Tf - 0.2];
   f_wat = [0.85; 0.85; 0.85];
   f_liq = f_wat ./ (1 + (fcp * (Tf - T_in)) .^ 2);
   f_ice = (f_wat - f_liq) * ro_liq / ro_ice;

   [T_out, f_ice_out, f_liq_out, f_wat_out, dLdT] = ...
      icemodel.column.liquid_fraction_function(T_in, f_ice, f_liq);
   [dLdT_ref, f_wat_ref] = icemodel.column.liquid_fraction_derivative( ...
      T_out, f_ice_out, f_liq_out);

   testCase.verifyEqual(T_out, T_in, 'AbsTol', 1e-10);
   testCase.verifyEqual(f_wat_out, f_wat_ref, 'AbsTol', 1e-12);
   testCase.verifyEqual(dLdT, dLdT_ref, 'RelTol', 1e-12);
end

function test_melttemp_caps_above_freezing(testCase)
   % icemodel.surface.physical_surface_temperature should clip supercooled
   % inputs at freezing but leave colder values untouched.

   Tf = icemodel.physicalConstant('Tf');
   testCase.verifyEqual( ...
      icemodel.surface.physical_surface_temperature(Tf + 3), Tf);
   testCase.verifyEqual( ...
      icemodel.surface.physical_surface_temperature(Tf - 2), Tf - 2);
end

function test_liqavail_drains_only_available_liquid(testCase)
   % available_liquid_water should respect the residual liquid floor and
   % only drain the liquid that is actually available above that floor.

   [h_resid, h_avail, h_drain, h_ice, h_liq, h_air] = ...
      icemodel.column.available_liquid_water( ...
      0.7, 0.2, 0.1, 274, 273.16, 0.02, 0.0, true, 1.0);

   testCase.verifyGreaterThan(h_resid, 0);
   testCase.verifyGreaterThan(h_avail, 0);
   testCase.verifyGreaterThan(h_drain, 0);
   testCase.verifyEqual(h_liq, h_resid, 'AbsTol', 1e-12);
   testCase.verifyEqual(h_ice + h_liq + h_air, 1.0, 'AbsTol', 1e-12);
end

function test_volbal_drains_excess_and_respects_total_volume(testCase)
   % enforce_control_volume_balance should preserve total volume while
   % routing any overflow into explicit excess terms.

   [h_ice, h_liq, h_air, x_ice, x_liq] = ...
      icemodel.column.enforce_control_volume_balance(0.8, 0.4, 0.05, 1.0);

   testCase.verifyEqual(h_ice + h_liq + h_air, 1.0, 'AbsTol', 1e-12);
   testCase.verifyGreaterThanOrEqual(x_liq, 0);
   testCase.verifyGreaterThanOrEqual(x_ice, 0);
end

function test_pevap_preserves_sign_and_scaling(testCase)
   % potential_surface_vapor_tendency should keep the latent-mass increment
   % proportional to the underlying evaporative power input.

   Qe = 50;
   dt = 900;
   dz = 0.04;

   [d_pevp, pevp] = ...
      icemodel.kernels.potential_surface_vapor_tendency(Qe, dt, dz);

   testCase.verifyGreaterThan(d_pevp, 0);
   testCase.verifyEqual(d_pevp, pevp * dt / dz, 'RelTol', 1e-12);
end

function test_mztransform_updates_melt_zone_consistently(testCase)
   % icemodel.column.meltzone_transform should move a melt-zone state
   % forward without skipping phase bounds or producing impossible phase
   % fractions.

   [ro_ice, ro_liq, Tf] = icemodel.physicalConstant('ro_ice', 'ro_liq', 'Tf');
   [TL, TH, f_ell_min, f_ell_max] = icemodel.column.meltzone_bounds();
   fcp = icemodel.parameterLookup('fcp');

   T_old = Tf - 0.5;
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - T_old)) ^ 2);
   f_liq_min = f_wat * f_ell_min;
   f_liq_max = f_wat * f_ell_max;
   dLdT = icemodel.column.liquid_fraction_derivative(T_old, [], [], f_wat);

   [T_new, f_ice_new, f_liq_new, ok] = ...
      icemodel.column.meltzone_transform(ro_liq * 0.002, T_old, f_liq, ...
      f_wat, dLdT, TL, TH, f_liq_min, f_liq_max, true, true, false);

   testCase.verifyTrue(ok);
   testCase.verifyGreaterThan(f_liq_new, f_liq);
   testCase.verifyLessThanOrEqual(f_liq_new, f_wat);
   testCase.verifyLessThanOrEqual(f_ice_new + f_liq_new * ro_liq / ro_ice, ...
      1 + 1e-9);
   testCase.verifyLessThanOrEqual(T_new, Tf);
end

function test_mztransform_allows_melt_zone_exit_to_frozen_branch(testCase)
   % A node that starts within the melt zone may legitimately freeze back
   % below TL during the corrector step without forcing a timestep retry, as
   % long as the predictor overshoots the melt-zone boundary only slightly.

   [ro_ice, ro_liq, Tf] = icemodel.physicalConstant('ro_ice', 'ro_liq', 'Tf');
   [TL, TH, f_ell_min, f_ell_max] = icemodel.column.meltzone_bounds();
   fcp = icemodel.parameterLookup('fcp');

   T_old = TL + 0.01;
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - T_old)) ^ 2);
   f_liq_min = f_wat * f_ell_min;
   f_liq_max = f_wat * f_ell_max;
   d_fliq = 0.97 * f_liq_min - f_liq;
   dLdT = icemodel.column.liquid_fraction_derivative(T_old, [], [], f_wat);

   [T_new, f_ice_new, f_liq_new, ok] = ...
      icemodel.column.meltzone_transform(ro_liq * d_fliq, T_old, f_liq, ...
      f_wat, dLdT, TL, TH, f_liq_min, f_liq_max, true, true, false);

   testCase.verifyTrue(ok);
   testCase.verifyLessThan(T_new, TL);
   testCase.verifyLessThan(f_liq_new, f_liq_min);
   testCase.verifyGreaterThan(f_liq_new, 0);
   testCase.verifyLessThanOrEqual(f_ice_new + f_liq_new * ro_liq / ro_ice, ...
      1 + 1e-9);
end

function test_mztransform_rejects_large_freeze_out_overshoot(testCase)
   % A large overshoot of the lower melt-zone boundary should be rejected so
   % the timestep can be shortened before trusting the transformed predictor.

   [ro_ice, ro_liq, Tf] = icemodel.physicalConstant('ro_ice', 'ro_liq', 'Tf');
   [TL, TH, f_ell_min, f_ell_max] = icemodel.column.meltzone_bounds();
   fcp = icemodel.parameterLookup('fcp');

   T_old = TL + 0.01;
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - T_old)) ^ 2);
   f_liq_min = f_wat * f_ell_min;
   f_liq_max = f_wat * f_ell_max;
   d_fliq = -1.5 * f_liq;
   dLdT = icemodel.column.liquid_fraction_derivative(T_old, [], [], f_wat);

   [T_new, ~, f_liq_new, ok] = icemodel.column.meltzone_transform( ...
      ro_liq * d_fliq, T_old, f_liq, f_wat, dLdT, TL, TH, f_liq_min, ...
      f_liq_max, true, true, false);

   testCase.verifyFalse(ok);
   testCase.verifyLessThan(T_new, 0);
   testCase.verifyLessThan(f_liq_new, f_liq_min);
end

function test_mztransform_rejects_phase_skip(testCase)
   % The transform should reject an attempted jump that skips across the
   % melt-zone bounds.

   [Tf] = icemodel.physicalConstant('Tf');
   [TL, TH, f_ell_min, f_ell_max] = icemodel.column.meltzone_bounds();
   fcp = icemodel.parameterLookup('fcp');
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - (TL - 1))) ^ 2);
   dLdT = 0;

   [~, ~, ~, ok] = icemodel.column.meltzone_transform(TH + 0.5, TL - 1.0, ...
      f_liq, f_wat, dLdT, TL, TH, f_wat * f_ell_min, f_wat * f_ell_max, ...
      false, true, false);

   testCase.verifyFalse(ok);
end

function test_gecoefs_applies_robin_top_boundary_adjustment(testCase)
   % icemodel.column.assemble_enthalpy_system should change the top-row
   % diagonal and source terms when the Robin boundary path is requested.

   JJ = 3;
   T = [268; 267; 266];
   f_ice = 0.9 * ones(JJ, 1);
   f_liq = 0.01 * ones(JJ, 1);
   dHdT = 1.8e6 * ones(JJ, 1);
   dLdT = 1e-3 * ones(JJ, 1);
   drovdT = 1e-6 * ones(JJ, 1);
   dH = zeros(JJ, 1);
   Sc = zeros(JJ, 1);
   k_eff = 2.0 * ones(JJ, 1);
   delz = [0.02; 0.04; 0.04; 0.02];
   fn = [1; 0.5; 0.5; 0];
   dz = 0.04 * ones(JJ, 1);
   dt = 900;
   Ts = 269;
   Fc = 10;
   Fp = -5;
   [TL, TH] = icemodel.column.meltzone_bounds();

   [~, aP_dir, ~, b_dir, ~, a1] = icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dLdT, drovdT, dH, Sc, zeros(JJ, 1), k_eff, ...
      delz, fn, dz, dt, TL, TH, Ts, Fc, Fp, 1);
   [~, aP_rob, ~, b_rob] = icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dLdT, drovdT, dH, Sc, zeros(JJ, 1), k_eff, ...
      delz, fn, dz, dt, TL, TH, Ts, Fc, Fp, 2);

   testCase.verifyEqual(aP_rob(1) - aP_dir(1), ...
      -a1 - Fp * a1 / (a1 - Fp), ...
      'RelTol', 1e-12);
   testCase.verifyNotEqual(b_rob(1), b_dir(1));
end

function test_iceablation_and_surface_runoff_budget(testCase)
   % Surface ablation helpers should report melt, sublimation, and runoff
   % in a self-consistent budget on a simple forcing state.

   opts = struct('smbmodel', 'skinmodel', 'skinfreeze', true);
   [surf_mlt, surf_frz, surf_sub, surf_con, surf_rof] = ...
      icemodel.surface.diagnose_surface_ablation(50, -20, 10, 0, 0, 0, ...
      0, 0, 3600, opts);

   testCase.verifyGreaterThan(surf_mlt, 0);
   testCase.verifyEqual(surf_frz, 0, 'AbsTol', 1e-12);
   testCase.verifyGreaterThan(surf_sub, 0);
   testCase.verifyEqual(surf_rof, surf_mlt, 'RelTol', 1e-12);
   testCase.verifyEqual(surf_con, 0, 'AbsTol', 1e-12);
end

function test_surface_runoff_and_ice_runoff_build_cumulative_series(testCase)
   % The runoff accumulators should produce monotonic cumulative series from
   % incremental latent-mass inputs.

   ice1 = struct('Qe', [-10; 5; 0], 'Qm', [20; -5; 10]);
   ice1 = icemodel.surface.diagnose_surface_runoff(ice1, 3600);
   testCase.verifyGreaterThan(ice1.runoff(end), 0);
   testCase.verifyGreaterThan(ice1.melt(end), 0);

   opts = struct('dz_thermal', 0.04, 'tlag', 1);
   ice2 = struct('df_liq', [0.1 -0.05 0.0; 0.0 0.05 -0.02]);
   ice1_run = icemodel.column.diagnose_column_runoff(struct(), ice2, opts);
   testCase.verifyGreaterThanOrEqual(min(diff(ice1_run.runoff)), 0);
   testCase.verifyGreaterThan(ice1_run.melt(end), 0);
end

function test_budget_surface_mass_balance_and_merge_thin_layers(testCase)
   % budget_surface_mass_balance followed by merge_thin_layers should
   % combine a too-thin surface layer while preserving the expected output
   % array sizes and merge diagnostic.

   Tf = icemodel.physicalConstant('Tf');

   T = [Tf - 2; Tf - 2.5; Tf - 3.0];
   f_ice = [0.05; 0.60; 0.65];
   f_liq = [0.01; 0.01; 0.01];
   xf_liq = f_liq;
   Sc = zeros(3, 1);
   Sp = zeros(3, 1);
   d_liq = zeros(3, 1);
   d_evp = zeros(3, 1);
   d_lyr = zeros(3, 1);

   [T_new, f_ice_new, f_liq_new, ~, ~] = ...
      icemodel.column.budget_surface_mass_balance(T, f_ice, f_liq, ...
      xf_liq, 0.0, d_liq, d_evp, 0.02, 0.1);
   [T_new, f_ice_new, f_liq_new, ~, ~, d_lyr_new, lcflag] = ...
      icemodel.column.merge_thin_layers(T_new, f_ice_new, f_liq_new, Sc, ...
      Sp, 0.04, 0.0, d_lyr, 0.1);

   testCase.verifyTrue(any(lcflag));
   testCase.verifyEqual(numel(T_new), 3);
   testCase.verifyEqual(numel(f_ice_new), 3);
   testCase.verifyEqual(numel(f_liq_new), 3);
   testCase.verifyGreaterThan(sum(d_lyr_new), 0);
end
