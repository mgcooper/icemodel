function tests = test_phase_and_state_kernels
   %TEST_PHASE_AND_STATE_KERNELS Verify phase-change and state-update kernels.
   tests = functiontests(localfunctions);
end

function test_liquid_fraction_function_and_deriv_are_self_consistent(testCase)
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
   % potential_surface_vapor_demand should keep the latent-mass increment
   % proportional to the underlying evaporative power input.

   Qe = 50;
   dt = 900;
   dz = 0.04;

   [d_pevp, pevp] = ...
      icemodel.kernels.potential_surface_vapor_demand(Qe, dt, dz);

   testCase.verifyGreaterThan(d_pevp, 0);
   testCase.verifyEqual(d_pevp, pevp * dt / dz, 'RelTol', 1e-12);
end

function test_mztransform_updates_melt_zone_consistently(testCase)
   % icemodel.column.meltzone_transform should move a melt-zone state
   % forward without skipping phase bounds or producing impossible phase
   % fractions.

   [ro_ice, ro_liq, Tf] = icemodel.physicalConstant('ro_ice', 'ro_liq', 'Tf');
   fcp = icemodel.parameterLookup('fcp');

   T_old = Tf - 0.5;
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - T_old)) ^ 2);
   [f_liq_min, f_liq_max] = icemodel.column.meltzone_bounds(f_wat);
   dLdT = icemodel.column.liquid_fraction_derivative(T_old, [], [], f_wat);

   [T_new, f_ice_new, f_liq_new, ok] = ...
      icemodel.column.meltzone_transform(ro_liq * 0.002, T_old, f_liq, ...
      f_wat, dLdT, f_liq_min, f_liq_max, true, false);

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
   [TL, ~] = icemodel.parameterLookup('TL', 'TH');
   fcp = icemodel.parameterLookup('fcp');

   T_old = TL + 0.01;
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - T_old)) ^ 2);
   [f_liq_min, f_liq_max] = icemodel.column.meltzone_bounds(f_wat);
   d_fliq = 0.97 * f_liq_min - f_liq;
   dLdT = icemodel.column.liquid_fraction_derivative(T_old, [], [], f_wat);

   [T_new, f_ice_new, f_liq_new, ok] = ...
      icemodel.column.meltzone_transform(ro_liq * d_fliq, T_old, f_liq, ...
      f_wat, dLdT, f_liq_min, f_liq_max, true, false);

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

   [~, ro_liq, Tf] = icemodel.physicalConstant('ro_ice', 'ro_liq', 'Tf');
   [TL, ~] = icemodel.parameterLookup('TL', 'TH');
   fcp = icemodel.parameterLookup('fcp');

   T_old = TL + 0.01;
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - T_old)) ^ 2);
   [f_liq_min, f_liq_max] = icemodel.column.meltzone_bounds(f_wat);
   d_fliq = -1.5 * f_liq;
   dLdT = icemodel.column.liquid_fraction_derivative(T_old, [], [], f_wat);

   [T_new, ~, f_liq_new, ok] = icemodel.column.meltzone_transform( ...
      ro_liq * d_fliq, T_old, f_liq, f_wat, dLdT, f_liq_min, ...
      f_liq_max, true, false);

   testCase.verifyFalse(ok);
   testCase.verifyLessThan(T_new, 0);
   testCase.verifyLessThan(f_liq_new, f_liq_min);
end

function test_mztransform_rejects_phase_skip(testCase)
   % The transform should reject an attempted jump that skips across the
   % melt-zone bounds.

   [Tf] = icemodel.physicalConstant('Tf');
   [TL, TH] = icemodel.parameterLookup('TL', 'TH');
   fcp = icemodel.parameterLookup('fcp');
   f_wat = 0.85;
   f_liq = f_wat / (1 + (fcp * (Tf - (TL - 1))) ^ 2);
   dLdT = 0;
   [f_liq_min, f_liq_max] = icemodel.column.meltzone_bounds(f_wat);

   [~, ~, ~, ok] = icemodel.column.meltzone_transform(TH + 0.5, TL - 1.0, ...
      f_liq, f_wat, dLdT, f_liq_min, f_liq_max, ...
      false, false);

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
   k_eff_faces = 1.0 ./ ((1.0 - fn) ./ [k_eff(1); k_eff] ...
      + fn ./ [k_eff; k_eff(end)]);
   q_deferred_faces = zeros(JJ + 1, 1);
   [~, aP_dir, ~, b_dir, ~, a1] = icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dLdT, drovdT, dH, Sc, zeros(JJ, 1), ...
      k_eff_faces, delz, dz, dt, Ts, Fc, Fp, 1, q_deferred_faces);
   [~, aP_rob, ~, b_rob] = icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dLdT, drovdT, dH, Sc, zeros(JJ, 1), ...
      k_eff_faces, delz, dz, dt, Ts, Fc, Fp, 2, q_deferred_faces);

   testCase.verifyEqual(aP_rob(1) - aP_dir(1), ...
      -a1 - Fp * a1 / (a1 - Fp), ...
      'RelTol', 1e-12);
   testCase.verifyNotEqual(b_rob(1), b_dir(1));
end

function test_assemble_enthalpy_system_adds_deferred_flux_convergence(testCase)
   % The deferred face flux must enter each cell as in minus out.

   JJ = 3;
   T = [268; 267; 266];
   f_ice = 0.8 * ones(JJ, 1);
   f_liq = zeros(JJ, 1);
   dHdT = 1.8e6 * ones(JJ, 1);
   dFdT = zeros(JJ, 1);
   drovdT = 1e-6 * ones(JJ, 1);
   dH = zeros(JJ, 1);
   Sc = zeros(JJ, 1);
   Sp = zeros(JJ, 1);
   k_eff_faces = 2.0 * ones(JJ + 1, 1);
   delz = [0.02; 0.04; 0.04; 0.02];
   dz = 0.04 * ones(JJ, 1);
   q_zero = zeros(JJ + 1, 1);
   q_deferred_faces = [0; 3; -2; 0];

   [~, ~, ~, b_zero] = icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dFdT, drovdT, dH, Sc, Sp, ...
      k_eff_faces, delz, dz, 900, 269, 0, 0, 1, q_zero);
   [~, ~, ~, b_flux] = icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dFdT, drovdT, dH, Sc, Sp, ...
      k_eff_faces, delz, dz, 900, 269, 0, 0, 1, q_deferred_faces);

   expected = q_deferred_faces(1:JJ) - q_deferred_faces(2:JJ+1);
   testCase.verifyEqual(b_flux - b_zero, expected, 'AbsTol', 10 * eps);
end

function test_subsurface_error_matches_zero_flux_and_adds_top_convergence(testCase)
   % Zero deferred flux must reproduce the matrix-and-source formula.
   % Nonzero flux adds the assembler's top-face convergence.

   [Lf, ro_liq] = icemodel.physicalConstant('Lf', 'ro_liq');
   T = [268.2; 267.8; 267.0];
   T_old = [268.0; 267.7; 267.0];
   f_ice = 0.8 * ones(3, 1);
   f_liq = zeros(3, 1);
   f_liq_old = zeros(3, 1);
   drovdT = 1e-6 * ones(3, 1);
   dHdT = 1.8e6 * ones(3, 1);
   Sc = [2; 0; 0];
   dz = 0.04 * ones(3, 1);
   dt = 900;
   a2 = 25;
   Fc = 10;
   Fp = -5;
   a1 = 40;
   q_zero = zeros(4, 1);

   returned_zero = icemodel.column.subsurface_linearization_error( ...
      T, T_old, f_ice, f_liq, f_liq_old, drovdT, dHdT, Sc, q_zero, ...
      dz, dt, a2, Fc, Fp, a1);

   % Assemble the zero-deferred-flux reference independently.
   L_top = icemodel.vapor.latent_enthalpy_switch(f_liq(1));
   denominator = dHdT(1) ...
      + L_top * drovdT(1) * (1.0 - f_ice(1) - f_liq(1));
   zero_flux_reference = (dt / dz(1) ...
      * (a2 * (T(2) - T(1)) ...
      + Fc + Fp * (Fc + a1 * T(1)) / (a1 - Fp) ...
      + Sc(1) * dz(1)) ...
      - ro_liq * Lf * (f_liq(1) - f_liq_old(1))) / denominator ...
      - (T(1) - T_old(1));
   testCase.verifyEqual( ...
      returned_zero, zero_flux_reference, 'RelTol', 1e-14);

   % Build the nonzero case through the production assembler. The state is
   % below the melt zone, so the assembled unknown is T and the top-row
   % residual can be evaluated directly from aN, aP, aS, and b.
   q_deferred_faces = [4; -3; 0; 0];
   dFdT = zeros(3, 1);
   dH = denominator * (T - T_old) ...
      + ro_liq * Lf * (f_liq - f_liq_old);
   k_eff_faces = [0.9; 0.8; 0.7; 0.6];
   delz = [0.02; 0.04; 0.04; 0.02];
   T_sfc = 269;
   bc = 2;
   [aN, aP, aS, b, ~, a1_assembled, a2_assembled] = ...
      icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dFdT, drovdT, dH, Sc, zeros(3, 1), ...
      k_eff_faces, delz, dz, dt, T_sfc, Fc, Fp, bc, ...
      q_deferred_faces);
   returned_flux = icemodel.column.subsurface_linearization_error( ...
      T, T_old, f_ice, f_liq, f_liq_old, drovdT, dHdT, Sc, ...
      q_deferred_faces, dz, dt, a2_assembled, Fc, Fp, a1_assembled);
   top_row_residual = b(1) + aN(1) * T_sfc + aS(1) * T(2) ...
      - aP(1) * T(1);
   assembled_reference = top_row_residual * dt / dz(1) / denominator;
   testCase.verifyEqual(returned_flux, assembled_reference, ...
      'RelTol', 1e-13);
end

function test_subsurface_error_supports_one_closed_cell(testCase)
   % A one-cell column has no interior south flux, but its two boundary-face
   % deferred terms must still enter as in-minus-out convergence.

   [Lf, ro_liq] = icemodel.physicalConstant('Lf', 'ro_liq');
   T = 268.2;
   T_old = 268.0;
   f_ice = 0.8;
   f_liq = 0.0;
   f_liq_old = 0.0;
   drovdT = 1e-6;
   dHdT = 1.8e6;
   Sc = 2;
   dz = 0.04;
   dt = 900;
   a2 = 25;
   Fc = 10;
   Fp = -5;
   a1 = 40;

   L_top = icemodel.vapor.latent_enthalpy_switch(f_liq);
   denominator = dHdT ...
      + L_top * drovdT * (1.0 - f_ice - f_liq);
   reference = (dt / dz ...
      * (Fc + Fp * (Fc + a1 * T) / (a1 - Fp) + Sc * dz) ...
      - ro_liq * Lf * (f_liq - f_liq_old)) / denominator ...
      - (T - T_old);

   returned_zero = icemodel.column.subsurface_linearization_error( ...
      T, T_old, f_ice, f_liq, f_liq_old, drovdT, dHdT, Sc, zeros(2, 1), ...
      dz, dt, a2, Fc, Fp, a1);
   testCase.verifyEqual(returned_zero, reference, 'RelTol', 1e-14);

   q_deferred_faces = [4; -3];
   returned_flux = icemodel.column.subsurface_linearization_error( ...
      T, T_old, f_ice, f_liq, f_liq_old, drovdT, dHdT, Sc, ...
      q_deferred_faces, dz, dt, a2, Fc, Fp, a1);
   expected_flux = reference + dt / dz ...
      * (q_deferred_faces(1) - q_deferred_faces(2)) / denominator;
   testCase.verifyEqual(returned_flux, expected_flux, 'RelTol', 1e-14);
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
   % array sizes and remesh budget counts.

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
   dz = 0.04;

   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [T_new, f_ice_new, f_liq_new, ~, ~, ~, ~, ~, ~, budget] = ...
      icemodel.column.budget_surface_mass_balance(T, f_ice, f_liq, ...
      xf_liq, 0.0, d_liq, d_evp, 0.0, zeros(3, 1), zeros(3, 1), 0.02, ...
      0.1, budget, dz);
   [T_new, f_ice_new, f_liq_new, ~, ~, d_lyr_new, budget] = ...
      icemodel.column.merge_thin_layers(T_new, f_ice_new, f_liq_new, Sc, ...
      Sp, dz, 0.0, d_lyr, 0.1, budget);

   % f_ice(1) = 0.05 starts below f_ice_min = 0.1, so the top cell is the one
   % merged; the observable is the budget's top-removal count, because
   % merge_thin_layers returns no per-cell eligibility mask.
   testCase.verifyGreaterThan(budget.mass_budget_top_deletion_count, 0);
   testCase.verifyEqual(numel(T_new), 3);
   testCase.verifyEqual(numel(f_ice_new), 3);
   testCase.verifyEqual(numel(f_liq_new), 3);
   testCase.verifyGreaterThan(sum(d_lyr_new), 0);
end
