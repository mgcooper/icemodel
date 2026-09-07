function tests = test_update_grain_radius
   %TEST_UPDATE_GRAIN_RADIUS Pin realized-exchange grain growth.
   %
   % See also: icemodel.column.update_grain_radius
   tests = functiontests(localfunctions);
end

function test_dry_growth_uses_realized_substep_flux(testCase)
   % Dry growth reproduces Jordan Eq. 33 from the substep's own instantaneous
   % face flux: interior faces from U_vap, and face 1 rebuilt from the
   % realized top-cell exchange d_vap.

   [r_eff, f_liq, U_vap, d_vap, dz1, dt] = grainFixture(4);
   [g1, r_max, Uv_max] = icemodel.parameterLookup('g1', 'r_max', 'Uv_max');
   ro_liq = icemodel.physicalConstant('ro_liq');
   JJ = numel(r_eff);

   U_vap_faces = U_vap;
   U_vap_faces(1) = abs(d_vap) * dz1 * ro_liq / dt;
   U_vap_nodes = min(0.5 * (abs(U_vap_faces(1:JJ)) ...
      + abs(U_vap_faces(2:JJ + 1))), Uv_max);
   diam = 2 * r_eff;
   expected = min(0.5 * (diam + dt * g1 .* U_vap_nodes ./ diam), r_max);

   returned = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, U_vap, d_vap, dz1, dt);
   testCase.verifyEqual(returned, expected, 'AbsTol', 0);

   % No realized exchange and no interior flux leaves dry grains unchanged.
   still = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, zeros(JJ + 1, 1), 0.0, dz1, dt);
   testCase.verifyEqual(still, r_eff, 'AbsTol', 0);
end

function test_dry_growth_caps_the_step_mean_flux(testCase)
   % The Uv_max cap limits large accumulated exchanges at every node.

   [r_eff, f_liq, ~, ~, dz1, dt] = grainFixture(3);
   [g1, r_max, Uv_max] = icemodel.parameterLookup('g1', 'r_max', 'Uv_max');
   ro_liq = icemodel.physicalConstant('ro_liq');
   JJ = numel(r_eff);

   % Drive every interior face, and the surface face through d_vap,
   % to twice Uv_max so every node's mean magnitude clamps at the cap.
   U_vap = 2 * Uv_max * ones(JJ + 1, 1);
   d_vap = 2 * Uv_max * dt / (dz1 * ro_liq);

   diam = 2 * r_eff;
   expected = min(0.5 * (diam + dt * g1 * Uv_max ./ diam), r_max);
   returned = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, U_vap, d_vap, dz1, dt);

   testCase.verifyEqual(returned, expected, 'AbsTol', 0);
end

function test_wet_growth_uses_both_liquid_branches(testCase)
   % These thresholds select dry, liquid-dependent, and capped wet growth.

   [r_eff, ~, ~, ~, dz1, dt] = grainFixture(6);
   f_liq = [9e-5; 1e-4; 0.02; 0.089; 0.09; 0.30];
   [g2, r_max] = icemodel.parameterLookup('g2', 'r_max');
   diam = 2 * r_eff;

   expected_diam = diam;
   wet_lo = 2:4;
   expected_diam(wet_lo) = expected_diam(wet_lo) + ...
      dt * g2 .* (f_liq(wet_lo) + 0.05) ./ diam(wet_lo);
   wet_hi = 5:6;
   expected_diam(wet_hi) = expected_diam(wet_hi) + ...
      dt * g2 * 0.14 ./ diam(wet_hi);
   expected = min(0.5 * expected_diam, r_max);

   returned = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, zeros(7, 1), 0.0, dz1, dt);
   testCase.verifyEqual(returned, expected, 'AbsTol', 0);
   testCase.verifyEqual(returned(5), returned(6), 'AbsTol', 0);
end

function test_growth_stops_at_the_radius_cap(testCase)
   % A long wet step reaches r_max, while a normal step remains below it.

   r_max = icemodel.parameterLookup('r_max');
   r_eff = 0.99 * r_max * ones(3, 1);
   f_liq = 0.30 * ones(3, 1);
   U_vap = zeros(4, 1);
   dz1 = 0.04;

   capped = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, U_vap, 0.0, dz1, 8.64e5);
   testCase.verifyEqual(capped, r_max * ones(3, 1), 'AbsTol', 0);

   uncapped = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, U_vap, 0.0, dz1, 900);
   testCase.verifyGreaterThan(uncapped, r_eff);
   testCase.verifyLessThan(uncapped, r_max * ones(3, 1));
end

function test_opposite_signed_substeps_add_instead_of_cancel(testCase)
   % update_grain_radius runs once per accepted substep on that substep's own
   % flux magnitude, not on a cross-substep signed accumulator. A sign
   % reversal between two chained substeps must therefore not cancel the
   % accumulated growth.

   [r_eff, f_liq, U_vap, ~, dz1, dt] = grainFixture(3);
   d_vap_up = 6e-4;
   d_vap_down = -6e-4;

   grown_up = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, U_vap, d_vap_up, dz1, dt);
   grown_down = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, U_vap, d_vap_down, dz1, dt);
   grown_zero = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, U_vap, 0.0, dz1, dt);

   % Equal-magnitude opposite-signed exchanges cause identical growth,
   % because the flux enters through abs(d_vap).
   testCase.verifyEqual(grown_up, grown_down, 'AbsTol', 0);
   testCase.verifyGreaterThan(grown_up(1), grown_zero(1));

   % Chaining the up substep and then the down substep must accumulate
   % growth rather than net back toward the zero-exchange result a signed
   % accumulator would produce.
   chained = icemodel.column.update_grain_radius( ...
      grown_up, f_liq, U_vap, d_vap_down, dz1, dt);
   testCase.verifyGreaterThan(chained(1), grown_up(1));
end

function test_interface_has_no_saturation_dependency(testCase)
   % Grain growth accepts only realized exchange and retains codegen support.

   testCase.verifyEqual(nargin(@icemodel.column.update_grain_radius), 6);
   testCase.verifyEqual(nargout(@icemodel.column.update_grain_radius), 1);

   source = fileread(which('icemodel.column.update_grain_radius'));
   testCase.verifyFalse(contains(source, 'icemodel.vapor.'));
   testCase.verifyTrue(contains(source, '%#codegen'));
end

function [r_eff, f_liq, U_vap, d_vap, dz1, dt] = grainFixture(JJ)
   %GRAINFIXTURE Return one dry column and an interior-face vapor flux.

   r_eff = 5e-4 * ones(JJ, 1);
   f_liq = zeros(JJ, 1);

   % Interior faces (2:JJ) carry a nonzero mass flux [kg m-2 s-1]. The
   % boundary faces (1 and JJ+1) are zero on input: update_grain_radius
   % overrides face 1 from d_vap, and the interior call always closes
   % face JJ+1.
   U_vap = zeros(JJ + 1, 1);
   U_vap(2:JJ) = linspace(5e-7, 1e-7, JJ - 1)';
   d_vap = 3e-4;
   dz1 = 0.04;
   dt = 900;
end
