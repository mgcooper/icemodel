function tests = test_update_grain_radius
   %TEST_UPDATE_GRAIN_RADIUS Pin realized-exchange grain growth.
   %
   % See also: icemodel.column.update_grain_radius
   tests = functiontests(localfunctions);
end

function test_dry_growth_uses_gross_realized_exchange(testCase)
   % Dry growth reproduces Jordan Eq. 33 from gross face throughput.

   [r_eff, f_liq, d_vap_faces, dt] = grainFixture(4);
   [g1, r_max, Uv_max] = icemodel.parameterLookup('g1', 'r_max', 'Uv_max');
   ro_liq = icemodel.physicalConstant('ro_liq');

   U_vap_faces = d_vap_faces * ro_liq / dt;
   U_vap_nodes = min(0.5 * (U_vap_faces(1:4) + U_vap_faces(2:5)), Uv_max);
   diam = 2 * r_eff;
   expected = min(0.5 * (diam + dt * g1 .* U_vap_nodes ./ diam), r_max);

   returned = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, d_vap_faces, dt);
   testCase.verifyEqual(returned, expected, 'AbsTol', 0);

   % No realized exchange leaves dry grains unchanged.
   still = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, zeros(5, 1), dt);
   testCase.verifyEqual(still, r_eff, 'AbsTol', 0);
end

function test_dry_growth_caps_the_step_mean_flux(testCase)
   % The Uv_max cap limits large accumulated exchanges at every node.

   [r_eff, f_liq, ~, dt] = grainFixture(3);
   [g1, r_max, Uv_max] = icemodel.parameterLookup('g1', 'r_max', 'Uv_max');
   ro_liq = icemodel.physicalConstant('ro_liq');
   d_vap_faces = 2 * Uv_max * dt / ro_liq * ones(4, 1);

   diam = 2 * r_eff;
   expected = min(0.5 * (diam + dt * g1 * Uv_max ./ diam), r_max);
   returned = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, d_vap_faces, dt);

   testCase.verifyEqual(returned, expected, 'AbsTol', 0);
end

function test_wet_growth_uses_both_liquid_branches(testCase)
   % These thresholds select dry, liquid-dependent, and capped wet growth.

   [r_eff, ~, ~, dt] = grainFixture(6);
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
      r_eff, f_liq, zeros(7, 1), dt);
   testCase.verifyEqual(returned, expected, 'AbsTol', 0);
   testCase.verifyEqual(returned(5), returned(6), 'AbsTol', 0);
end

function test_growth_stops_at_the_radius_cap(testCase)
   % A long wet step reaches r_max, while a normal step remains below it.

   r_max = icemodel.parameterLookup('r_max');
   r_eff = 0.99 * r_max * ones(3, 1);
   f_liq = 0.30 * ones(3, 1);
   d_vap_faces = zeros(4, 1);

   capped = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, d_vap_faces, 8.64e5);
   testCase.verifyEqual(capped, r_max * ones(3, 1), 'AbsTol', 0);

   uncapped = icemodel.column.update_grain_radius( ...
      r_eff, f_liq, d_vap_faces, 900);
   testCase.verifyGreaterThan(uncapped, r_eff);
   testCase.verifyLessThan(uncapped, r_max * ones(3, 1));
end

function test_gross_reversals_add_instead_of_cancel(testCase)
   % Opposite substep exchanges must contribute through their gross sum.

   [r_eff, f_liq, ~, dt] = grainFixture(3);
   net = zeros(4, 1);
   gross = zeros(4, 1);
   gross(1) = 4e-7;

   from_net = icemodel.column.update_grain_radius(r_eff, f_liq, net, dt);
   from_gross = icemodel.column.update_grain_radius(r_eff, f_liq, gross, dt);

   testCase.verifyEqual(from_net, r_eff, 'AbsTol', 0);
   testCase.verifyGreaterThan(from_gross(1), r_eff(1));
end

function test_interface_has_no_saturation_dependency(testCase)
   % Grain growth accepts only realized exchange and retains codegen support.

   testCase.verifyEqual(nargin(@icemodel.column.update_grain_radius), 4);
   testCase.verifyEqual(nargout(@icemodel.column.update_grain_radius), 1);

   source = fileread(which('icemodel.column.update_grain_radius'));
   testCase.verifyFalse(contains(source, 'icemodel.vapor.'));
   testCase.verifyTrue(contains(source, '%#codegen'));
end

function [r_eff, f_liq, d_vap_faces, dt] = grainFixture(JJ)
   %GRAINFIXTURE Return one dry column and gross face exchange.

   r_eff = 5e-4 * ones(JJ, 1);
   f_liq = zeros(JJ, 1);
   d_vap_faces = linspace(5e-7, 0, JJ + 1)';
   dt = 900;
end
