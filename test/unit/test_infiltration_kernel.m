function tests = test_infiltration_kernel
   tests = functiontests(localfunctions);
end

function test_zero_input_zero_output(testCase)
   % With q_top = 0, ripe column at Tf, no flux, no state change.
   N = 20;
   dz = 0.05;
   dt = 60;
   Tf = icemodel.physicalConstant('Tf');
   f_ice = 0.30 * ones(N, 1);
   f_liq = 0.05 * ones(N, 1);  % above residual but below sat
   T = Tf * ones(N, 1);
   q_top = 0;
   [f_liq2, f_ice2, T2, diag] = icemodel.column.infiltration( ...
      f_liq, f_ice, T, dz, dt, q_top);
   verifyEqual(testCase, diag.inflow_total, 0, 'AbsTol', 0);
   % Ripe column: liquid will drain by gravity even with q_top=0.
   verifyTrue(testCase, all(f_ice2 == f_ice), 'ice unchanged at Tf');
   verifyTrue(testCase, all(T2 <= Tf + eps), 'no superheating');
   % Mass balance: inflow + d_storage = outflow (no refreezing at Tf).
   d_liq = sum(f_liq2 - f_liq) * dz;
   verifyEqual(testCase, -d_liq, diag.outflow_total, 'AbsTol', 1e-12);
end

function test_mass_balance_ripe_column(testCase)
   % Inflow + storage_change_liq = outflow + storage_change_ice (zero on ripe).
   N = 20;
   dz = 0.05;
   dt = 60;
   Tf = icemodel.physicalConstant('Tf');
   f_ice = 0.30 * ones(N, 1);
   f_liq = 0.05 * ones(N, 1);
   T = Tf * ones(N, 1);
   q_top = 1e-6;
   [f_liq2, f_ice2, T2, diag] = icemodel.column.infiltration( ...
      f_liq, f_ice, T, dz, dt, q_top); %#ok<ASGLU>
   d_liq = sum(f_liq2 - f_liq) * dz;
   residual = (diag.inflow_total - diag.outflow_total) - d_liq;
   verifyEqual(testCase, residual, 0, 'AbsTol', 1e-9);
   % No refreezing at Tf -> ice unchanged.
   verifyEqual(testCase, f_ice2, f_ice, 'AbsTol', 1e-15);
end

function test_cold_content_refreezing(testCase)
   % Inflow into a sub-freezing column refreezes -> f_liq drops, T rises,
   % f_ice grows; total water-equivalent mass conserved.
   N = 20;
   dz = 0.05;
   dt = 60;
   [Tf, ro_ice, ro_liq] = icemodel.physicalConstant('Tf', 'ro_ice', 'ro_liq');
   f_ice = 0.30 * ones(N, 1);
   f_liq = zeros(N, 1);
   T = (Tf - 5) * ones(N, 1);
   q_top = 1e-6;
   [f_liq2, f_ice2, T2, diag] = icemodel.column.infiltration( ...
      f_liq, f_ice, T, dz, dt, q_top);
   verifyTrue(testCase, sum(f_ice2) > sum(f_ice), 'ice should grow');
   verifyTrue(testCase, mean(T2) > mean(T), 'column should warm');

   % Total water-equivalent mass balance.
   ro_iwe = ro_ice / ro_liq;
   we_init = (sum(f_liq) + sum(f_ice) * ro_iwe) * dz;
   we_final = (sum(f_liq2) + sum(f_ice2) * ro_iwe) * dz;
   residual = (we_final - we_init) - (diag.inflow_total - diag.outflow_total);
   verifyEqual(testCase, residual, 0, 'AbsTol', 1e-9);
end

function test_overfilled_layer_clamps_to_the_water_equivalent_bound(testCase)
   % The ioverfill clamp had no test. A layer that starts above the bound must
   % come back to exactly ro_ice/ro_liq * (1 - f_ice). That is the pore volume
   % in water equivalent, and the clamp goes no further.

   N = 5;
   dz = 0.05;
   dt = 60;
   Tf = icemodel.physicalConstant('Tf');
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');

   f_ice = 0.30 * ones(N, 1);
   T = Tf * ones(N, 1);

   % Start cell 3 above the bound, so the clamp must pull it back.
   bound = ro_ice / ro_liq * (1 - f_ice);
   f_liq = 0.05 * ones(N, 1);
   f_liq(3) = bound(3) + 0.02;
   testCase.assertGreaterThan(f_liq(3), bound(3));

   % A short step so drainage cannot move the overfilled layer before the
   % clamp reaches it. Over a long step the layer drains instead, and the
   % final value comes from the flux, not the clamp.
   [f_liq_short, ~, ~, ~] = icemodel.column.infiltration( ...
      f_liq, f_ice, T, dz, 1e-3, 0);

   % The clamped layer lands on the bound exactly, not below it.
   testCase.verifyEqual(f_liq_short(3), bound(3), 'AbsTol', 1e-15);

   % Over a realistic step the invariant still holds: no layer sits above
   % the bound, whether the clamp or the drainage got it there.
   [f_liq2, ~, ~, ~] = icemodel.column.infiltration( ...
      f_liq, f_ice, T, dz, dt, 0);
   testCase.verifyLessThanOrEqual(f_liq2, bound + 1e-15);

   % The bound is the pore volume scaled to water equivalent. It sits below
   % the pore volume by (1 - ro_ice/ro_liq) * (1 - f_ice).
   shortfall = (1 - f_ice(3)) - bound(3);
   testCase.verifyEqual(shortfall, ...
      (1 - ro_ice / ro_liq) * (1 - f_ice(3)), 'RelTol', 1e-14);
end

function test_n_sub_increases_with_q_top(testCase)
   % Larger top flux -> larger characteristic speed -> more substeps.
   N = 10;
   dz = 0.05;
   dt = 60;
   Tf = icemodel.physicalConstant('Tf');
   f_ice = 0.30 * ones(N, 1);
   f_liq = 0.05 * ones(N, 1);
   T = Tf * ones(N, 1);
   [~, ~, ~, d_small] = icemodel.column.infiltration( ...
      f_liq, f_ice, T, dz, dt, 1e-7);
   [~, ~, ~, d_large] = icemodel.column.infiltration( ...
      f_liq, f_ice, T, dz, dt, 1e-5);
   verifyGreaterThanOrEqual(testCase, d_large.n_sub, d_small.n_sub);
end
