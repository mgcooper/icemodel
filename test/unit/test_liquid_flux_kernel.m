function tests = test_liquid_flux_kernel
   tests = functiontests(localfunctions);
end

function test_zero_flux_at_or_below_residual(testCase)
   % Below the capillary-residual threshold, no liquid drains.
   f_ice = 0.30 * ones(5, 1);
   f_res_pore = 0.07;
   f_res = f_res_pore * (1 - f_ice);
   f_liq = f_res - 1e-6;
   [q, dq] = icemodel.column.liquid_flux(f_liq, f_ice, f_res_pore=f_res_pore);
   verifyEqual(testCase, q, zeros(5, 1), 'AbsTol', 1e-15);
   verifyEqual(testCase, dq, zeros(5, 1), 'AbsTol', 1e-15);
end

function test_monotone_in_liquid(testCase)
   % q is monotone non-decreasing in f_liq above the residual.
   f_ice = 0.30;
   f_res_pore = 0.07;
   f_res = f_res_pore * (1 - f_ice);
   f_liq_grid = linspace(f_res + 1e-6, 0.5, 25)';
   q = zeros(size(f_liq_grid));
   for i = 1:numel(f_liq_grid)
      q(i) = icemodel.column.liquid_flux(f_liq_grid(i), f_ice, f_res_pore=f_res_pore);
   end
   verifyTrue(testCase, all(diff(q) > 0), ...
      'flux must increase strictly with f_liq above residual');
end

function test_dq_dfliq_matches_finite_difference_colbeck(testCase)
   % Closed-form derivative agrees with finite-difference for Colbeck path.
   f_ice = 0.30;
   f_res_pore = 0.07;
   f_liq0 = 0.20;
   h = 1e-7;
   [~, dq_analytic] = icemodel.column.liquid_flux(f_liq0, f_ice, f_res_pore=f_res_pore);
   q_plus  = icemodel.column.liquid_flux(f_liq0 + h, f_ice, f_res_pore=f_res_pore);
   q_minus = icemodel.column.liquid_flux(f_liq0 - h, f_ice, f_res_pore=f_res_pore);
   dq_fd = (q_plus - q_minus) / (2 * h);
   verifyEqual(testCase, dq_analytic, dq_fd, 'RelTol', 1e-4);
end

function test_dq_dfliq_matches_finite_difference_shimizu(testCase)
   % Closed-form derivative agrees with finite-difference for Shimizu /
   % van Genuchten path.
   f_ice = 0.30;
   f_res_pore = 0.07;
   grainsz = 2e-3;
   f_liq0 = 0.20;
   h = 1e-7;
   [~, dq_analytic] = icemodel.column.liquid_flux(f_liq0, f_ice, ...
      f_res_pore=f_res_pore, k_sat_method="shimizu1970", grainsz=grainsz);
   q_plus  = icemodel.column.liquid_flux(f_liq0 + h, f_ice, ...
      f_res_pore=f_res_pore, k_sat_method="shimizu1970", grainsz=grainsz);
   q_minus = icemodel.column.liquid_flux(f_liq0 - h, f_ice, ...
      f_res_pore=f_res_pore, k_sat_method="shimizu1970", grainsz=grainsz);
   dq_fd = (q_plus - q_minus) / (2 * h);
   verifyEqual(testCase, dq_analytic, dq_fd, 'RelTol', 1e-4);
end

function test_no_flux_when_pore_space_zero(testCase)
   % Solid ice (no pore space) drains nothing regardless of f_liq input.
   f_ice = 1.0 * ones(3, 1);
   f_liq = 0.0 * ones(3, 1);
   q = icemodel.column.liquid_flux(f_liq, f_ice, f_res_pore=0.07);
   verifyEqual(testCase, q, zeros(3, 1), 'AbsTol', 1e-15);
end

function test_q_units_are_m_per_s(testCase)
   % q at near-saturation is bounded by the saturated hydraulic conductivity.
   f_ice = 0.30;
   f_res_pore = 0.07;
   f_por = 1 - f_ice;
   k_sat = icemodel.column.saturated_hydraulic_conductivity(f_ice, 0);
   f_liq = f_por - 1e-6;  % near saturation
   q = icemodel.column.liquid_flux(f_liq, f_ice, f_res_pore=f_res_pore);
   verifyLessThanOrEqual(testCase, q, k_sat * (1 + 1e-6));
end

function test_default_f_res_pore_uses_parameterLookup(testCase)
   % With f_res_pore omitted, the kernel uses the parameterLookup default
   % (Colbeck-benchmark value, currently 0.07). Confirm by passing the
   % same value explicitly and comparing.
   f_ice = 0.30;
   f_liq = 0.10;
   q_default = icemodel.column.liquid_flux(f_liq, f_ice);
   q_explicit = icemodel.column.liquid_flux(f_liq, f_ice, ...
      f_res_pore=icemodel.parameterLookup('f_res_pore'));
   verifyEqual(testCase, q_default, q_explicit, 'AbsTol', 0);
end

function test_shimizu_requires_grainsz(testCase)
   % Calling shimizu1970 without grainsz must error with the canonical id.
   f_ice = 0.30;
   f_liq = 0.10;
   verifyError(testCase, ...
      @() icemodel.column.liquid_flux(f_liq, f_ice, k_sat_method="shimizu1970"), ...
      'icemodel:column:saturated_hydraulic_conductivity:missingGrainsz');
end
