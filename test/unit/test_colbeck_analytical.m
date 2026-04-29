function tests = test_colbeck_analytical
   tests = functiontests(localfunctions);
end

function test_exp1_finite_payload(testCase)
   sol = icemodel.verification.colbeck.analyticalSolution("exp1");
   verifyEqual(testCase, numel(sol.time_seconds), 600);
   verifyTrue(testCase, all(isfinite(sol.snow_liquid_water_storage_m)));
   verifyTrue(testCase, all(isfinite(sol.bottom_outflow_mps)));
   verifyEqual(testCase, sol.parameters.kind, "ripe");
end

function test_exp1_steady_saturation_below_one(testCase)
   sol = icemodel.verification.colbeck.analyticalSolution("exp1");
   verifyLessThan(testCase, sol.parameters.S_inf, 1);
   verifyGreaterThan(testCase, sol.parameters.S_inf, 0);
   verifyGreaterThan(testCase, sol.parameters.c_shock_m_per_s, 0);
end

function test_exp1_front_reaches_bottom_within_window(testCase)
   sol = icemodel.verification.colbeck.analyticalSolution("exp1");
   def = icemodel.verification.colbeck.caseDefinition();
   verifyLessThan(testCase, sol.parameters.t_arrival_s, def.time.rain_window_s);
end

function test_exp1_mass_balance(testCase)
   sol = icemodel.verification.colbeck.analyticalSolution("exp1");
   def = icemodel.verification.colbeck.caseDefinition();
   q_top = def.boundary.q_top_m_per_s;
   t = sol.time_seconds;
   t_off = def.time.rain_window_s;
   inflow = q_top * min(t, t_off);
   outflow_cum = cumtrapz(t, sol.bottom_outflow_mps);
   storage_initial = def.experiments(1).f_liq * def.grid.total_depth_m;
   storage_change = sol.snow_liquid_water_storage_m - storage_initial;
   residual = storage_change - (inflow - outflow_cum);
   verifyLessThan(testCase, max(abs(residual)), 1e-3, ...
      'analytical solution must be mass-conservative to within 1 mm');
end

function test_exp1_outflow_zero_before_arrival(testCase)
   sol = icemodel.verification.colbeck.analyticalSolution("exp1");
   t_arr = sol.parameters.t_arrival_s;
   pre_arrival = sol.time_seconds < t_arr - 1;
   verifyTrue(testCase, all(sol.bottom_outflow_mps(pre_arrival) == 0), ...
      'bottom outflow must be zero before the shock reaches the bottom');
end

function test_exp2_finite_payload(testCase)
   % Cold-snow Clark 2017 wetting-front solution returns real values.
   sol = icemodel.verification.colbeck.analyticalSolution("exp2");
   verifyEqual(testCase, sol.parameters.kind, "cold");
   verifyTrue(testCase, all(isfinite(sol.snow_liquid_water_storage_m)));
   verifyTrue(testCase, all(isfinite(sol.bottom_outflow_mps)));
   verifyGreaterThan(testCase, max(sol.snow_liquid_water_storage_m), 0);
end

function test_exp3_finite_payload(testCase)
   sol = icemodel.verification.colbeck.analyticalSolution("exp3");
   % exp3 is cold partial-duration with merged shock; the analytical
   % path falls back to the fine-grid kinematic-wave PDE reference.
   verifyTrue(testCase, ismember(sol.parameters.kind, ...
      ["cold", "cold_partial_pde"]));
   verifyTrue(testCase, all(isfinite(sol.snow_liquid_water_storage_m)));
   verifyTrue(testCase, all(isfinite(sol.bottom_outflow_mps)));
end

function test_cold_outflow_zero_before_front_arrival(testCase)
   sol = icemodel.verification.colbeck.analyticalSolution("exp2");
   t_arr = sol.parameters.t_arrival_s;
   pre = sol.time_seconds < t_arr - 1;
   verifyTrue(testCase, all(sol.bottom_outflow_mps(pre) == 0), ...
      'cold-snow bottom outflow must be zero before the wetting front arrives');
end

function test_unknown_experiment_errors(testCase)
   verifyError(testCase, ...
      @() icemodel.verification.colbeck.analyticalSolution("expX"), ...
      'icemodel:colbeck:unknownExperiment');
end
