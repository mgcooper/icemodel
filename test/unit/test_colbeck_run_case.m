function tests = test_colbeck_run_case
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   testCase.TestData.def = icemodel.verification.colbeck.caseDefinition();
end

function test_numerical_bundle_shape(testCase)
   bundle = icemodel.verification.colbeck.runCase(struct(), kind="numerical");
   verifyEqual(testCase, bundle.format, 'experiment_bundle');
   verifyEqual(testCase, bundle.metadata.solver_kind, "numerical_icemodel");
   for name = ["exp1", "exp2", "exp3"]
      tt = bundle.experiments.(char(name));
      verifyEqual(testCase, height(tt), 600);
      verifyTrue(testCase, all(ismember( ...
         ["snow_liquid_water_storage_m", "bottom_outflow_mps"], ...
         string(tt.Properties.VariableNames))));
      verifyEqual(testCase, tt.Properties.DimensionNames{1}, 'Time');
   end
end

function test_analytical_bundle_shape(testCase)
   bundle = icemodel.verification.colbeck.runCase(struct(), kind="analytical");
   verifyEqual(testCase, bundle.format, 'experiment_bundle');
   verifyEqual(testCase, bundle.metadata.solver_kind, "analytical_icemodel");
   for name = ["exp1", "exp2", "exp3"]
      tt = bundle.experiments.(char(name));
      verifyTrue(testCase, all(isfinite(tt.snow_liquid_water_storage_m)));
      verifyTrue(testCase, all(isfinite(tt.bottom_outflow_mps)));
   end
end

function test_numerical_finite(testCase)
   bundle = icemodel.verification.colbeck.runCase(struct(), kind="numerical");
   for name = ["exp1", "exp2", "exp3"]
      tt = bundle.experiments.(char(name));
      verifyTrue(testCase, all(isfinite(tt.snow_liquid_water_storage_m)));
      verifyTrue(testCase, all(isfinite(tt.bottom_outflow_mps)));
      verifyGreaterThanOrEqual(testCase, ...
         min(tt.snow_liquid_water_storage_m), 0);
   end
end

function test_numerical_matches_analytical_for_exp1_during_rain(testCase)
   % During the rain window the kinematic-wave shock is well resolved by
   % the first-order upwind scheme; storage and outflow agree at the mm
   % level. The post-rain rarefaction is less well resolved (the upwind
   % scheme has numerical diffusion which smears the rarefaction); that
   % regression is exercised by compareSolutions, not asserted here.
   num = icemodel.verification.colbeck.runCase(struct(), ...
      kind="numerical", experiment_names="exp1");
   ana = icemodel.verification.colbeck.runCase(struct(), ...
      kind="analytical", experiment_names="exp1");
   def = icemodel.verification.colbeck.caseDefinition();
   rain_idx = (1:def.time.rain_window_s / def.time.dt_s).';
   diff = num.experiments.exp1.snow_liquid_water_storage_m(rain_idx) ...
        - ana.experiments.exp1.snow_liquid_water_storage_m(rain_idx);
   verifyLessThan(testCase, sqrt(mean(diff.^2)), 1e-3, ...
      'exp1 storage during rain RMSE > 1 mm');
end

function test_experiment_names_subset(testCase)
   bundle = icemodel.verification.colbeck.runCase(struct(), ...
      kind="numerical", experiment_names=["exp1"]);
   verifyEqual(testCase, fieldnames(bundle.experiments), {'exp1'});
end
