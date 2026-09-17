function tests = test_candidate_from_icemodel_output
   %TEST_CANDIDATE_FROM_ICEMODEL_OUTPUT Test the icemodel output adapter.
   %
   % icemodel.verification.candidateFromIcemodelOutput maps model ICE1/ICE2
   % fields to the verification schema. These tests give it model-like
   % structs without the synthetic-snow hook, so they need no verification
   % archive and no model run.
   tests = functiontests(localfunctions);
end

function test_esm_site_candidate_maps_model_snow_fields(testCase)
   % A model that stores snow_depth_m, snow_density_kg_m3, and Tsfc gets
   % depth, SWE, and Celsius surface temperature columns.

   [ice1, ice2, opts] = modelLikeOutput();
   manifest = struct("case_type", "esm_site", "comparison_variables", ...
      ["snow_depth_m", "swe_kg_m2", "surface_temp_C"]);

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);
   returned = candidate.data;

   % The timetable keeps the model output times.
   testCase.verifyClass(returned.Time, 'datetime');
   testCase.verifyEqual(returned.Time, ice1.Time);

   % Depth passes through as finite numbers [m].
   testCase.verifyTrue(isnumeric(returned.snow_depth_m) ...
      && all(isfinite(returned.snow_depth_m)));
   testCase.verifyEqual(returned.snow_depth_m, [0.10; 0.20; 0.30]);

   % SWE [kg m-2] is depth [m] times bulk snow density [kg m-3].
   expected = [25; 60; 105];
   testCase.verifyEqual(returned.swe_kg_m2, expected, 'AbsTol', 1e-12);

   % icemodel sets the freezing point at the triple point, 273.16 K.
   expected = [-10; -9; -8];
   testCase.verifyEqual(returned.surface_temp_C, expected, 'AbsTol', 1e-12);
end

function test_esm_site_candidate_accepts_snow_depth_name(testCase)
   % A model that stores snow_depth instead of snow_depth_m gets the same
   % depth and SWE columns.

   [ice1, ice2, opts] = modelLikeOutput();
   ice1.snow_depth = ice1.snow_depth_m;
   ice1 = rmfield(ice1, "snow_depth_m");
   manifest = struct("case_type", "esm_site", "comparison_variables", ...
      ["snow_depth_m", "swe_kg_m2"]);

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);

   returned = candidate.data.snow_depth_m;
   testCase.verifyEqual(returned, [0.10; 0.20; 0.30]);
   returned = candidate.data.swe_kg_m2;
   testCase.verifyEqual(returned, [25; 60; 105], 'AbsTol', 1e-12);
end

function test_esm_site_candidate_samples_column_temperature_by_depth(testCase)
   % soil_temp_<k>_C samples ice2.Tice at the k-th manifest depth on the
   % opts.dz_thermal grid and converts it to Celsius.

   [ice1, ice2, opts] = modelLikeOutput();
   manifest = struct("case_type", "esm_site", ...
      "comparison_variables", "soil_temp_2_C", ...
      "observation_variables", struct("soil_depths_m", [0.00, 0.08]));

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);

   % Depth 0.08 m is the face between nodes 2 and 3 of the 0.04 m mesh. A
   % face belongs to the node below it, so the value comes from node 3, which
   % holds 270.16, 270.66, and 271.16 K.
   returned = candidate.data.soil_temp_2_C;
   expected = [-3.0; -2.5; -2.0];
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
end

function test_esm_site_candidate_omits_depths_outside_the_column(testCase)
   % The fixture column spans 0 to 0.16 m. The model has no value at 5 m,
   % above the surface, or at a NaN depth. The candidate has no column for
   % those soil temperatures.

   [ice1, ice2, opts] = modelLikeOutput();
   manifest = struct("case_type", "esm_site", ...
      "comparison_variables", ["soil_temp_1_C", "soil_temp_2_C", ...
      "soil_temp_3_C"], ...
      "observation_variables", struct("soil_depths_m", [5.0, -0.04, NaN]));

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);

   returned = candidate.data.Properties.VariableNames;
   testCase.verifyEmpty(returned);
end

function [ice1, ice2, opts] = modelLikeOutput()
   %MODELLIKEOUTPUT Build ICE1, ICE2, and OPTS shaped like a model run.
   %
   % ICE1 holds three hourly rows. ICE2.Tice is the depth-by-time column
   % temperature [K] on a 0.04 m grid, the field name that icemodel writes.

   time = datetime(2000, 1, 1, 0, 0, 0) + hours(0:2);
   ice1 = struct( ...
      "Time", time(:), ...
      "snow_depth_m", [0.10; 0.20; 0.30], ...
      "snow_density_kg_m3", [250; 300; 350], ...
      "Tsfc", [263.16; 264.16; 265.16]);
   ice2 = struct("Tice", [ ...
      272.16 272.36 272.56; ...
      271.16 271.46 271.76; ...
      270.16 270.66 271.16; ...
      269.16 269.86 270.56]);
   opts = struct("smbmodel", "icemodel", "sitename", "wfj", ...
      "simyears", 2000, "dz_thermal", 0.04);
end
