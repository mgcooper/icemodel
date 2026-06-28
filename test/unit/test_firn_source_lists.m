function tests = test_firn_source_lists
   %TEST_FIRN_SOURCE_LISTS Verify manifest source-list derivation.
   tests = functiontests(localfunctions);
end

function test_models_are_eval_sources_when_staged(testCase)
   % MAR/MERRA/RACMO can all be comparison targets. MAR/MERRA are also forcing
   % sources only when they carry met files; RACMO remains eval-only.
   colocation = struct();
   colocation.promice = leg(true, "met_kanu_promice.mat");
   colocation.sumup = struct('staged', true);
   colocation.gcnet = leg(true, "met_dye2_gcnet.mat");
   colocation.mar = leg(true, "met_kanu_mar.mat");
   colocation.merra = leg(true, "met_kanu_merra.mat");
   colocation.racmo = struct('staged', true, 'data_files', "kanu_racmo.mat");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   % Native/met-capable sources with met files should be forcing sources.
   expected_forcing = ["promice"; "gcnet"; "mar"; "merra"];
   testCase.verifyEqual(returned_forcing, expected_forcing);

   % Staged observation/model legs should be eval sources; RACMO is eval only.
   expected_eval = ["promice_obs"; "sumup_obs"; "gcnet_obs"; ...
      "mar"; "merra"; "racmo"];
   testCase.verifyEqual(returned_eval, expected_eval);
end

function test_skipped_model_legs_are_not_sources(testCase)
   % A recorded skipped leg is useful provenance, but it must not appear in
   % forcing_sources or eval_sources because no data exist to compare against.
   colocation = struct();
   colocation.sumup = struct('staged', true);
   colocation.mar = struct('staged', false, 'reason', "missing source");
   colocation.racmo = struct('staged', false, 'reason', "missing source");

   [returned_forcing, returned_eval] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);

   % Skipped or data-only legs should not create forcing sources.
   expected_forcing = strings(0, 1);
   testCase.verifyEqual(returned_forcing, expected_forcing);

   % The staged SUMup observation leg is the only comparison source here.
   expected_eval = "sumup_obs";
   testCase.verifyEqual(returned_eval, expected_eval);
end

function returned = leg(staged, met_file)
   %LEG Build one staged met-capable colocation leg.
   returned = struct('staged', staged, 'met_files', met_file);
end
