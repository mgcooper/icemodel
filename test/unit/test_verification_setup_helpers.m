function tests = test_verification_setup_helpers
   %TEST_VERIFICATION_SETUP_HELPERS Verify shared family-adapter transforms.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Install the package path without creating staged verification artifacts.
   [~, ~, ~, ~, cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
end

function teardownOnce(testCase)
   % Release the shared test environment after direct helper checks.
   testCase.TestData.cleanup = [];
end

function test_state_case_entry_preserves_dry_run_entry(testCase)
   % Dry-run adapters supply complete metadata-only entries that must survive.
   expected = struct('case_id', 'dry', 'forcing_sources', {{'preview'}}, ...
      'eval_sources', {{'preview_obs'}}, 'sentinel', 17);
   state = struct('entry', expected, 'dry_run', true, ...
      'colocation', struct('promice', stagedPromiceLeg()));

   actual = icemodel.verification.setup.stateCaseEntry(state);

   testCase.verifyEqual(actual, expected);
end

function test_state_case_entry_refreshes_colocation_sources(testCase)
   % Persisted adapters must replace stale legs and both derived source lists.
   old_colocation = struct('mar', struct('staged', false));
   entry = struct('case_id', 'persisted', 'colocation', old_colocation, ...
      'forcing_sources', {{'stale'}}, 'eval_sources', {{'stale'}});
   colocation = struct('promice', stagedPromiceLeg());
   state = struct('entry', entry, 'dry_run', false, ...
      'colocation', colocation);

   actual = icemodel.verification.setup.stateCaseEntry(state);

   testCase.verifyEqual(actual.colocation, colocation);
   testCase.verifyEqual(actual.forcing_sources, {'promice'});
   testCase.verifyEqual(actual.eval_sources, {'promice_obs'});
end

function test_sumup_comparison_variables_select_nonempty_groups(testCase)
   % Present, empty, and absent groups exercise every discovery outcome while
   % retaining the canonical density-temperature-SMB order.
   observations = struct('density', table(1), ...
      'subsurface_temperature', table(), 'smb', table(2));

   actual = ...
      icemodel.verification.setup.sumupComparisonVariables(observations);

   testCase.verifyEqual(actual, ["density"; "smb"]);
   testCase.verifyEqual( ...
      icemodel.verification.setup.sumupComparisonVariables(struct()), ...
      strings(0, 1));
end

function leg = stagedPromiceLeg()
   %STAGEDPROMICELEG Return one leg usable for both forcing and evaluation.
   leg = struct('staged', true, 'eval_staged', true, ...
      'met_files', "promice.mat", 'data_files', "promice_data.mat", ...
      'source_id', "promice");
end
