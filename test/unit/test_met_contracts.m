function tests = test_met_contracts
   %TEST_MET_CONTRACTS Verify met loading, processing, and saved-result paths.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build one synthetic workspace and reuse it across the met contract
   % tests so the fixture cost stays out of each individual assertion.

   testCase.TestData.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      [2015; 2016], configure=true, nsteps=24, dt_seconds=3600);
end

function teardown(testCase)
   % Remove the shared synthetic workspace after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_loadmet_concatenates_years_and_computes_exchange(testCase)
   % Multi-year loads should concatenate years and derive exchange terms
   % like De exactly once through the namespaced loading path.

   workspace = testCase.TestData.workspace;
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', [2015 2016]);

   [met, opts_out] = icemodel.loadmet(opts);

   testCase.verifyEqual(unique(year(met.Time))', [2015 2016]);
   testCase.verifyEqual(height(met), 48);
   testCase.verifyEqual(numel(opts_out.metfname), 2);
   testCase.verifyTrue(isvariable('De', met));
   testCase.verifyTrue(all(isfinite(met.De)));
   testCase.verifyGreaterThan(min(met.De), 0);
end

function test_loadmet_swaps_inline_modis_from_metfile(testCase)
   % Inline MODIS-style data embedded in the met file should override the
   % requested output variable when userdata is selected.

   workspace = testCase.TestData.workspace;
   [met_src, ~] = icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'include_modis', true, ...
      'metdir', workspace.metdir);

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   met = icemodel.loadmet(opts);

   testCase.verifyEqual(met.albedo, met_src.MODIS, 'AbsTol', 1e-12);
end

function test_loadmet_swaps_external_userdata_file(testCase)
   % External userdata files should override the selected met variable when
   % the requested timetable contains the expected replacement field.

   workspace = testCase.TestData.workspace;
   opts_base = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016);
   met_base = icemodel.loadmet(opts_base);

   userdata_values = min(max(met_base.albedo - 0.08, 0.05), 0.95);
   icemodel.test.fixtures.writeSyntheticUserdataFile( ...
      workspace.userdatadir, 2016, ...
      'sitename', workspace.sitename, ...
      'userdata', 'modis', ...
      'varname', 'modis', ...
      'Time', met_base.Time, ...
      'values', userdata_values);

   opts_swap = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');
   met_swap = icemodel.loadmet(opts_swap);

   testCase.verifyEqual(met_swap.albedo, userdata_values, 'AbsTol', 1e-12);
end

function test_loadmet_errors_when_userdata_file_lacks_Data(testCase)
   % Corrupt userdata files should fail loudly instead of silently falling
   % back to the met data.

   workspace = testCase.TestData.workspace;
   filepath = fullfile(workspace.userdatadir, 'kanm_modis_2016.mat');
   warn_state = warning('off', 'MATLAB:load:variableNotFound');
   cleanup = onCleanup(@() warning(warn_state));
   bogus = 1;
   save(filepath, 'bogus');

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', 2016, userdata='modis', uservars='albedo');

   try
      icemodel.loadmet(opts);
      testCase.verifyFail('expected loadmet to fail when userdata lacks Data');
   catch ME
      testCase.verifyTrue(contains(ME.message, ...
         'userdata file does not contain timetable "Data"'));
   end
   clear cleanup
end

function test_processmet_supports_native_and_hourly_cadence(testCase)
   % PROCESSMET should preserve native cadence when requested and collapse
   % to hourly outputs when the hourly contract is requested instead.

   [met_native, ~] = icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'nsteps', 8, 'dt_seconds', 900);

   native = icemodel.processmet(met_native, newTimeStep="native");
   hourly = icemodel.processmet(met_native, newTimeStep="hourly");

   testCase.verifyEqual(height(native), 8);
   testCase.verifyEqual(height(hourly), 2);
   testCase.verifyTrue(all(ismember( ...
      {'tsfc', 'swu', 'swn', 'lwu', 'lwn', 'netr'}, ...
      string(hourly.Properties.VariableNames))));
end

function test_loadresults_defaults_to_output_years(testCase)
   % LOADRESULTS should prefer opts.output_years when the saved run spans
   % spinup and retained output years.

   workspace = testCase.TestData.workspace;
   pathoutput = fullfile(workspace.outputdir, 'kanm', 'skinmodel', ...
      'loadresults');
   if exist(pathoutput, 'dir') ~= 7
      mkdir(pathoutput);
   end
   if exist(fullfile(pathoutput, '2016'), 'dir') ~= 7
      mkdir(fullfile(pathoutput, '2016'));
   end
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'skinmodel', [2015 2016], ...
      testname='loadresults', saveflag=true, n_spinup_years=1, ...
      output_profile='minimal', solver=1, dt=3600, pathoutput=pathoutput);
   [~, ~, opts] = icemodel.test.helpers.runSmbModel(opts);

   [ice1, ~, met] = icemodel.loadresults(opts);

   testCase.verifyEqual(unique(year(ice1.Time))', 2016);
   testCase.verifyEqual(unique(year(met.Time))', 2016);
end

function test_postprocess_explicit_met_return_is_hourly(testCase)
   % The explicit met output from POSTPROCESS should stay aligned with the
   % hourly postprocessed model outputs.

   localws = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      localws));

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      localws, 'skinmodel', 2016, dt=900, solver=1, testname='postprocess');
   [ice1_raw, ice2_raw, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);

   [ice1_pp, ~, met_pp] = icemodel.postprocess(ice1_raw, ice2_raw, opts, ...
      met.swd, met.lwd, met.albedo, met.Time);

   testCase.verifyEqual(height(ice1_pp), 24);
   testCase.verifyEqual(height(met_pp), 24);
   clear cleanup
end

function test_retimeHourlyFixedStep_matches_legacy_hourly_mean(testCase)
   % The fixed-step hourly aggregation helper should reproduce the legacy
   % timetable RETIME path exactly on aligned 15-minute model output.

   localws = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      localws));

   % Run one synthetic quarter-hour case and convert the raw surface output
   % into the same timetable form POSTPROCESS retimes.
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      localws, 'skinmodel', 2016, dt=900, solver=1, ...
      testname='postprocess_retime_equivalence');
   [ice1_raw, ~, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);
   ice1_tt = rawIce1ToTimetable(ice1_raw, met.Time);

   % Compare the legacy timetable RETIME result to the fixed-step helper.
   legacy = legacyHourlyMean(ice1_tt);
   fixed = icemodel.internal.retimeHourlyFixedStep(ice1_tt);

   testCase.verifyEqual(fixed.Time, legacy.Time);
   testCase.verifyEqual(double(fixed{:,:}), double(legacy{:,:}), ...
      'AbsTol', 5e-5);
   clear cleanup
end

function ice1_tt = rawIce1ToTimetable(ice1_raw, time)
   %RAWICE1TOTIMETABLE Convert raw model output to POSTPROCESS timetable form.

   % Match POSTPROCESS by converting logical convergence flags to single
   % before building the timetable-backed hourly series.
   if isfield(ice1_raw, 'Tice_converged')
      ice1_raw.Tice_converged = single(ice1_raw.Tice_converged);
   end
   if isfield(ice1_raw, 'Tsfc_converged')
      ice1_raw.Tsfc_converged = single(ice1_raw.Tsfc_converged);
   end

   % Convert the raw struct into the exact timetable shape retimed by
   % POSTPROCESS.
   time.TimeZone = 'UTC';
   ice1_tt = struct2table(ice1_raw);
   ice1_tt = table2timetable(ice1_tt, 'RowTimes', time);
end

function ice1_hourly = legacyHourlyMean(ice1_tt)
   %LEGACYHOURLYMEAN Reproduce the timetable RETIME branch from POSTPROCESS.

   % Match the legacy path exactly, including the leap-day removal.
   ice1_hourly = retime(ice1_tt, 'hourly', 'mean');
   ice1_hourly = ice1_hourly(~(month(ice1_hourly.Time) == 2 ...
      & day(ice1_hourly.Time) == 29), :);
end
