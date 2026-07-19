function tests = test_sumup_observation_matching
   %TEST_SUMUP_OBSERVATION_MATCHING Verify reusable SUMup comparison contracts.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   %SETUPONCE Add project code and allocate an isolated evaluation root.
   test_path = mfilename('fullpath');
   project_root = fileparts(fileparts(fileparts(test_path)));
   addpath(genpath(fullfile(project_root, 'icemodel')))
   testCase.TestData.tmpdir = tempname;
   mkdir(testCase.TestData.tmpdir)
end

function teardownOnce(testCase)
   %TEARDOWNONCE Remove only the isolated fixtures created by this test file.
   if isfolder(testCase.TestData.tmpdir)
      rmdir(testCase.TestData.tmpdir, 's')
   end
end

function test_interval_smb_requires_complete_support_and_retains_ids(testCase)
   % Duplicate interval support must retain each source id, while an internal
   % model gap must reject the whole accumulated observation interval.
   t0 = datetime(2020, 1, 1, 'TimeZone', 'UTC');
   smb = table(["first"; "second"], repmat(t0, 2, 1), ...
      repmat(t0 + hours(3), 2, 1), [0.006; 0.006], ...
      'VariableNames', {'observation_id', 'start_date', 'end_date', 'smb'});
   smb.Properties.VariableUnits = {'', '', '', 'm w.e.'};
   Time = (t0:hours(1):t0 + hours(2))';
   candidate = timetable(Time, repmat(0.002, 3, 1), ...
      'VariableNames', {'smb'});
   candidate.Properties.VariableUnits = {'m w.e. h-1'};

   matches = icemodel.verification.matchObservations( ...
      struct('smb', smb), struct('smb', candidate));

   testCase.verifyEqual(matches.intervals.observation_id, ...
      ["first"; "second"])
   testCase.verifyEqual(matches.intervals.modeled, [0.006; 0.006], ...
      'AbsTol', 1e-12)
   testCase.verifyEqual(matches.intervals.unit, repmat("m w.e.", 2, 1))
   testCase.verifyEqual(string(matches.metrics.status), "ok")

   % A globally hourly candidate with a missing internal hour does not cover
   % the half-open four-hour observation interval.
   gapped_smb = table("gap", t0, t0 + hours(4), 0.008, ...
      'VariableNames', {'observation_id', 'start_date', 'end_date', 'smb'});
   gapped_smb.Properties.VariableUnits = {'', '', '', 'm w.e.'};
   gapped_time = t0 + hours([0; 1; 3; 4]);
   gapped_candidate = timetable(gapped_time, repmat(0.002, 4, 1), ...
      'VariableNames', {'smb'});
   gapped_candidate.Properties.VariableUnits = {'m w.e. h-1'};
   gapped = icemodel.verification.matchObservations( ...
      struct('smb', gapped_smb), struct('smb', gapped_candidate));

   testCase.verifyEmpty(gapped.intervals)
   testCase.verifyEqual(string(gapped.metrics.status), "no_overlap")

   % A synthetic malformed candidate with duplicate model postings must fail
   % closed because counting both would inflate an accumulated interval. This
   % defensive API test does not represent staged data: the canonical artifact
   % audit separately requires every time-series axis to be duplicate-free.
   duplicate_time = t0 + hours([0; 1; 1; 2]);
   duplicate_candidate = timetable(duplicate_time, ...
      repmat(0.002, 4, 1), 'VariableNames', {'smb'});
   duplicate_candidate.Properties.VariableUnits = {'m w.e. h-1'};
   duplicate = icemodel.verification.matchObservations( ...
      struct('smb', smb(1, :)), struct('smb', duplicate_candidate));

   testCase.verifyEmpty(duplicate.intervals)
   testCase.verifyEqual(string(duplicate.metrics.status), "no_overlap")
end

function test_profiles_match_exact_identity_and_date_without_pooling(testCase)
   % Multiple physical profiles on one date must pair by exact identity.
   t = datetime(2015, 6, 15, 12, 0, 0, 'TimeZone', 'UTC');
   observations = profileTable(["a"; "a"; "b"; "b"], ...
      repmat(t, 4, 1), [0; 1; 0; 1], [300; 500; 600; 800], "density");
   candidate = profileTable(["a"; "a"; "b"; "b"], ...
      repmat(t, 4, 1), [0; 1; 0; 1], [310; 510; 620; 820], "density");

   matches = icemodel.verification.matchObservations( ...
      struct('density', observations), struct('density', candidate));

   testCase.verifyEqual(matches.profiles.profile_id, ...
      ["a"; "a"; "b"; "b"])
   testCase.verifyEqual(matches.profiles.modeled, [310; 510; 620; 820])
   testCase.verifyEqual(string(matches.metrics.status), ["ok"; "ok"])

   % One model identity on a date may serve every observed physical profile.
   % A full subdaily history must select the exact/nearest timestamp rather
   % than pool repeated depth rows from separate model states.
   earlier = dateshift(t, 'start', 'day');
   one_model = profileTable(repmat("icemodel", 4, 1), ...
      [earlier; earlier; t; t], [0; 1; 0; 1], ...
      [100; 200; 400; 600], "density");
   shared = icemodel.verification.matchObservations( ...
      struct('density', observations), struct('density', one_model));
   testCase.verifyEqual(height(shared.profiles), 4)
   testCase.verifyEqual(shared.profiles.modeled, [400; 600; 400; 600])

   % Multiple model identities without an exact observation identity remain
   % explicit diagnostics instead of being concatenated or choosing the first.
   ambiguous_model = profileTable(["x"; "x"; "y"; "y"], ...
      repmat(t, 4, 1), [0; 1; 0; 1], [1; 2; 3; 4], "density");
   ambiguous = icemodel.verification.matchObservations( ...
      struct('density', observations), struct('density', ambiguous_model));
   testCase.verifyEmpty(ambiguous.profiles)
   testCase.verifyEqual(string(ambiguous.metrics.status), ...
      repmat("ambiguous_candidate_profile", 2, 1))
end

function test_timetable_row_times_define_temperature_profile_dates(testCase)
   % SUMup temperature timetables must group by their Time row dimension.
   t1 = datetime(2012, 5, 1, 'TimeZone', 'UTC');
   t2 = datetime(2013, 5, 1, 'TimeZone', 'UTC');
   Time = [t1; t1; t2; t2];
   profile_id = ["core"; "core"; "core"; "core"];
   depth = [0; 1; 0; 1];
   subsurface_temperature = [-15; -14; -13; -12];
   observations = timetable(Time, profile_id, depth, ...
      subsurface_temperature);
   observations.Properties.VariableUnits = {'', 'm', 'degC'};
   candidate = profileTable(repmat("icemodel", 4, 1), Time, depth, ...
      [-14; -13; -12; -11], "subsurface_temperature");

   matches = icemodel.verification.matchObservations( ...
      struct('subsurface_temperature', observations), ...
      struct('subsurface_temperature', candidate));

   testCase.verifyEqual(unique(matches.profiles.match_id, 'stable'), ...
      ["core|2012-05-01"; "core|2013-05-01"])
   testCase.verifyEqual(height(matches.profiles), 4)
   testCase.verifyEqual(string(matches.metrics.status), ["ok"; "ok"])
end

function test_profile_interpolation_uses_only_common_depth_support(testCase)
   % Candidate values interpolate to observed depths without extrapolation.
   observations = table([-1; 0; 1; 2; 3], [10; 20; 30; 40; 50], ...
      'VariableNames', {'depth', 'density'});
   candidate = table([0; 2], [100; 200], ...
      'VariableNames', {'depth', 'density'});
   observations.Properties.VariableUnits = {'m', 'kg m-3'};
   candidate.Properties.VariableUnits = {'m', 'kg m-3'};

   matches = icemodel.verification.matchObservations( ...
      struct('density', observations), struct('density', candidate));

   testCase.verifyEqual(matches.profiles.depth, [0; 1; 2])
   testCase.verifyEqual(matches.profiles.observed, [20; 30; 40])
   testCase.verifyEqual(matches.profiles.modeled, [100; 150; 200])
end

function test_legacy_profile_bundle_shapes_remain_comparable(testCase)
   % A no-depth temperature timetable remains a time series, while a depth
   % timetable with distinct row tags remains one undated profile.
   Time = datetime(2010, 1, 1:3, 'TimeZone', 'UTC')';
   observed_temperature = timetable(Time, [-10; -9; -8], ...
      'VariableNames', {'subsurface_temperature'});
   modeled_temperature = timetable(Time, [-9; -8; -7], ...
      'VariableNames', {'subsurface_temperature'});
   observed_temperature.Properties.VariableUnits = {'degC'};
   modeled_temperature.Properties.VariableUnits = {'degC'};
   temperature = icemodel.verification.matchObservations( ...
      struct('subsurface_temperature', observed_temperature), ...
      struct('subsurface_temperature', modeled_temperature));
   testCase.verifyEqual(height(temperature.metrics), 1)
   testCase.verifyEqual(temperature.metrics.rmse, 1, 'AbsTol', 1e-12)
   testCase.verifyTrue(all(isnan(temperature.profiles.depth)))

   observed_density = table([0; 2], [300; 500], ...
      'VariableNames', {'depth', 'density'});
   model_profile = timetable(Time, [0; 1; 2], [300; 400; 500], ...
      'VariableNames', {'midpoint', 'density'});
   density = icemodel.verification.matchObservations( ...
      struct('density', observed_density), ...
      struct('density', model_profile));
   testCase.verifyEqual(string(density.metrics.status), "ok")
   testCase.verifyEqual(density.metrics.rmse, 0, 'AbsTol', 1e-12)
   testCase.verifyEqual(density.profiles.depth, [0; 2])
end

function test_profile_failures_remain_explicit_soft_diagnostics(testCase)
   % Missing variables, absent dates, and disjoint depth support must remain
   % distinguishable soft statuses for model-development triage.
   t = datetime(2015, 6, 15, 'TimeZone', 'UTC');
   observations = profileTable(repmat("core", 2, 1), repmat(t, 2, 1), ...
      [0; 1], [300; 500], "density");

   missing = icemodel.verification.matchObservations( ...
      struct('density', observations), struct());
   testCase.verifyEqual(string(missing.metrics.status), ...
      "missing_candidate_variable")

   wrong_date = profileTable(repmat("icemodel", 2, 1), ...
      repmat(t + days(1), 2, 1), [0; 1], [310; 510], "density");
   dated = icemodel.verification.matchObservations( ...
      struct('density', observations), struct('density', wrong_date));
   testCase.verifyEqual(string(dated.metrics.status), ...
      "missing_candidate_date")

   disjoint = profileTable(repmat("icemodel", 2, 1), repmat(t, 2, 1), ...
      [10; 11], [310; 510], "density");
   no_overlap = icemodel.verification.matchObservations( ...
      struct('density', observations), struct('density', disjoint));
   testCase.verifyEqual(string(no_overlap.metrics.status), "no_overlap")
end

function test_candidate_adapter_preserves_all_dated_profile_columns(testCase)
   % The model adapter must retain every depth-by-time density/temperature
   % column with canonical names, UTC dates, depth, and physical units.
   Time = [datetime(2014, 1, 1); datetime(2015, 1, 1)];
   ice1 = struct('Time', Time, 'smb', [0.1; 0.2]);
   ice2 = struct( ...
      'ro_sno', [300, 310; 400, 410; 500, 510], ...
      'T', [260, 261; 262, 263; 264, 265]);
   opts = struct('smbmodel', 'icemodel', 'sitename', 'fixture', ...
      'simyears', 2014:2015, 'dz_thermal', 0.5);
   manifest = struct('case_type', 'firn_observational', ...
      'comparison_variables', ...
      ["density", "subsurface_temperature", "smb"], ...
      'observation_variables', struct());

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);

   density = candidate.data.density;
   temperature = candidate.data.subsurface_temperature;
   Tf = icemodel.physicalConstant('Tf');
   testCase.verifyEqual(height(density), 6)
   testCase.verifyEqual(density.density, ice2.ro_sno(:))
   testCase.verifyEqual(numel(unique(density.datetime)), 2)
   testCase.verifyEqual(density.Properties.VariableUnits, ...
      {'', '', 'm', 'kg m-3'})
   testCase.verifyEqual(temperature.subsurface_temperature, ...
      ice2.T(:) - Tf)
   testCase.verifyEqual(temperature.Properties.VariableUnits, ...
      {'', '', 'm', 'degC'})
   testCase.verifyEqual(density.datetime.TimeZone, 'UTC')

   % A one-column state remains one dated profile rather than being replicated
   % over the complete model time vector.
   one_column = ice2;
   one_column.ro_sno = one_column.ro_sno(:, 1);
   one_column.T = one_column.T(:, 1);
   one = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, one_column, opts, manifest);
   testCase.verifyEqual(height(one.data.density), 3)
   testCase.verifyEqual(unique(one.data.density.datetime), ...
      datetime(2014, 1, 1, 'TimeZone', 'UTC'))

   % A multi-column history that cannot map one-to-one to Time is invalid.
   bad = ice2;
   bad.ro_sno = repmat(bad.ro_sno(:, 1), 1, 3);
   testCase.verifyError(@() ...
      icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, bad, opts, manifest), ...
      'icemodel:verification:candidateFromIcemodelOutput:profileTimeSize')
end

function test_comparecase_delegates_to_public_matcher(testCase)
   % comparecase must expose the same metrics/aligned rows as a direct call.
   eval_root = fullfile(testCase.TestData.tmpdir, 'eval');
   case_root = fullfile(eval_root, 'sumup', 'matcher');
   mkdir(case_root)
   t = datetime(2015, 6, 15, 'TimeZone', 'UTC');
   density = profileTable(repmat("core", 2, 1), repmat(t, 2, 1), ...
      [0; 1], [300; 500], "density");
   smb = table("annual", datetime(2015, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2015, 1, 2, 'TimeZone', 'UTC'), 0.024, ...
      'VariableNames', {'observation_id', 'start_date', 'end_date', 'smb'});
   smb.Properties.VariableUnits = {'', '', '', 'm w.e.'};
   targets = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('density', density, 'smb', smb), ...
      'metadata', struct('source', 'matcher fixture'));
   save(fullfile(case_root, 'observations.mat'), 'targets')
   writeMatcherManifest(eval_root);

   model_density = profileTable(repmat("icemodel", 2, 1), ...
      repmat(t, 2, 1), [0; 1], [310; 510], "density");
   Time = (datetime(2015, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2015, 1, 1, 23, 0, 0, 'TimeZone', 'UTC'))';
   model_smb = timetable(Time, repmat(0.001, 24, 1), ...
      'VariableNames', {'smb'});
   model_smb.Properties.VariableUnits = {'m w.e. h-1'};
   candidate = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('density', model_density, 'smb', model_smb), ...
      'metadata', struct('source', 'model fixture'));

   direct = icemodel.verification.matchObservations( ...
      targets, candidate, variables=["density", "smb"]);
   result = icemodel.verification.comparecase("matcher", ...
      evaluation_data_root=eval_root, dataset_family="sumup", ...
      candidate=candidate, make_plot=false);

   testCase.verifyEqual(result.metrics, direct.metrics)
   testCase.verifyEqual(result.aligned, direct.aligned)
end

function test_focused_profile_and_interval_plots_use_matched_rows(testCase)
   % Focused plotting must work without a separate validation-only plot path.
   t = datetime(2015, 6, 15, 'TimeZone', 'UTC');
   density = profileTable(repmat("core", 2, 1), repmat(t, 2, 1), ...
      [0; 1], [300; 500], "density");
   model_density = profileTable(repmat("icemodel", 2, 1), ...
      repmat(t, 2, 1), [0; 1], [310; 510], "density");
   profile_match = icemodel.verification.matchObservations( ...
      struct('density', density), struct('density', model_density), ...
      make_plot=true, plot_kind="profile", visible="off");
   testCase.addTeardown(@() close(profile_match.figure))
   testCase.verifyTrue(isgraphics(profile_match.figure, 'figure'))
   testCase.verifyGreaterThanOrEqual( ...
      numel(findobj(profile_match.figure, 'Type', 'Line')), 2)

   smb = table("day", t, t + hours(2), 0.002, ...
      'VariableNames', {'observation_id', 'start_date', 'end_date', 'smb'});
   smb.Properties.VariableUnits = {'', '', '', 'm w.e.'};
   Time = [t; t + hours(1)];
   model_smb = timetable(Time, [0.001; 0.001], ...
      'VariableNames', {'smb'});
   model_smb.Properties.VariableUnits = {'m w.e. h-1'};
   interval_match = icemodel.verification.matchObservations( ...
      struct('smb', smb), struct('smb', model_smb), ...
      make_plot=true, plot_kind="interval", visible="off");
   testCase.addTeardown(@() close(interval_match.figure))
   testCase.verifyTrue(isgraphics(interval_match.figure, 'figure'))
   testCase.verifyGreaterThanOrEqual( ...
      numel(findobj(interval_match.figure, 'Type', 'Line')), 2)
end

function profiles = profileTable(ids, timestamps, depths, values, varname)
   %PROFILETABLE Build one canonical dated profile fixture with physical units.
   profiles = table(string(ids), timestamps, depths, values, ...
      'VariableNames', {'profile_id', 'datetime', 'depth', char(varname)});
   if varname == "subsurface_temperature"
      unit = 'degC';
   else
      unit = 'kg m-3';
   end
   profiles.Properties.VariableUnits = {'', '', 'm', unit};
end

function writeMatcherManifest(eval_root)
   %WRITEMATCHERMANIFEST Write a minimal SUMup manifest for comparecase parity.
   values = { ...
      'matcher'
      'firn_observational'
      'MATCHER'
      'MATCHER'
      'accumulation'
      {'firn'}
      'none'
      struct('lat_wgs84', 67, 'lon_wgs84', -48, ...
      'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', 1000)
      struct('start', '2015-01-01 00:00:00', ...
      'end', '2015-06-15 00:00:00')
      'matcher/observations.mat'
      {}
      {'sumup_obs'}
      {'density', 'smb'}
      struct('present_groups', {{'density', 'smb'}}, ...
      'native_variables', {{'density', 'smb'}})
      struct('sumup', struct('kind', 'firn_profile_obs', 'staged', true))
      'irregular'
      'public matcher delegation fixture'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "sumup", "", "", "tiny", "today", entry, ...
      struct('site', {}, 'reason', {}));
   icemodel.verification.setup.writeManifest( ...
      fullfile(eval_root, 'sumup', 'manifest.json'), manifest);
end
