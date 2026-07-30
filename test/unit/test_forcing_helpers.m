function tests = test_forcing_helpers
   %TEST_FORCING_HELPERS Verify the icemodel.forcing.helpers contracts.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Bootstrap the repo test environment so exactremap is available when
   % this file is run directly.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.bootstrap_cleanup = cleanup;
end

function setup(testCase)
   % Per-test temporary output directory for the write helpers.

   testCase.TestData.outdir = string(tempname);
   mkdir(testCase.TestData.outdir)
end

function teardown(testCase)
   % Remove the temporary output directory.

   if isfolder(testCase.TestData.outdir)
      rmdir(testCase.TestData.outdir, 's')
   end
end

function test_forcing_namelist_includes_staged_native_labels(testCase)
   % Runtime validators must accept labels written by verification importers.
   names = icemodel.namelists.forcings();
   expected = ["imau", "retmip", "promice", "gcnet", "esm_snowmip"];

   testCase.verifyTrue(all(ismember(expected, names)));
end

function test_runtime_namelists_include_staged_verification_cases(testCase)
   % Public point-run validators must accept staged verification case ids and
   % userdata labels produced by the dataset importers.
   sites = ["kanu", "dye2_long", "dye2_2016", "summit", "fa", ...
      "s21", "s22", "s23", "humphrey"];
   userdata = ["promice", "gcnet", "imau", "retmip", "esm_snowmip"];

   testCase.verifyTrue(all(ismember(sites, icemodel.namelists.sitename())));
   testCase.verifyTrue(all(ismember(userdata, icemodel.namelists.userdata())));
end

function test_verificationSourceDir_uses_repo_data_root(testCase)
   % Blank verification source roots resolve under the top-level repo data tree.
   expected = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'retmip', 'samimi'));

   actual = icemodel.forcing.helpers.verificationSourceDir( ...
      "", ["retmip", "samimi"]);

   testCase.verifyEqual(actual, expected);
   testCase.verifyEqual(icemodel.forcing.helpers.verificationSourceDir( ...
      "/tmp/source", "imau"), "/tmp/source");
end

function test_resampleMetTimestep_preserves_window(testCase)
   % The helper represents every source interval through its exclusive end.
   Time = (datetime(2012, 1, 1, 3, 0, 0):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0))';
   met = timetable(Time, [260; 261; 262], [100; 110; 120], ...
      'VariableNames', {'tair', 'swd'});

   returned = icemodel.forcing.helpers.resampleMetTimestep(met, "15m");

   testCase.verifyEqual(returned.Time(1), met.Time(1));
   testCase.verifyEqual(returned.Time(end), met.Time(end) + minutes(45));
   testCase.verifyEqual(seconds(median(diff(returned.Time))), 900);
   testCase.verifyEqual(height(returned), 4 * height(met));
   testCase.verifyEqual(returned.swd, repelem(met.swd, 4));
end

function test_resampleMetTimestep_preserves_gaps_wdir_and_metadata(testCase)
   % Source holds preserve direction values and leave NaN/omitted support missing.
   Time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours([0; 1; 2; 6; 7]);
   swd = [100; 110; NaN; 160; 170];
   wdir = [350; 10; NaN; 90; 100];
   met = timetable(Time, swd, wdir);
   met.Properties.VariableUnits = {'W m-2', 'degrees'};
   met.Properties.VariableDescriptions = {'downward shortwave', 'wind direction'};
   met.Properties.UserData = struct('source_token', "keep");
   met = addprop(met, 'Site', 'table');
   met.Properties.CustomProperties.Site = "test";

   returned = icemodel.forcing.helpers.resampleMetTimestep(met, "15m");
   first_hour = returned.Time >= Time(1) & returned.Time < Time(2);
   outage = returned.Time >= Time(3) & returned.Time < Time(4);

   testCase.verifyEqual(returned.wdir(first_hour), ...
      repmat(met.wdir(1), nnz(first_hour), 1));
   testCase.verifyEqual(returned.wdir(returned.Time == Time(1) + minutes(30)), ...
      350, 'AbsTol', 1e-9);
   testCase.verifyTrue(all(isnan(returned.swd(outage))));
   testCase.verifyTrue(all(isnan(returned.wdir(outage))));
   testCase.verifyEqual(returned.Properties.VariableUnits, ...
      met.Properties.VariableUnits);
   testCase.verifyEqual(returned.Properties.VariableDescriptions, ...
      met.Properties.VariableDescriptions);
   testCase.verifyEqual(returned.Properties.CustomProperties.Site, "test");
   testCase.verifyEqual(string(returned.Properties.UserData.source_token), ...
      "keep");
   testCase.verifyEqual( ...
      returned.Properties.UserData.met_resample_source_time_gap_count, 1);
   testCase.verifyEqual(string( ...
      returned.Properties.UserData.met_resample_time_semantics), ...
      "interval_start");
end

function test_resampleMetTimestep_native_invalid_and_idempotent(testCase)
   % Blank cadence is exact; unsupported cadence errors; guarded output is stable.
   met = makeSyntheticMet(datetime(2016, 1, 1, TimeZone="UTC"), 4);
   native = icemodel.forcing.helpers.resampleMetTimestep(met, "");
   first = icemodel.forcing.helpers.resampleMetTimestep(met, "15m");
   second = icemodel.forcing.helpers.resampleMetTimestep(first, "15m");

   testCase.verifyEqual(native, met);
   testCase.verifyEqual(second, first);
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.resampleMetTimestep(met, "30m"), ...
      'icemodel:forcing:resampleMetTimestep:unsupportedDt');
end

function test_resampleMetTimestep_conserves_interval_means_for_model(testCase)
   % [0,100] are forward interval means, so no part of 100 enters hour one.
   t0 = datetime(2016, 1, 1, 0, 0, 0, TimeZone="UTC");
   met = makeSyntheticMet(t0, 3);
   met.swd = [0; 100; 200];

   returned = icemodel.forcing.helpers.resampleMetTimestep(met, "15m");
   opts = struct('z0_bulk', 1e-3, 'z_tair', 2, 'z_wind', 10);
   consumed = nan(4, 1);
   for n = 1:4
      [~, consumed(n)] = icemodel.timestepping.getforcings( ...
         returned, n, true, opts);
   end

   testCase.verifyEqual(consumed, zeros(4, 1));
   testCase.verifyEqual(sum(consumed * 900), 0);
   testCase.verifyEqual(returned.swd(5:8), 100 * ones(4, 1));
   testCase.verifyEqual( ...
      string(returned.Properties.UserData.met_resample_policy), ...
      "interval_start_zero_order_hold");
end

function test_step_observation_heights_selects_dynamic_or_scalar_values(testCase)
   % The model boundary must scalarize dynamic PROMICE heights without
   % changing scalar forcing-family heights or requiring unused RH geometry.
   dynamic = struct('z_tair', [2.5; 2.4], ...
      'z_wind', [2.5; 2.4], 'z_relh', [2.5; 2.4], 'z0_bulk', 1e-3);
   returned = icemodel.surface.step_observation_heights(dynamic, 2);
   testCase.verifyEqual(returned.z_tair, 2.4);
   testCase.verifyEqual(returned.z_wind, 2.4);
   testCase.verifyEqual(returned.z_relh, 2.4);

   scalar = rmfield(dynamic, 'z_relh');
   scalar.z_tair = 2;
   scalar.z_wind = 10;
   returned = icemodel.surface.step_observation_heights(scalar, 2);
   testCase.verifyEqual(returned.z_tair, 2);
   testCase.verifyEqual(returned.z_wind, 10);
   testCase.verifyFalse(isfield(returned, 'z_relh'));
end

function test_resampleMetTimestep_accepts_two_row_hourly_cadence(testCase)
   % A compact contiguous hourly smoke source is a known application cadence.
   met = makeSyntheticMet( ...
      datetime(2016, 1, 1, TimeZone="UTC"), 2);

   returned = icemodel.forcing.helpers.resampleMetTimestep(met, "15m");

   testCase.verifyEqual(height(returned), 8);
   testCase.verifyEqual(returned.tair, repelem(met.tair, 4));
end

function test_resampleMetTimestep_rejects_ambiguous_two_row_cadence(testCase)
   % Ninety- and 120-minute lone steps cannot prove cadence versus an outage.
   met = makeSyntheticMet( ...
      datetime(2016, 1, 1, TimeZone="UTC"), 2);

   for step = [90, 120]
      candidate = met;
      candidate.Time(2) = candidate.Time(1) + minutes(step);
      testCase.verifyError(@() ...
         icemodel.forcing.helpers.resampleMetTimestep(candidate, "15m"), ...
         'icemodel:forcing:resampleMetTimestep:ambiguousSourceCadence');
   end
end

function test_resampleMetTimestep_rejects_fractional_native_gap(testCase)
   % A 75-minute step is 15-minute aligned but not a whole hourly omission.
   met = makeSyntheticMet( ...
      datetime(2016, 1, 1, TimeZone="UTC"), 4);
   met.Time(4) = met.Time(3) + minutes(75);

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.resampleMetTimestep(met, "15m"), ...
      'icemodel:forcing:resampleMetTimestep:irregularSourceTime');
end

function test_resampleMetTimestep_rejects_removed_context_option(testCase)
   % Cadence context was an internal draft API and must not remain callable.
   source = makeSyntheticMet( ...
      datetime(2016, 1, 1, TimeZone="UTC"), 3);

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.resampleMetTimestep( ...
      source, "15m", source_cadence_seconds=7200), ...
      'MATLAB:TooManyInputs');
end

function test_resampleMetTimestep_preserves_single_sample(testCase)
   % One row has no inferable support and must not trigger interpolation/extrapolation.
   source = makeSyntheticMet( ...
      datetime(2016, 1, 1, 0, 0, 0, TimeZone="UTC"), 1);

   returned = icemodel.forcing.helpers.resampleMetTimestep(source, "15m");
   native = icemodel.forcing.helpers.resampleMetTimestep(source, "");

   testCase.verifyEqual(returned.Time, source.Time);
   testCase.verifyEqual(returned{:, :}, source{:, :});
   testCase.verifyEqual(string( ...
      returned.Properties.UserData.met_resample_policy), ...
      "single_sample_unchanged");
   testCase.verifyEqual(native, source);
end

function test_rcmSourceCoverage_uses_repo_forcing_root(testCase)
   % Blank RCM source dirs should probe repo-local data/forcing, independent of
   % the active ICEMODEL_DATA_PATH configuration.
   old_data_path = getenv('ICEMODEL_DATA_PATH');
   cleanup_env = onCleanup(@() setenv('ICEMODEL_DATA_PATH', old_data_path));
   setenv('ICEMODEL_DATA_PATH', char(testCase.TestData.outdir));
   mar_dir = fullfile(icemodel.internal.fullpath('data'), 'forcing', 'mar');
   forcing_dir = fileparts(mar_dir);
   had_mar_dir = isfolder(mar_dir);
   had_forcing_dir = isfolder(forcing_dir);
   icemodel.helpers.ensureDirExists(mar_dir);
   marker = fullfile(mar_dir, 'MARv3.11-test-1999.nc');
   fid = fopen(marker, 'w');
   cleaner = onCleanup(@() cleanupMarkerTree( ...
      marker, mar_dir, had_mar_dir, forcing_dir, had_forcing_dir));
   fclose(fid);

   coverage = icemodel.verification.setup.rcmSourceCoverage( ...
      "mar", struct('mar', ""));

   testCase.verifyTrue(any(coverage.mar.years == 1999));
   clear cleaner cleanup_env
end

function test_optionalDate_treats_nat_as_unset(testCase)
   % Runtime opts use NaT for omitted bounds; helpers should fall back instead of
   % filtering every sample out.
   fallback = datetime(2012, 1, 1, 'TimeZone', 'UTC');

   returned = icemodel.forcing.helpers.optionalDate(NaT, fallback);

   testCase.verifyEqual(returned, fallback);
end

function test_projectLocation_adds_epsg3413_fields(testCase)
   % Source and manifest metadata should derive projected coordinates through
   % one helper so observations.mat and manifest site_location stay aligned.
   location = struct('lat_wgs84', 67, 'lon_wgs84', -47, 'elev_m', 1840);

   returned = icemodel.forcing.helpers.projectLocation(location);

   testCase.verifyTrue(isfield(returned, 'x_epsg3413'));
   testCase.verifyTrue(isfield(returned, 'y_epsg3413'));
   testCase.verifyTrue(all(isfinite([returned.x_epsg3413, ...
      returned.y_epsg3413])));
end

function test_normalizeLocations_preserves_scalar_and_expands_point_list( ...
      testCase)
   % Gridded builders must share one scalar-versus-batch interpretation while
   % treating a polygon as one opaque location.
   point = [67, -47];
   [locations, batch] = ...
      icemodel.forcing.helpers.normalizeLocations(point);
   testCase.verifyFalse(batch);
   testCase.verifySize(locations, [1 1]);
   testCase.verifyEqual(locations{1}, point);

   points = [67, -47; 68, -48; 69, -49];
   [locations, batch] = ...
      icemodel.forcing.helpers.normalizeLocations(points);
   testCase.verifyTrue(batch);
   testCase.verifySize(locations, [1 3]);
   testCase.verifyEqual(vertcat(locations{:}), points);

   polygon = polyshape([0 1 1 0], [0 0 1 1]);
   [locations, batch] = ...
      icemodel.forcing.helpers.normalizeLocations(polygon);
   testCase.verifyFalse(batch);
   testCase.verifySize(locations, [1 1]);
   testCase.verifyEqual(locations{1}.Vertices, polygon.Vertices);
end

function test_attachLocationMetadata_preserves_optional_slope(testCase)
   % One attachment helper must preserve source coordinates and real terrain
   % slope while retaining NaN for families that do not supply slope metadata.
   Time = datetime(2020, 1, 1, 'TimeZone', 'UTC');
   Data = timetable(273.15, 'RowTimes', Time, ...
      'VariableNames', {'tair'});
   location = struct('lat_wgs84', 67, 'lon_wgs84', -47, ...
      'x_epsg3413', 11, 'y_epsg3413', 22, 'elev_m', 1840, ...
      'slope', 0.015);

   with_slope = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, location);
   custom = with_slope.Properties.CustomProperties;
   testCase.verifyEqual(custom.X, 11);
   testCase.verifyEqual(custom.Y, 22);
   testCase.verifyEqual(custom.Lat, 67);
   testCase.verifyEqual(custom.Lon, -47);
   testCase.verifyEqual(custom.Elev, 1840);
   testCase.verifyEqual(custom.Slope, 0.015);
   testCase.verifyEqual(string(custom.ScalarUnits), ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"]);

   without_slope = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, rmfield(location, 'slope'));
   testCase.verifyTrue(isnan( ...
      without_slope.Properties.CustomProperties.Slope));
end

function test_columnizeMetadata_transposes_row_vectors_recursively(testCase)
   % Saved metadata should be column-oriented for interactive inspection, while
   % scalar text and true matrices keep their shape.
   metadata = struct( ...
      'labels', ["a", "b", "c"], ...
      'values', [1 2 3], ...
      'text', 'row text', ...
      'matrix', [1 2; 3 4], ...
      'nested', struct('cells', {{'x', 'y'}}), ...
      'records', struct('name', {"one", "two"}, 'value', {1, 2}), ...
      'empty_records', reshape(struct('name', {}, 'value', {}), 1, 0), ...
      'checks', struct('numnan', array2table([1 2], ...
         'VariableNames', ["tair", "swd"])));

   returned = icemodel.forcing.helpers.columnizeMetadata(metadata);

   testCase.verifySize(returned.labels, [3 1]);
   testCase.verifySize(returned.values, [3 1]);
   testCase.verifySize(returned.text, [1 8]);
   testCase.verifySize(returned.matrix, [2 2]);
   testCase.verifySize(returned.nested.cells, [2 1]);
   testCase.verifySize(returned.records, [2 1]);
   testCase.verifySize(returned.empty_records, [0 1]);
   testCase.verifySize(returned.checks.numnan, [2 2]);
   testCase.verifyEqual(string(returned.checks.numnan.variable), ...
      ["tair"; "swd"]);
   testCase.verifyEqual(returned.checks.numnan.value, [1; 2]);
end

%% metfilename

function test_metTimestepSuffix_closed_registry(testCase)
   % Writer and runtime naming share one closed cadence-to-tag registry.
   testCase.verifyEqual( ...
      icemodel.forcing.helpers.metTimestepSuffix(900), "15m");
   testCase.verifyEqual( ...
      icemodel.forcing.helpers.metTimestepSuffix(1800), "30m");
   testCase.verifyEqual( ...
      icemodel.forcing.helpers.metTimestepSuffix(3600), "1hr");
   testCase.verifyEqual( ...
      icemodel.forcing.helpers.metTimestepSuffix("30m"), "30m");
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.metTimestepSuffix(1200), ...
      'icemodel:forcing:metTimestepSuffix:unsupportedTimestep');
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.metTimestepSuffix("hourly"), ...
      'icemodel:forcing:metTimestepSuffix:unsupportedTimestep');
end

function test_metfilename_window_form(testCase)
   % Window form encodes YYYYMMDD start/end stamps and the dt suffix.

   name = icemodel.forcing.helpers.metfilename("kanm", "mar", ...
      datetime(2015, 10, 1), datetime(2016, 9, 30), 3600);
   testCase.verifyEqual(name, "met_kanm_mar_20151001_20160930_1hr.mat");
end

function test_metfilename_yearly_form(testCase)
   % Per-year legacy form encodes the year and all proven model-met cadences.

   name = icemodel.forcing.helpers.metfilename("kanm", "mar", 2016, [], 900);
   testCase.verifyEqual(name, "met_kanm_mar_2016_15m.mat");

   name = icemodel.forcing.helpers.metfilename("kanm", "mar", 2016, [], "1hr");
   testCase.verifyEqual(name, "met_kanm_mar_2016_1hr.mat");

   name = icemodel.forcing.helpers.metfilename( ...
      "dye2_2016", "retmip", 2016, [], 1800);
   testCase.verifyEqual(name, "met_dye2_2016_retmip_2016_30m.mat");
end

function test_metfilename_rejects_bad_dt(testCase)
   % Unsupported timesteps must error, mirroring createMetFileNames.

   testCase.verifyError(@() icemodel.forcing.helpers.metfilename( ...
      "kanm", "mar", 2016, [], 1200), ...
      'icemodel:forcing:metfilename:unsupportedTimestep');
end

function test_metfilename_roundtrip_with_createMetFileNames(testCase)
   % The writer-side names must match what the reader-side
   % icemodel.createMetFileNames builds from equivalent opts.

   opts = struct('sitename', 'kanm', 'forcings', 'mar', ...
      'simyears', 2016, 'dt', 3600);
   expected = icemodel.createMetFileNames(opts);
   actual = icemodel.forcing.helpers.metfilename("kanm", "mar", 2016, [], 3600);
   testCase.verifyEqual(char(actual), expected{1});

   opts.startdate = datetime(2015, 10, 1);
   opts.enddate = datetime(2016, 9, 30);
   expected = icemodel.createMetFileNames(opts);
   actual = icemodel.forcing.helpers.metfilename("kanm", "mar", ...
      opts.startdate, opts.enddate, 3600);
   testCase.verifyEqual(char(actual), expected{1});

   % The proven RetMIP/Samimi native override round-trips through the same API.
   opts = struct('sitename', 'dye2_2016', 'forcings', 'retmip', ...
      'simyears', 2016, 'dt', 1800, ...
      'startdate', datetime(2016, 5, 2), ...
      'enddate', datetime(2016, 5, 3));
   expected = icemodel.createMetFileNames(opts);
   actual = icemodel.forcing.helpers.metfilename("dye2_2016", "retmip", ...
      opts.startdate, opts.enddate, 1800);
   testCase.verifyEqual(char(actual), expected{1});
end

function test_promice_filled_simyears_resolve_window_artifact(testCase)
   % Simyear-only runtime requests must find the one published filled window.
   root = testCase.TestData.outdir;
   met_dir = fullfile(root, 'met', 'promice_filled');
   mkdir(met_dir);
   staged_name = 'met_kanm_promice_filled_20190101_20221231_15m.mat';
   marker = true;
   save(fullfile(met_dir, staged_name), 'marker');
   opts = struct('sitename', 'kanm', 'forcings', 'promice_filled', ...
      'simyears', [2020, 2021], 'dt', 900, 'pathinput', root);

   % Prefer an enclosing published artifact, then the exact requested window.
   names = icemodel.createMetFileNames(opts);
   testCase.verifyEqual(names, {staged_name});
   delete(fullfile(met_dir, staged_name));
   names = icemodel.createMetFileNames(opts);
   testCase.verifyEqual(names, ...
      {'met_kanm_promice_filled_20200101_20211231_15m.mat'});
   opts.dt = 3600;
   testCase.verifyError(@() icemodel.createMetFileNames(opts), ...
      'icemodel:createMetFileNames:unsupportedPromiceFilledCadence');
end

function test_precipitation_consistency_checks_mass_and_shape(testCase)
   % The shared predicate must reject missing, negative, and unbalanced phases.
   valid = icemodel.forcing.helpers.precipitationConsistency( ...
      [0; 0.1; NaN; 0.1; 0.1], ...
      [0; 0.04; 0; -0.01; 0.06], ...
      [0; 0.06; 0; 0.11; 0.06]);
   testCase.verifyEqual(valid, [true; true; false; false; false]);
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.precipitationConsistency(0, [0; 0], 0), ...
      'icemodel:forcing:precipitationConsistency:sizeMismatch');
end

function test_precipitation_validity_checks_partial_and_complete_splits(testCase)
   % Missing phases are honest, but every finite partial split must still
   % admit a nonnegative complement and every complete split must balance.
   valid = icemodel.forcing.helpers.precipitationValidity( ...
      [0.1; 0.1; NaN; -0.1; 0.1; 0.1; 0.1; NaN], ...
      [0.04; 0.11; -0.01; NaN; 0.04; NaN; 0.04; NaN], ...
      [NaN; NaN; NaN; NaN; 0.06; 0.11; 0.05; NaN]);
   testCase.verifyEqual(valid, ...
      [true; false; false; false; true; false; false; true]);
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.precipitationValidity(0, [0; 0], 0), ...
      'icemodel:forcing:precipitationValidity:sizeMismatch');
end

function test_promice_registered_as_forcing_source(testCase)
   % "promice" is an accepted forcings name (the generic PROMICE/GC-Net AWS
   % station-met source) and resolves through createMetFileNames to a
   % met_<site>_promice file, while the legacy met_<site>_<site> form still
   % works. The forcings validator accepts it and rejects an unknown name.

   testCase.verifyTrue(ismember("promice", icemodel.namelists.forcings()));
   testCase.verifyWarningFree( ...
      @() icemodel.validators.mustBeForcingName("promice"));
   testCase.verifyError( ...
      @() icemodel.validators.mustBeForcingName("not_a_forcing"), ...
      'icemodel:validators:mustBeForcingName');

   % createMetFileNames resolves met_<site>_promice (per-year and window forms).
   opts = struct('sitename', 'kanm', 'forcings', 'promice', ...
      'simyears', 2016, 'dt', 3600);
   names = icemodel.createMetFileNames(opts);
   testCase.verifyEqual(names{1}, 'met_kanm_promice_2016_1hr.mat');

   opts.startdate = datetime(2015, 10, 1);
   opts.enddate = datetime(2016, 9, 30);
   names = icemodel.createMetFileNames(opts);
   testCase.verifyEqual(names{1}, ...
      'met_kanm_promice_20151001_20160930_1hr.mat');

   % Legacy met_<site>_<site> is unaffected.
   legacy = icemodel.createMetFileNames(struct('sitename', 'kanm', ...
      'forcings', 'kanm', 'simyears', 2016, 'dt', 3600));
   testCase.verifyEqual(legacy{1}, 'met_kanm_kanm_2016_1hr.mat');
end

function test_promice_forcings_sets_boom_obs_heights(testCase)
   % PROMICE height defaults are unresolved NaN sentinels at setopts time:
   % loadmet resolves them through the runtime POLICY A3 chain (measured ->
   % interpolated -> nominal 2.6 m), so setopts itself must never bake in
   % a height.

   opts = icemodel.setopts('icemodel', 'kanm', 2015, 'promice');
   testCase.verifyEqual(opts.forcings, 'promice');
   testCase.verifyTrue(isnan(opts.z_tair));
   testCase.verifyTrue(isnan(opts.z_wind));
   testCase.verifyTrue(isnan(opts.z_relh));

   mar = icemodel.setopts('icemodel', 'kanm', 2015, 'mar');
   testCase.verifyEqual(mar.z_wind, 10.0);
end

function test_native_verification_forcings_set_obs_heights(testCase)
   % Native verification labels use their source measurement-height policies.
   imau = icemodel.setopts('icemodel', 's21', 2014, 'imau');
   testCase.verifyEqual(imau.z_tair, 2.0);
   testCase.verifyEqual(imau.z_wind, 10.0);
   testCase.verifyEqual(imau.z_relh, imau.z_tair);

   gcnet_retmip = icemodel.setopts('icemodel', 'summit', 2001, 'retmip');
   testCase.verifyEqual(gcnet_retmip.z_tair, 2.0);
   testCase.verifyEqual(gcnet_retmip.z_wind, 10.0);

   dy2_alias = icemodel.setopts('icemodel', 'DY2', 2001, 'retmip');
   testCase.verifyEqual(dy2_alias.z_tair, 2.0);
   testCase.verifyEqual(dy2_alias.z_wind, 10.0);

   samimi_retmip = icemodel.setopts('icemodel', 'dye2_2016', 2016, 'retmip');
   testCase.verifyEqual(samimi_retmip.z_tair, 2.0);
   testCase.verifyEqual(samimi_retmip.z_wind, 10.0);

   promice_retmip = icemodel.setopts('icemodel', 'kanu', 2016, 'retmip');
   testCase.verifyEqual(promice_retmip.z_tair, 2.6);
   testCase.verifyEqual(promice_retmip.z_wind, 2.6);

   esm = icemodel.setopts('icemodel', 'cdp', 2005, 'esm_snowmip');
   testCase.verifyEqual(esm.z_tair, 2.0);
   testCase.verifyEqual(esm.z_wind, 2.0);
end

%% variableUnits

function test_variableUnits_canonical_and_precip_rate(testCase)
   % The shared map labels precipitation channels with the canonical m s-1
   % rate, diagnostic mass fluxes with mWE/h, fluxes with W m-2, and the
   % tice string with K, in the order requested.

   names = ["tair", "swd", "ppt", "snow", "rain", "snowf", "rainf", ...
      "precip", "runoff", "melt", "smb", "evap", "swe", "albedo", ...
      "wspd", "psfc", "rh", "tice1", "tice8", "snow_depth", "ablation"];
   units = string(icemodel.forcing.helpers.variableUnits(names));

   precip = ["ppt", "snow", "rain", "snowf", "rainf", "precip"];
   testCase.verifyTrue(all(units(ismember(names, precip)) == "m s-1"));

   massflux = ["runoff", "melt", "smb", "evap", "snowf_subl"];
   snowf_subl_unit = string(icemodel.forcing.helpers.variableUnits("snowf_subl"));
   testCase.verifyEqual(snowf_subl_unit, "mWE/h");
   testCase.verifyTrue(all(units(ismember(names, massflux)) == "mWE/h"));

   testCase.verifyEqual(units(names == "tair"), "K");
   testCase.verifyEqual(units(names == "swd"), "W m-2");
   testCase.verifyEqual(units(names == "swe"), "kg m-2");
   testCase.verifyEqual(units(names == "wspd"), "m s-1");
   testCase.verifyEqual(units(names == "tice1"), "K");
   testCase.verifyEqual(units(names == "tice8"), "K");
   testCase.verifyEqual(units(names == "snow_depth"), "m");
end

function test_variableUnits_errors_on_unmapped_channel(testCase)
   % An emitted channel missing from the map must error, not ship unlabeled.

   testCase.verifyError(@() icemodel.forcing.helpers.variableUnits( ...
      ["tair", "not_a_real_channel"]), ...
      'icemodel:forcing:variableUnits:unmappedChannel');
end

function test_variableUnits_wraps_single_source(testCase)
   % variableUnits is a thin wrapper: its unit must equal the unit field of
   % the single canonical source for every core channel, with no duplicated
   % list of its own.

   names = ["tair", "swd", "lwd", "albedo", "wspd", "rh", "psfc", "ppt", ...
      "tsfc", "shf", "lhf", "swe", "tice1", "dtice1"];
   vu = string(icemodel.forcing.helpers.variableUnits(names));
   info = icemodel.netcdf.defaults.variable(names);
   testCase.verifyEqual(vu, string({info.unit}));
end

%% single canonical variable-metadata source

function test_variable_consistent_triplet_for_core_channels(testCase)
   % The single source returns a consistent {standard_name, long_name, unit}
   % for the core forcing/eval channels.

   info = icemodel.netcdf.defaults.variable("tair");
   testCase.verifyEqual(info.standard_name, 'air_temperature');
   testCase.verifyEqual(info.unit, 'K');
   testCase.verifyTrue(info.is_cf);

   info = icemodel.netcdf.defaults.variable("swd");
   testCase.verifyEqual(info.standard_name, ...
      'surface_downwelling_shortwave_flux_in_air');
   testCase.verifyEqual(info.unit, 'W m-2');

   info = icemodel.netcdf.defaults.variable("ppt");
   testCase.verifyEqual(info.unit, 'm s-1');
   testCase.verifyEqual(info.standard_name, 'precipitation_flux');
end

function test_variable_temperature_observation_units_are_celsius(testCase)
   % Generic observation temperature columns are staged in degrees C; forcing
   % channels with Kelvin values use canonical names such as tair or tsfc.

   info = icemodel.netcdf.defaults.variable("temperature");
   testCase.verifyEqual(info.unit, 'degC');
   testCase.verifyEqual(info.long_name, 'temperature');
   testCase.verifyEqual(info.standard_name, '');
end

function test_variable_pattern_channels(testCase)
   % Indexed tice/dtice plus source-family observation aliases resolve by
   % pattern so staged tables remain self-describing without enumerating every
   % thermistor / soil-temperature / snow-depth channel.

   info = icemodel.netcdf.defaults.variable(["tice1", "tice12", ...
      "dtice3", "soil_temp_5_C", "snd_auto_m"]);
   testCase.verifyEqual(string({info.unit}), ["K", "K", "m", "degC", "m"]);
   testCase.verifyEqual(info(1).standard_name, 'land_ice_temperature');
   testCase.verifyEqual(info(3).standard_name, 'depth');
   testCase.verifyEqual(info(5).standard_name, 'surface_snow_thickness');
end

function test_variable_marks_non_cf_channels(testCase)
   % Channels with no official CF name are marked is_cf=false, standard ''.

   info = icemodel.netcdf.defaults.variable(["swe", "modis", "rh_ice"]);
   testCase.verifyEqual(string({info.standard_name}), ["", "", ""]);
   testCase.verifyFalse(any([info.is_cf]));
end

function test_variable_validatecf_passes_for_all_cf_names(testCase)
   % Every standard_name we claim as CF must be in the official table
   % (validated programmatically against the cached/fixture CF table).

   map = icemodel.netcdf.defaults.variables();
   names = string(map.keys());
   names = names(:)';
   testCase.verifyWarningFree(@() ...
      icemodel.netcdf.defaults.variable(names, validatecf=true));
end

function test_variable_errors_on_unknown_channel(testCase)
   testCase.verifyError(@() ...
      icemodel.netcdf.defaults.variable("not_a_real_channel"), ...
      'icemodel:netcdf:variable:unknownChannel');
end

%% CF standard-name table

function test_cfStandardNames_fixture_parses(testCase)
   % The committed CF fixture parses to a string set containing the core
   % standard names the canonical map uses.

   repo_root = icemodel.internal.fullpath();
   fixture = fullfile(repo_root, "test", "fixtures", ...
      "cf-standard-names", "cf-standard-name-table.xml");
   testCase.assumeTrue(isfile(fixture), 'CF fixture missing');

   names = icemodel.netcdf.defaults.cfStandardNames(file=fixture);
   testCase.verifyClass(names, 'string');
   testCase.verifyTrue(all(ismember( ...
      ["air_temperature", "surface_temperature", "precipitation_flux", ...
      "land_ice_temperature"], names)));
end

function test_cfStandardNames_offline_fallback_resolves_repo_fixture(testCase)
   % A failed web fetch must resolve the committed fixture inside this repo.

   stub_dir = string(tempname);
   mkdir(stub_dir)
   cleanup_dir = onCleanup(@() rmdir(stub_dir, 's'));
   writeFailingWebread(stub_dir);
   addpath(stub_dir, '-begin')
   cleanup_path = onCleanup(@() rmpath(stub_dir));
   clear icemodel.netcdf.defaults.cfStandardNames

   testCase.verifyWarning(@() ...
      icemodel.netcdf.defaults.cfStandardNames(refresh=true), ...
      'icemodel:netcdf:cfStandardNames:usingFixture');

   clear icemodel.netcdf.defaults.cfStandardNames
   clear cleanup_path cleanup_dir
end

function writeFailingWebread(folder)
   %WRITEFAILINGWEBREAD Shadow webread so the offline fallback is deterministic.

   filename = fullfile(folder, 'webread.m');
   fid = fopen(filename, 'w');
   assert(fid > 0, 'could not create the webread test stub')
   closer = onCleanup(@() fclose(fid));
   fprintf(fid, 'function varargout = webread(varargin)\n');
   fprintf(fid, 'error(''test:offline'', ''offline test stub'')\n');
   fprintf(fid, 'end\n');
   clear closer
end

%% stampMetadata embedding

function test_stampMetadata_embeds_all_three(testCase)
   % Stamping a timetable sets VariableUnits, VariableDescriptions, and the
   % StandardNames custom property, aligned to the variable order.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   tt = timetable(Time, [1;2;3], [4;5;6], 'VariableNames', {'tair', 'swd'});
   tt = icemodel.forcing.helpers.stampMetadata(tt);

   testCase.verifyEqual(tt.Properties.DimensionNames{1}, 'Time');
   testCase.verifyEqual(string(tt.Properties.VariableUnits), ...
      ["K", "W m-2"]);
   testCase.verifyEqual(string(tt.Properties.VariableDescriptions(1)), ...
      "near-surface air temperature");
   testCase.verifyEqual(tt.Properties.CustomProperties.StandardNames, ...
      ["air_temperature", "surface_downwelling_shortwave_flux_in_air"]);
end

function test_stampMetadata_allows_unknown_table_columns_when_not_strict(testCase)
   % Verification tables may carry source keys next to science channels. The
   % known science columns still get units while unknown labels stay blank.

   T = table([320; 330], ["a"; "b"], ...
      'VariableNames', {'density', 'source_label'});

   returned = icemodel.forcing.helpers.stampMetadata(T, strict=false);

   testCase.verifyEqual(string(returned.Properties.VariableUnits), ...
      ["kg m-3", ""]);
   testCase.verifyEqual( ...
      returned.Properties.CustomProperties.StandardNames(1), ...
      "");
end

function test_data2met_carries_embedded_metadata(testCase)
   % A met timetable built by data2met carries the embedded properties from
   % the single source (units, descriptions, and CF standard names).

   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met = icemodel.forcing.data2met(met, validate=false);

   names = string(met.Properties.VariableNames);
   testCase.verifyEqual(string(met.Properties.VariableUnits), ...
      string(icemodel.forcing.helpers.variableUnits(names)));
   testCase.verifyTrue(isprop(met.Properties.CustomProperties, ...
      "StandardNames"));
   testCase.verifyEqual( ...
      numel(met.Properties.CustomProperties.StandardNames), numel(names));
   testCase.verifyEqual( ...
      met.Properties.CustomProperties.StandardNames(names == "tair"), ...
      "air_temperature");
end

function test_data2met_defaults_to_nan_contract_placeholders(testCase)
   % Missing required channels default to explicit NaN placeholders.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   Data = timetable(Time, [265; 266; 267], [0; 0; 0], ...
      'VariableNames', {'tair', 'ppt'});

   met = icemodel.forcing.data2met(Data);

   required = icemodel.forcing.helpers.metvariables();
   testCase.verifyTrue(all(ismember(required, string(met.Properties.VariableNames))));
   testCase.verifyTrue(all(isnan(met.swd)));
   testCase.verifyEqual(string(met.Properties.VariableUnits( ...
      string(met.Properties.VariableNames) == "ppt")), "m s-1");
end

function test_data2met_strict_mode_rejects_missing_contract(testCase)
   % fillwithmissing=false is the opt-in strict contract for direct callers.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   Data = timetable(Time, [265; 266; 267], [0; 0; 0], ...
      'VariableNames', {'tair', 'ppt'});

   f = @() icemodel.forcing.data2met(Data, fillwithmissing=false);

   testCase.verifyError(f, ...
      'icemodel:forcing:validatemet:missingVariables');
end

function test_data2metCollection_preserves_scalar_native_output(testCase)
   % A scalar Data timetable stays scalar and returns finalized met metadata.

   Data = makeSyntheticMet(datetime(2016, 1, 1), 3);
   Data.Properties.UserData = struct('source_family', "synthetic");

   [met, metadata] = icemodel.forcing.helpers.data2metCollection(Data);

   testCase.verifyClass(met, 'timetable');
   testCase.verifyEqual(met.Time, Data.Time);
   testCase.verifyEqual(met.tair, Data.tair);
   testCase.verifyEqual(metadata.source_family, "synthetic");
   testCase.verifyTrue(metadata.fillwithmissing);
   testCase.verifyEqual(metadata.met_variables, ...
      string(met.Properties.VariableNames(:)));
   testCase.verifyEqual(met.Properties.UserData, metadata);
end

function test_data2metCollection_preserves_cell_shape_and_resamples(testCase)
   % A point collection preserves met and metadata shapes through resampling.

   first = makeSyntheticMet(datetime(2016, 1, 1), 3);
   second = makeSyntheticMet(datetime(2017, 1, 1), 3);
   first.Properties.UserData = struct('site', "first");
   second.Properties.UserData = struct('site', "second");
   Data = reshape({first, second}, 2, 1);

   [met, metadata] = icemodel.forcing.helpers.data2metCollection( ...
      Data, dt_out="15m", fillwithmissing=false);

   testCase.verifyClass(met, 'cell');
   testCase.verifySize(met, [2, 1]);
   testCase.verifyClass(metadata, 'struct');
   testCase.verifySize(metadata, [2, 1]);
   testCase.verifyEqual(height(met{1}), 12);
   testCase.verifyEqual(height(met{2}), 12);
   testCase.verifyEqual(met{1}.Time(1:4), ...
      first.Time(1) + minutes((0:15:45)'));
   testCase.verifyEqual(met{2}.tair(1:4), ...
      repmat(second.tair(1), 4, 1));
   testCase.verifyEqual([metadata.site]', ["first"; "second"]);
   testCase.verifyFalse(any([metadata.fillwithmissing]));
   for k = 1:numel(met)
      testCase.verifyEqual(metadata(k).met_variables, ...
         string(met{k}.Properties.VariableNames(:)));
      testCase.verifyEqual(met{k}.Properties.UserData, metadata(k));
   end
end

function test_data2metCollection_rejects_nonstruct_metadata(testCase)
   % Finalized met metadata must use the canonical scalar-struct schema.

   Data = makeSyntheticMet(datetime(2016, 1, 1), 3);
   Data.Properties.UserData = "invalid";

   f = @() icemodel.forcing.helpers.data2metCollection(Data);

   testCase.verifyError(f, ...
      'icemodel:forcing:data2metCollection:invalidMetadata');
end

function test_completeMetVariables_can_add_split_precip_placeholders(testCase)
   % RetMIP-style native products need ppt/rainf/snowf declared even when the
   % split channels are intentionally missing and left for runtime substitution.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   met = timetable(Time, [265; 266; 267], ...
      'VariableNames', {'tair'});

   returned = icemodel.forcing.helpers.completeMetVariables(met, ...
      include_split_precip=true);

   expected = ["ppt", "rainf", "snowf"];
   testCase.verifyTrue(all(ismember(expected, ...
      string(returned.Properties.VariableNames))));
   testCase.verifyTrue(all(isnan(returned.ppt)));
   testCase.verifyTrue(all(isnan(returned.rainf)));
   testCase.verifyTrue(all(isnan(returned.snowf)));
end

%% windFromComponents

function test_windFromComponents_cardinal_directions(testCase)
   % Meteorological convention: direction the wind blows FROM, clockwise
   % from north, in (0, 360].

   u = [0; -1; 0; 1];   % toward: north, west, south, east
   v = [-1; 0; 1; 0];
   [wspd, wdir] = icemodel.forcing.helpers.windFromComponents(u, v);

   testCase.verifyEqual(wspd, ones(4, 1), 'AbsTol', 1e-12);
   testCase.verifyEqual(wdir, [360; 90; 180; 270], 'AbsTol', 1e-12);
end

%% metchecks

function test_metchecks_counts_and_fills_gaps(testCase)
   % NaN counts reflect the pre-fill state; gaps fill linearly with
   % nearest-value end fill.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 4, 0, 0))';
   tair = [270; 271; NaN; 273; NaN];
   met = timetable(Time, tair);

   [met, checks] = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(checks.numnan.tair, 2);
   testCase.verifyEqual(met.tair, [270; 271; 272; 273; 273], 'AbsTol', 1e-12);
end

function test_metchecks_preserves_all_nan_placeholders(testCase)
   % An all-NaN placeholder should be counted, but never filled or clamped into
   % fabricated physical values.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   swd = [NaN; NaN; NaN];
   rh = [NaN; NaN; NaN];
   met = timetable(Time, swd, rh);

   [returned, checks] = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(checks.numnan.swd, 3);
   testCase.verifyEqual(checks.numnan.rh, 3);
   testCase.verifyTrue(all(isnan(returned.swd)));
   testCase.verifyTrue(all(isnan(returned.rh)));
end

function test_metchecks_clamps_ranges(testCase)
   % Recognized variables clamp to the documented physical ranges.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 3, 0, 0))';
   albedo = [Inf; 1.2; 0.5; 0.01];
   rh = [Inf; 105; 50; 2];
   wspd = [-Inf; 0; 5; 5];
   tsfc = [Inf; 280; 270; 271];   % kelvin (min > 100)
   met = timetable(Time, albedo, rh, wspd, tsfc);

   met = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(met.albedo, [0.98; 0.98; 0.5; 0.05], ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(met.rh, [99.99; 99.99; 50; 5], ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(met.wspd, [0.1; 0.1; 5; 5], 'AbsTol', 1e-12);
   testCase.verifyEqual(met.tsfc, [273.16; 273.16; 270; 271], ...
      'AbsTol', 1e-12);
end

function test_metchecks_fills_wdir_circularly(testCase)
   % A gap between 350 and 10 degrees must fill near north (360/0), not
   % at the linear midpoint 180. This is the intentional fix relative to
   % the legacy runoff metchecks.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   wdir = [350; NaN; 10];
   met = timetable(Time, wdir);

   met = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(met.wdir(2), 360, 'AbsTol', 1e-9);
   testCase.verifyEqual(met.wdir([1 3]), [350; 10], 'AbsTol', 1e-9);
end

function test_metchecks_preserves_all_nan_placeholder(testCase)
   % All-NaN placeholder channels stay missing; metchecks only records counts.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   met = timetable(Time, nan(3, 1), 'VariableNames', {'swd'});

   [met, checks] = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(checks.numnan.swd, 3);
   testCase.verifyTrue(all(isnan(met.swd)));
end

%% dailyToHourly

function test_dailyToHourly_holds_endpoints_and_validates_bounds(testCase)
   % Interpolation remains linear inside native support, holds both endpoints,
   % and rejects a bounded diagnostic whose source is already nonphysical.

   t_daily = [datetime(2016, 1, 1); datetime(2016, 1, 2)];
   t_hourly = (datetime(2015, 12, 31, 23, 0, 0):hours(1): ...
      datetime(2016, 1, 2, 23, 0, 0))';
   daily = [0; 24];

   hourly = icemodel.forcing.helpers.dailyToHourly(daily, t_daily, t_hourly);

   expected = [0; (0:24)'; 24 * ones(23, 1)];
   testCase.verifyEqual(hourly, expected, 'AbsTol', 1e-9);
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.dailyToHourly( ...
      [0; 1.1], t_daily, t_hourly, bounds=[0 1]), ...
      'icemodel:forcing:dailyToHourly:sourceOutOfBounds');

   % Conservative remapping may introduce an ulp-scale boundary overshoot;
   % that numerical noise is snapped without relaxing the physical contract.
   bounded = icemodel.forcing.helpers.dailyToHourly( ...
      [-128 * eps; 1 + 128 * eps], t_daily, t_daily, bounds=[0 1]);
   testCase.verifyEqual(bounded, [0; 1]);
end

%% applyMarDailyQualityControl

function test_applyMarDailyQualityControl_preserves_repairs_and_is_idempotent( ...
      testCase)
   % Matching hourly structure survives; mismatched/missing/partial days use
   % the daily rate; a missing daily reference is retained as unverified.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 3, 0, 0, 0, 'TimeZone', 'UTC'))';
   hourly_shape = (1:24)' / sum(1:24);
   runoff = [hourly_shape; zeros(24, 1); 9];
   smb = zeros(numel(Time), 1);
   smb(30) = NaN;
   smb(end) = 7;
   melt = (1:numel(Time))';
   T = timetable(Time, runoff, smb, melt);
   replacements = struct( ...
      'runoff', [ones(24, 1) / 24; ones(24, 1) / 12; 3 / 24], ...
      'smb', [zeros(48, 1); NaN]);

   [returned, metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      T, replacements, sector=1);

   testCase.verifyEqual(returned.runoff(1:24), hourly_shape);
   testCase.verifyEqual(returned.runoff(25:end), ...
      replacements.runoff(25:end));
   testCase.verifyEqual(returned.smb(1:48), zeros(48, 1));
   testCase.verifyEqual(returned.smb(end), 7);
   testCase.verifyEqual(returned.melt, melt);
   testCase.verifyEqual(returned.Time, Time);
   testCase.verifyEqual(string(metadata.mar_qc_method), ...
      "daily_constrained_hourly");
   testCase.verifyEqual(string(metadata.mar_qc_status), "applied");
   testCase.verifyEqual(metadata.mar_qc_channels, ["runoff"; "smb"]);
   testCase.verifyEqual(metadata.mar_qc_sector, 1);
   testCase.verifyEqual(string(metadata.mar_qc_sector_name), ...
      "permanent_ice");
   testCase.verifyEqual(metadata.mar_qc_complete_utc_day_count, 2);
   testCase.verifyEqual(metadata.mar_qc_partial_utc_day_count, 1);
   testCase.verifyEqual(metadata.mar_qc_runoff_day_status, uint8([1; 2; 2]));
   testCase.verifyEqual(metadata.mar_qc_smb_day_status, uint8([1; 2; 3]));
   testCase.verifyEqual(metadata.mar_qc_runoff_daily_reference_mwe, [1; 2; 3]);
   testCase.verifyEqual(metadata.mar_qc_replaced_runoff_count, 25);
   testCase.verifyEqual(metadata.mar_qc_replaced_smb_count, 1);

   [second, second_metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      returned, replacements, metadata, sector=1);
   testCase.verifyEqual(second, returned);
   testCase.verifyEqual(second_metadata, metadata);

   % Sector identity is part of the currentness contract even when a second
   % sector happens to provide the same daily values and leaves data unchanged.
   [other_sector, other_metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      returned, replacements, metadata, sector=2);
   testCase.verifyEqual(other_sector, returned);
   testCase.verifyEqual(other_metadata.mar_qc_sector, 2);
   testCase.verifyEqual(string(other_metadata.mar_qc_sector_name), "tundra");
end

function test_applyMarDailyQualityControl_tolerance_missing_channel_and_fallback( ...
      testCase)
   % The tolerance boundary is inclusive, a missing channel is added, and a
   % reduced sector-2 source remains hourly with every day unverified.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 1, 23, 0, 0, 'TimeZone', 'UTC'))';
   T = timetable(Time, ones(24, 1), 'VariableNames', {'melt'});
   replacements = struct('runoff', ones(24, 1) / 24);

   [returned, metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      T, replacements, sector=2);

   testCase.verifyEqual(returned.runoff, replacements.runoff);
   testCase.verifyEqual(returned.melt, T.melt);
   testCase.verifyEqual(metadata.mar_qc_channels, "runoff");
   testCase.verifyEqual(metadata.mar_qc_replaced_runoff_count, 24);
   testCase.verifyEqual(metadata.mar_qc_replaced_smb_count, 0);
   testCase.verifyEqual(string(metadata.mar_qc_sector_name), "tundra");
   testCase.verifyEqual(metadata.mar_qc_runoff_day_status, uint8(2));

   boundary = timetable(Time, ones(24, 1) / 24 + 1e-3 / 24, ...
      'VariableNames', {'runoff'});
   [preserved, at_boundary] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      boundary, replacements, abs_tolerance_mwe_day=1e-3, ...
      rel_tolerance=0);
   testCase.verifyEqual(preserved.runoff, boundary.runoff);
   testCase.verifyEqual(at_boundary.mar_qc_runoff_day_status, uint8(1));

   beyond = boundary;
   beyond.runoff(1) = beyond.runoff(1) + 1e-10;
   [repaired, over_boundary] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      beyond, replacements, abs_tolerance_mwe_day=1e-3, ...
      rel_tolerance=0);
   testCase.verifyEqual(repaired.runoff, replacements.runoff);
   testCase.verifyEqual(over_boundary.mar_qc_runoff_day_status, uint8(2));

   native = timetable(Time, linspace(0, 1, 24)', linspace(-1, 1, 24)', ...
      'VariableNames', {'runoff', 'smb'});
   [unchanged, empty_metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      native, struct(), sector=2);
   testCase.verifyEqual(unchanged, native);
   testCase.verifyEqual(string(empty_metadata.mar_qc_status), ...
      "not_applicable");
   testCase.verifyEmpty(empty_metadata.mar_qc_channels);
   testCase.verifyEqual(string(empty_metadata.mar_qc_fallback), ...
      "hourly_RUH_SMBH_retained_native_daily_unavailable");
   testCase.verifyEqual(string(empty_metadata.mar_qc_basis), ...
      "MAR native daily RU/SMB unavailable; retained hourly RUH/SMBH");
   testCase.verifyEqual(empty_metadata.mar_qc_runoff_day_status, uint8(3));
   testCase.verifyEqual(string(empty_metadata.mar_qc_sector_name), "tundra");
end

%% alignMarDailyMetadata

function test_alignMarDailyMetadata_clips_multiyear_hourly_ledgers(testCase)
   % A leap-day complete interior remains exact while both clipped boundary
   % days become unverified and every window-local summary is recomputed.
   source_days = (datetime(2019, 12, 30, 'TimeZone', 'UTC'):caldays(1): ...
      datetime(2020, 3, 2, 'TimeZone', 'UTC'))';
   metadata = marAlignmentFixture(source_days);
   retained = (datetime(2020, 2, 28, 12, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 3, 1, 11, 0, 0, 'TimeZone', 'UTC'))';
   middle = find(source_days == datetime(2020, 2, 29, 'TimeZone', 'UTC'));

   aligned = icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, source_days, retained);

   testCase.verifyEqual(aligned.mar_qc_runoff_day_status, uint8([3; 1; 3]));
   testCase.verifyEqual(aligned.mar_qc_smb_day_status, uint8([3; 2; 3]));
   testCase.verifyTrue(all(isnan( ...
      aligned.mar_qc_runoff_daily_reference_mwe([1 3]))));
   testCase.verifyEqual(aligned.mar_qc_runoff_daily_reference_mwe(2), ...
      metadata.mar_qc_runoff_daily_reference_mwe(middle));
   testCase.verifyEqual(aligned.mar_qc_complete_utc_day_count, 1);
   testCase.verifyEqual(aligned.mar_qc_partial_utc_day_count, 2);
   testCase.verifyEqual(aligned.mar_qc_preserved_runoff_day_count, 1);
   testCase.verifyEqual(aligned.mar_qc_replaced_runoff_day_count, 0);
   testCase.verifyEqual(aligned.mar_qc_unverified_runoff_day_count, 2);
   testCase.verifyEqual(aligned.mar_qc_preserved_smb_day_count, 0);
   testCase.verifyEqual(aligned.mar_qc_replaced_smb_day_count, 1);
   testCase.verifyEqual(aligned.mar_qc_unverified_smb_day_count, 2);

   % The complete-day ME ledger remains a mismatch with its exact source values;
   % boundary references/residuals are deliberately erased as unverified.
   testCase.verifyEqual(aligned.mar_diagnostic_melt_day_status, ...
      uint8([3; 2; 3]));
   testCase.verifyTrue(all(isnan( ...
      aligned.mar_diagnostic_melt_daily_reference_mwe([1 3]))));
   testCase.verifyEqual( ...
      aligned.mar_diagnostic_melt_daily_reference_mwe(2), ...
      metadata.mar_diagnostic_melt_daily_reference_mwe(middle));
   testCase.verifyEqual(aligned.mar_diagnostic_melt_residual_mwe_day(2), ...
      metadata.mar_diagnostic_melt_residual_mwe_day(middle));
   testCase.verifyEqual(string( ...
      aligned.mar_diagnostic_melt_validation_status), "mismatch");
   testCase.verifyEqual(aligned.mar_diagnostic_melt_validated_day_count, 0);
   testCase.verifyEqual(aligned.mar_diagnostic_melt_mismatch_day_count, 1);
   testCase.verifyEqual(aligned.mar_diagnostic_melt_unverified_day_count, 2);
   testCase.verifyEqual( ...
      aligned.mar_diagnostic_melt_max_abs_error_mwe_day, 0.25);

   % Cumulative history and unrelated provenance are outside the window ledger.
   testCase.verifyEqual(aligned.mar_qc_replaced_runoff_count, 17);
   testCase.verifyEqual(aligned.mar_qc_replaced_smb_count, 23);
   testCase.verifyEqual(aligned.source_files, metadata.source_files);
   testCase.verifyEqual(aligned.sentinel, metadata.sentinel);
end

function test_alignMarDailyMetadata_handles_15m_idempotence_and_reduced_melt( ...
      testCase)
   % One complete 15-minute day preserves its source ledgers exactly and the
   % already-aligned second pass is byte-idempotent at the struct-value level.
   day = datetime(2020, 2, 29, 'TimeZone', 'UTC');
   metadata = marAlignmentFixture(day);
   retained = (day:minutes(15):day + hours(23) + minutes(45))';
   first = icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, day, retained);
   second = icemodel.forcing.helpers.alignMarDailyMetadata( ...
      first, day, retained);
   testCase.verifyEqual(first, metadata);
   testCase.verifyEqual(second, first);

   % Empty optional ME/MEH ledgers remain explicitly unavailable while the
   % ordinary RU/SMB ledger still aligns on a complete hourly day.
   reduced = marAlignmentFixture(day);
   reduced.mar_diagnostic_melt_day_status = zeros(0, 1, 'uint8');
   reduced.mar_diagnostic_melt_daily_reference_mwe = zeros(0, 1);
   reduced.mar_diagnostic_melt_residual_mwe_day = zeros(0, 1);
   hourly = (day:hours(1):day + hours(23))';
   reduced = icemodel.forcing.helpers.alignMarDailyMetadata( ...
      reduced, day, hourly);
   testCase.verifyEmpty(reduced.mar_diagnostic_melt_day_status);
   testCase.verifyEqual(string( ...
      reduced.mar_diagnostic_melt_validation_status), "not_available");
   testCase.verifyEqual(reduced.mar_diagnostic_melt_validated_day_count, 0);
   testCase.verifyEqual(reduced.mar_diagnostic_melt_mismatch_day_count, 0);
   testCase.verifyEqual(reduced.mar_diagnostic_melt_unverified_day_count, 0);
   testCase.verifyTrue(isnan( ...
      reduced.mar_diagnostic_melt_max_abs_error_mwe_day));

   % A one-row retained axis has no inferable full-day support and is therefore
   % conservatively classified as partial rather than guessed complete.
   singleton = icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, day, day + hours(12));
   testCase.verifyEqual(singleton.mar_qc_complete_utc_day_count, 0);
   testCase.verifyEqual(singleton.mar_qc_partial_utc_day_count, 1);
   testCase.verifyEqual(singleton.mar_qc_runoff_day_status, uint8(3));
end

function test_alignMarDailyMetadata_fails_closed_on_ambiguous_inputs(testCase)
   % Wrong source lengths, schemas, codes, types, axes, and cadences must report
   % deterministic errors instead of silently positional-clipping provenance.
   days = (datetime(2012, 1, 1, 'TimeZone', 'UTC'):caldays(1): ...
      datetime(2012, 1, 2, 'TimeZone', 'UTC'))';
   retained = (days(1):hours(1):days(2) + hours(23))';
   metadata = marAlignmentFixture(days);

   wrong_length = metadata;
   wrong_length.mar_qc_smb_daily_reference_mwe(end) = [];
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      wrong_length, days, retained), ...
      'icemodel:forcing:alignMarDailyMetadata:ledgerLength');
   partial_schema = rmfield(metadata, ...
      'mar_diagnostic_melt_residual_mwe_day');
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      partial_schema, days, retained), ...
      'icemodel:forcing:alignMarDailyMetadata:ledgerSchema');
   wrong_type = metadata;
   wrong_type.mar_qc_runoff_day_status = "invalid";
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      wrong_type, days, retained), ...
      'icemodel:forcing:alignMarDailyMetadata:ledgerType');
   wrong_code = metadata;
   wrong_code.mar_qc_runoff_day_status(1) = 4;
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      wrong_code, days, retained), ...
      'icemodel:forcing:alignMarDailyMetadata:statusCode');
   integer_reference = metadata;
   integer_reference.mar_qc_runoff_daily_reference_mwe = int16([1; 2]);
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      integer_reference, days, retained), ...
      'icemodel:forcing:alignMarDailyMetadata:ledgerType');

   noon_days = days + hours(12);
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, noon_days, retained), ...
      'icemodel:forcing:alignMarDailyMetadata:sourceDayAxis');
   duplicate = retained;
   duplicate(2) = duplicate(1);
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, days, duplicate), ...
      'icemodel:forcing:alignMarDailyMetadata:retainedTimeAxis');
   outside = retained + calyears(1);
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, days, outside), ...
      'icemodel:forcing:alignMarDailyMetadata:retainedDayOutsideSource');
   irregular = days(1) + minutes([0; 20; 50]);
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, days, irregular), ...
      'icemodel:forcing:alignMarDailyMetadata:retainedCadence');
   missing = retained;
   missing(1) = NaT(1, 1, 'TimeZone', 'UTC');
   testCase.verifyError(@() icemodel.forcing.helpers.alignMarDailyMetadata( ...
      metadata, days, missing), ...
      'icemodel:forcing:alignMarDailyMetadata:missingTime');
end

function metadata = marAlignmentFixture(source_days)
   %MARALIGNMENTFIXTURE Create full RU/SMB and ME/MEH daily provenance vectors.
   n = numel(source_days);
   metadata = struct();
   metadata.mar_qc_runoff_day_status = ones(n, 1, 'uint8');
   metadata.mar_qc_runoff_daily_reference_mwe = (1:n)';
   metadata.mar_qc_smb_day_status = repmat(uint8(2), n, 1);
   metadata.mar_qc_smb_daily_reference_mwe = -(1:n)';
   metadata.mar_qc_complete_utc_day_count = n;
   metadata.mar_qc_partial_utc_day_count = 0;
   metadata.mar_qc_preserved_runoff_day_count = n;
   metadata.mar_qc_replaced_runoff_day_count = 0;
   metadata.mar_qc_unverified_runoff_day_count = 0;
   metadata.mar_qc_preserved_smb_day_count = 0;
   metadata.mar_qc_replaced_smb_day_count = n;
   metadata.mar_qc_unverified_smb_day_count = 0;
   metadata.mar_qc_replaced_runoff_count = 17;
   metadata.mar_qc_replaced_smb_count = 23;

   % Make every ME/MEH day a mismatch with a finite residual so clipping must
   % erase boundary evidence and recompute the retained maximum.
   metadata.mar_diagnostic_melt_day_status = repmat(uint8(2), n, 1);
   metadata.mar_diagnostic_melt_daily_reference_mwe = (101:100 + n)';
   metadata.mar_diagnostic_melt_residual_mwe_day = 0.25 * ones(n, 1);
   metadata.mar_diagnostic_melt_validation_status = "mismatch";
   metadata.mar_diagnostic_melt_validated_day_count = 0;
   metadata.mar_diagnostic_melt_mismatch_day_count = n;
   metadata.mar_diagnostic_melt_unverified_day_count = 0;
   metadata.mar_diagnostic_melt_max_abs_error_mwe_day = 0.25;
   metadata.source_files = ["one.nc"; "two.nc"];
   metadata.sentinel = struct('policy', "preserve", 'value', 42);
end

%% applyMarSnowDepthQualityControl

function test_applyMarSnowDepthQualityControl_isolates_middle_and_flags_edge( ...
      testCase)
   % A bracketed 2018 reset is isolated without masking 2017/2019; a large
   % endpoint jump into 2020 is retained but explicitly unverified.
   Time = (datetime(2016, 1, 1):days(1):datetime(2020, 12, 31))';
   base = containers.Map([2016 2017 2018 2019 2020], ...
      [1.0 1.1 13.0 1.2 4.0]);
   snowd = zeros(numel(Time), 1);
   for yyyy = 2016:2020
      rows = year(Time) == yyyy;
      snowd(rows) = base(yyyy) + 0.02 * sin((1:nnz(rows))' / 20);
   end
   T = timetable(Time, snowd);

   [screened, metadata] = ...
      icemodel.forcing.helpers.applyMarSnowDepthQualityControl(T);

   testCase.verifyTrue(all(isnan(screened.snowd(year(Time) == 2018))));
   testCase.verifyTrue(all(isfinite(screened.snowd(year(Time) == 2017))));
   testCase.verifyTrue(all(isfinite(screened.snowd(year(Time) == 2019))));
   testCase.verifyEqual(metadata.mar_snowd_masked_years, 2018);
   testCase.verifyEqual(metadata.mar_snowd_unverified_years, 2020);
   testCase.verifyEqual(string(metadata.mar_snowd_source), "SHSN2");
   testCase.verifyEqual(string(metadata.mar_snowd_shsn3_policy), ...
      "not_used_total_multilayer_snow_firn_thickness");

   corrupted = screened;
   first = find(year(Time) == 2018, 1);
   corrupted.snowd(first) = 99;
   [second, second_metadata] = ...
      icemodel.forcing.helpers.applyMarSnowDepthQualityControl( ...
      corrupted, metadata);
   testCase.verifyTrue(isnan(second.snowd(first)));
   testCase.verifyEqual(second_metadata.mar_snowd_masked_sample_count, ...
      metadata.mar_snowd_masked_sample_count + 1);
end

function test_applyMarSnowDepthQualityControl_handles_absent_and_short_series( ...
      testCase)
   % The helper records no-channel and one-year cases without inventing masks.
   Time = (datetime(2016, 1, 1):days(1):datetime(2016, 1, 3))';
   [without_snow, absent] = ...
      icemodel.forcing.helpers.applyMarSnowDepthQualityControl( ...
      timetable(Time, ones(3, 1), 'VariableNames', {'tair'}));
   testCase.verifyEqual(without_snow.tair, ones(3, 1));
   testCase.verifyEqual(string(absent.mar_snowd_qc_status), ...
      "not_applicable");

   [short, insufficient] = ...
      icemodel.forcing.helpers.applyMarSnowDepthQualityControl( ...
      timetable(Time, ones(3, 1), 'VariableNames', {'snowd'}));
   testCase.verifyEqual(short.snowd, ones(3, 1));
   testCase.verifyEqual(string(insufficient.mar_snowd_qc_status), ...
      "insufficient_context");
   testCase.verifyEmpty(insufficient.mar_snowd_masked_years);
end

%% applyMerraTimeSupport

function test_applyMerraTimeSupport_holds_tavg3_and_is_idempotent(testCase)
   % Proven three-hour source rows replace legacy hourly ramps without source I/O.
   Time = (datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 23, 0, 0, TimeZone="UTC"))';
   source = (0:7)';
   ramp = interp1((0:3:21)', source, (0:23)', 'linear', 'extrap');
   Data = timetable(Time, ramp, 10 + ramp, 20 + ramp, 30 + ramp, ...
      VariableNames={'runoff', 'albedo', 'snowd', 'swe'});
   metadata = struct( ...
      'keep', "yes", ...
      'merra_tavg3_source_grid_policy', ...
      'native_glc_timestamp_inventory', ...
      'merra_tavg3_expected_source_row_count', numel(source), ...
      'merra_tavg3_source_row_count', numel(source), ...
      'merra_tavg3_source_time_gap_count', 0, ...
      'merra_tavg3_missing_source_times', ...
      NaT(0, 1, 'TimeZone', 'UTC'));

   [corrected, stamped, diagnostics] = ...
      icemodel.forcing.helpers.applyMerraTimeSupport(Data, metadata);
   expected = repelem(source, 3);

   testCase.verifyEqual(corrected.runoff, expected);
   testCase.verifyEqual(corrected.albedo, 10 + expected);
   testCase.verifyEqual(corrected.snowd, 20 + expected);
   testCase.verifyEqual(corrected.swe, 30 + expected);
   testCase.verifyGreaterThan(diagnostics.replaced_count, 0);
   testCase.verifyTrue(diagnostics.metadata_changed);
   testCase.verifyEqual(string(stamped.keep), "yes");
   testCase.verifyEqual(string(stamped.merra_source_time_coordinate), ...
      "native_at_reader");

   [second, second_metadata, second_diagnostics] = ...
      icemodel.forcing.helpers.applyMerraTimeSupport(corrected, stamped);
   testCase.verifyEqual(second, corrected);
   testCase.verifyEqual(second_metadata, stamped);
   testCase.verifyEqual(second_diagnostics.replaced_count, 0);
   testCase.verifyFalse(second_diagnostics.metadata_changed);
end

function test_applyMerraTimeSupport_rejects_native_centers(testCase)
   % A leaked :30 coordinate cannot be source-light relabeled without context.
   Time = datetime(2012, 1, 1, 0, 30, 0, TimeZone="UTC") + hours((0:3)');
   Data = timetable(Time, ones(numel(Time), 1), VariableNames="runoff");

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.applyMerraTimeSupport(Data, struct()), ...
      'icemodel:forcing:applyMerraTimeSupport:notIntervalStart');
end

%% applyRacmoPrecipitationQualityControl

function test_applyRacmoPrecipitationQualityControl_is_physical_and_idempotent(testCase)
   % Zero the observed numerical undershoot class without changing legitimate
   % positive, zero, missing, time, or unrelated values.
   Time = (datetime(2012, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 1, 3, 0, 0, 'TimeZone', 'UTC'))';
   ppt = [0; -1.4465385902440174e-9; 4.450231790542602e-6; NaN];
   smb = [1; 2; 3; 4];
   Data = timetable(Time, ppt, smb);

   [returned, metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl(Data);

   testCase.verifyEqual(returned.ppt, [0; 0; ppt(3); NaN]);
   testCase.verifyEqual(returned.smb, smb);
   testCase.verifyEqual(returned.Time, Time);
   testCase.verifyEqual(string(metadata.racmo_ppt_qc_method), ...
      "negative_to_zero");
   testCase.verifyEqual(string(metadata.racmo_ppt_qc_status), "applied");
   testCase.verifyEqual(metadata.racmo_ppt_qc_replaced_count, 1);
   testCase.verifyEqual(metadata.racmo_ppt_qc_input_minimum, ppt(2));
   testCase.verifyEqual(metadata.racmo_ppt_qc_output_minimum, 0);
   testCase.verifyTrue(contains(string(metadata.racmo_ppt_qc_stage), ...
      "spatial_sampling"));

   [second, second_metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
      returned, metadata);
   testCase.verifyEqual(second, returned);
   testCase.verifyEqual(second_metadata, metadata);

   % Incomplete or incompatible prior provenance is never trusted.
   [~, fresh_metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
      returned, struct('racmo_ppt_qc_method', "other"));
   testCase.verifyEqual(string(fresh_metadata.racmo_ppt_qc_status), ...
      "unchanged");
   testCase.verifyEqual(fresh_metadata.racmo_ppt_qc_replaced_count, 0);
   incompatible = metadata;
   incompatible.racmo_ppt_qc_method = "other";
   [~, fresh_metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
      returned, incompatible);
   testCase.verifyEqual(string(fresh_metadata.racmo_ppt_qc_status), ...
      "unchanged");

   % A matching method does not rescue a different basis, impossible count, or
   % output minimum that disagrees with the current repaired payload.
   incompatible = metadata;
   incompatible.racmo_ppt_qc_basis = "different physical rule";
   [~, fresh_metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
      returned, incompatible);
   testCase.verifyEqual(string(fresh_metadata.racmo_ppt_qc_status), ...
      "unchanged");
   incompatible = metadata;
   incompatible.racmo_ppt_qc_replaced_count = -1;
   [~, fresh_metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
      returned, incompatible);
   testCase.verifyEqual(fresh_metadata.racmo_ppt_qc_replaced_count, 0);
   incompatible = metadata;
   incompatible.racmo_ppt_qc_output_minimum = 1;
   [~, fresh_metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
      returned, incompatible);
   testCase.verifyEqual(fresh_metadata.racmo_ppt_qc_output_minimum, 0);
   incompatible = metadata;
   incompatible.racmo_ppt_qc_status = struct();
   [~, fresh_metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
      returned, incompatible);
   testCase.verifyEqual(string(fresh_metadata.racmo_ppt_qc_status), ...
      "unchanged");
end

function test_applyRacmoPrecipitationQualityControl_marks_reduced_source(testCase)
   % A reduced fixture without precipitation is unchanged and explicitly not
   % applicable rather than receiving a fabricated ppt channel.
   Time = datetime(2012, 1, 1, 'TimeZone', 'UTC');
   Data = timetable(Time, 1, 'VariableNames', {'smb'});

   [returned, metadata] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl(Data);

   testCase.verifyEqual(returned, Data);
   testCase.verifyEqual(string(metadata.racmo_ppt_qc_status), ...
      "not_applicable");
   testCase.verifyEqual(metadata.racmo_ppt_qc_replaced_count, 0);
   testCase.verifyTrue(isnan(metadata.racmo_ppt_qc_input_minimum));
end

%% promiceShortwave

function test_intervalMaximumSolarElevation_matches_support_samples(testCase)
   % Interval geometry uses the same start/quarter/end support as the
   % PROMICE darkness rule, for both hourly and half-hourly postings.
   times = [ ...
      datetime(2020, 3, 1, 0, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2020, 3, 1, 1, 0, 0, 'TimeZone', 'UTC')];
   for interval = [hours(1), minutes(30)]
      sample_times = times + (0:0.25:1) .* interval;
      expected = max(icemodel.forcing.helpers.solarElevation( ...
         sample_times, 66.48, -42.50), [], 2);
      returned = ...
         icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
         times, 66.48, -42.50, interval);
      testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   end
end

function test_promiceShortwave_prefers_corrected_and_clamps_raw_fallback(testCase)
   % Corrected finite values win; raw values fill correction gaps; only a
   % remaining negative public fallback is zeroed; residual nonfinite source
   % values normalize to NaN and genuine missing stays NaN.
   Time = (datetime(2020, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 1, 1, 3, 0, 0, 'TimeZone', 'UTC'))';
   swd = [-5; 0; 20; Inf];
   swd_cor = [NaN; -3; NaN; Inf];
   swu = [-2; 0; 5; Inf];
   swu_cor = [NaN; -1; NaN; Inf];
   aws = timetable(Time, swd, swd_cor, swu, swu_cor);

   [public_swd, public_swu, metadata, masks] = ...
      icemodel.forcing.helpers.promiceShortwave(aws);

   testCase.verifyEqual(public_swd, [0; 0; 20; NaN]);
   testCase.verifyEqual(public_swu, [0; 0; 5; NaN]);
   testCase.verifyTrue(metadata.swd_source_present);
   testCase.verifyTrue(metadata.swd_corrected_source_present);
   testCase.verifyTrue(metadata.swd_observations_present);
   testCase.verifyEqual(metadata.swd_corrected_used_count, 1);
   testCase.verifyEqual(metadata.swd_raw_fallback_count, 2);
   testCase.verifyEqual(metadata.swd_raw_negative_count, 1);
   testCase.verifyEqual(metadata.swd_corrected_negative_count, 1);
   testCase.verifyEqual(metadata.swd_raw_fallback_negative_count, 1);
   testCase.verifyEqual(metadata.swd_negative_clamped_count, 2);
   testCase.verifyFalse(metadata.swd_darkness_fill_enabled);
   testCase.verifyEqual(metadata.swd_darkness_zero_filled_count, 0);
   testCase.verifyEqual(metadata.swd_remaining_missing_count, 1);
   testCase.verifyEqual(metadata.swd_raw_minimum, -5);
   testCase.verifyTrue(contains(string(metadata.swd_policy), "dsr_cor"));
   testCase.verifyEqual(metadata.swu_corrected_negative_count, 1);
   testCase.verifyEqual(metadata.swu_raw_fallback_negative_count, 1);
   testCase.verifyEqual(metadata.swu_negative_clamped_count, 2);
   testCase.verifyEqual(masks.swd_corrected_used, ...
      [false; true; false; false]);
   testCase.verifyEqual(masks.swd_raw_fallback, ...
      [true; false; true; false]);
   testCase.verifyEqual(masks.swd_negative_clamped, ...
      [true; true; false; false]);
end

function test_promiceShortwave_fills_only_deep_civil_night_missing(testCase)
   % A missing complete deep-night hour becomes physical zero. A finite source
   % value is never edited, including positive twilight radiation; missing
   % twilight/daylight hours stay NaN. The August timestamp also guards NOAA's
   % leap-year denominator: fixed-365 would incorrectly classify it as deep.
   Time = [datetime(2008, 1, 1, 3, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2008, 1, 1, 4, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2008, 1, 1, 12, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2008, 3, 1, 12, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2008, 8, 11, 2, 0, 0, 'TimeZone', 'UTC')];
   Time.TimeZone = 'America/New_York';
   swd = [NaN; 8; 4; NaN; NaN];
   swd_cor = nan(size(swd));
   swu = [NaN; 2; 1; NaN; NaN];
   swu_cor = nan(size(swu));
   aws = timetable(Time, swd, swd_cor, swu, swu_cor);

   [public_swd, public_swu, metadata, masks] = ...
      icemodel.forcing.helpers.promiceShortwave(aws, ...
      fill_darkness=true, latitude=67.09695086829153, ...
      longitude=-49.94707712063845);

   testCase.verifyEqual(public_swd, [0; 8; 4; NaN; NaN]);
   testCase.verifyEqual(public_swu, [0; 2; 1; NaN; NaN]);
   testCase.verifyEqual(metadata.swd_corrected_used_count, 0);
   testCase.verifyEqual(metadata.swd_raw_fallback_count, 2);
   testCase.verifyTrue(metadata.swd_darkness_fill_enabled);
   testCase.verifyEqual(metadata.swd_darkness_threshold_degrees, -6);
   testCase.verifyEqual(metadata.swd_darkness_zero_filled_count, 1);
   testCase.verifyEqual(metadata.swd_remaining_missing_count, 2);
   testCase.verifyEqual(metadata.swu_darkness_zero_filled_count, 1);
   testCase.verifyEqual(masks.swd_darkness_fill, ...
      [true; false; false; false; false]);
   testCase.verifyEqual(masks.swu_darkness_fill, ...
      [true; false; false; false; false]);
   testCase.verifyTrue(contains(string( ...
      metadata.swd_darkness_fill_method), "NOAA"));
   testCase.verifyTrue(contains(string( ...
      metadata.swd_darkness_fill_method), "whole-hour"));
   testCase.verifyTrue(contains(string(metadata.swd_policy), ...
      "twilight and daylight missing samples remain missing"));
end

function test_promiceShortwave_darkness_fill_requires_location(testCase)
   % Opt-in solar classification cannot silently guess an AWS location.
   Time = datetime(2020, 1, 1, 'TimeZone', 'UTC');
   aws = timetable(Time, 0, 'VariableNames', {'swd'});

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.promiceShortwave(aws, fill_darkness=true), ...
      'icemodel:forcing:helpers:promiceShortwave:missingLocation');
end

function test_promiceShortwave_uses_file_support_for_surgical_slice(testCase)
   % An outage-only slice has no finite local observations, but explicit
   % whole-source-file support must still permit only its proven deep-night zero.
   Time = datetime(2021, 3, 20, [3; 8], 0, 0, 'TimeZone', 'UTC');
   swd = nan(2, 1);
   swu = nan(2, 1);
   aws = timetable(Time, swd, swu);

   [public_swd, public_swu, metadata] = ...
      icemodel.forcing.helpers.promiceShortwave(aws, ...
      fill_darkness=true, latitude=67.09695086829153, ...
      longitude=-49.94707712063845, ...
      swd_source_file_observations_present=true, ...
      swu_source_file_observations_present=true);

   testCase.verifyEqual(public_swd, [0; NaN]);
   testCase.verifyEqual(public_swu, [0; NaN]);
   testCase.verifyFalse(metadata.swd_observations_present);
   testCase.verifyTrue(metadata.swd_source_file_observations_present);
   testCase.verifyEqual(metadata.swd_corrected_used_count, 0);
   testCase.verifyEqual(metadata.swd_raw_fallback_count, 0);
   testCase.verifyEqual(metadata.swd_darkness_zero_filled_count, 1);
   testCase.verifyEqual(metadata.swd_remaining_missing_count, 1);
   testCase.verifyTrue(contains(string(metadata.swd_policy), ...
      "requested window"));
   testCase.verifyFalse(contains(string(metadata.swd_policy), "placeholder"));
end

function test_promiceShortwave_supports_raw_only_source(testCase)
   % Older source products without corrected variables use raw measurements and
   % remove only finite negative sensor offsets.
   Time = (datetime(2020, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 1, 1, 2, 0, 0, 'TimeZone', 'UTC'))';
   swd = [-1; 0; 2];
   swu = [-0.5; 0; 1];
   aws = timetable(Time, swd, swu);

   [public_swd, public_swu, metadata] = ...
      icemodel.forcing.helpers.promiceShortwave(aws);

   testCase.verifyEqual(public_swd, [0; 0; 2]);
   testCase.verifyEqual(public_swu, [0; 0; 1]);
   testCase.verifyFalse(metadata.swd_corrected_source_present);
   testCase.verifyEqual(metadata.swd_corrected_used_count, 0);
   testCase.verifyTrue(contains(string(metadata.swd_policy), ...
      "negative sensor-offset"));
end

function test_promiceShortwave_marks_missing_channels_as_placeholders(testCase)
   % Sparse stations retain all-missing or absent radiation as explicit
   % placeholders and do not invent finite darkness observations.
   Time = (datetime(2020, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 1, 1, 1, 0, 0, 'TimeZone', 'UTC'))';
   aws = timetable(Time, nan(2, 1), [270; 271], ...
      'VariableNames', {'swd', 'tair'});

   [public_swd, public_swu, metadata] = ...
      icemodel.forcing.helpers.promiceShortwave(aws, ...
      fill_darkness=true, latitude=67.1, longitude=-49.9);

   testCase.verifyTrue(all(isnan(public_swd)));
   testCase.verifyTrue(all(isnan(public_swu)));
   testCase.verifyTrue(metadata.swd_source_present);
   testCase.verifyFalse(metadata.swu_source_present);
   testCase.verifyFalse(metadata.swd_observations_present);
   testCase.verifyEqual(metadata.swd_darkness_zero_filled_count, 0);
   testCase.verifyEqual(metadata.swd_remaining_missing_count, 2);
   testCase.verifyTrue(isnan(metadata.swd_raw_minimum));
   testCase.verifyTrue(all(contains(string(metadata.swd_policy), ...
      ["swd", "placeholder"])));
end

%% solarElevation

function test_solarElevation_interprets_times_as_utc_instants(testCase)
   % The shared NOAA helper preserves zoned instants and distinguishes
   % equatorial midnight from solar noon without an external ephemeris.
   utc_time = datetime(2021, 6, 21, [0; 12], 0, 0, 'TimeZone', 'UTC');
   local_time = utc_time;
   local_time.TimeZone = 'America/New_York';

   utc_elevation = icemodel.forcing.helpers.solarElevation( ...
      utc_time, 0, 0);
   local_elevation = icemodel.forcing.helpers.solarElevation( ...
      local_time, 0, 0);

   testCase.verifyEqual(local_elevation, utc_elevation, 'AbsTol', 1e-12);
   testCase.verifyLessThan(utc_elevation(1), -60);
   testCase.verifyGreaterThan(utc_elevation(2), 60);
end

%% interpRcm

function test_interpRcm_recovers_planar_field(testCase)
   % Linear interpolation of a planar field is exact inside the hull,
   % and the timestep loop (interpolant-value reuse) must scale values.

   [X, Y] = meshgrid(0:4, 0:4);
   x = X(:);
   y = Y(:);
   v1 = 2*x + 3*y;
   v = [v1, 2*v1];   % cells x time
   xq = [1.5; 2.5];
   yq = [2.5; 0.5];

   V = icemodel.forcing.helpers.interpRcm(x, y, v, xq, yq, method="linear");

   expected = [2*xq + 3*yq, 2*(2*xq + 3*yq)]';   % time x cells
   testCase.verifyEqual(V, expected, 'AbsTol', 1e-9);
end

%% gridLocation point methods

function test_gridLocation_nearest_picks_nearest_cell(testCase)
   % "nearest" returns the single nearest cell; collapse yields its value.

   [X, Y] = ndgrid(0:6, 0:6);
   loc = [2.3, 3.7];   % nearest cell is (2, 4)
   [start, count, collapse, ~, loctype] = ...
      icemodel.forcing.helpers.gridLocation(X, Y, loc, "nearest");

   testCase.verifyEqual(loctype, "point");
   testCase.verifyEqual(count, [1 1]);
   testCase.verifyEqual([X(start(1), start(2)), Y(start(1), start(2))], [2 4]);

   block = 2 * X(start(1), start(2)) + 3 * Y(start(1), start(2));   % 1 cell
   testCase.verifyEqual(collapse(block), 16, 'AbsTol', 1e-12);
end

function test_gridLocation_nearest_excludes_masked_cell(testCase)
   % Point selection must not choose a closer non-ice cell when a validity
   % mask is supplied by a gridded source such as RACMO.

   [X, Y] = ndgrid(0:4, 0:4);
   validmask = true(size(X));
   validmask(3, 4) = false;   % [2, 3], the masked nearest cell
   loc = [2.1, 3.4];
   [start, count, collapse] = ...
      icemodel.forcing.helpers.gridLocation(X, Y, loc, "nearest", ...
      validmask=validmask);

   testCase.verifyEqual(count, [1 1]);
   testCase.verifyEqual([X(start(1), start(2)), Y(start(1), start(2))], ...
      [2 4]);
   testCase.verifyEqual(collapse(7), 7);
end

function test_gridLocation_nearest_rejects_empty_mask(testCase)
   % A supplied mask with no eligible point cells must fail explicitly rather
   % than silently falling back to an off-mask cell.

   [X, Y] = ndgrid(0:2, 0:2);
   testCase.verifyError(@() icemodel.forcing.helpers.gridLocation( ...
      X, Y, [1 1], "nearest", validmask=false(size(X))), ...
      'icemodel:forcing:gridLocation:noValidPointCells');
end

function test_gridLocation_nearest_rejects_nonlocal_valid_cell(testCase)
   % A mask must not silently snap an off-domain point to a distant valid cell.

   [X, Y] = ndgrid(0:4, 0:4);
   validmask = false(size(X));
   validmask(end, end) = true;
   testCase.verifyError(@() icemodel.forcing.helpers.gridLocation( ...
      X, Y, [0.1 0.1], "nearest", validmask=validmask, maxdistance=2), ...
      'icemodel:forcing:gridLocation:pointOutsideValidDomain');
end

function test_gridLocation_natural_reproduces_planar_field(testCase)
   % "natural" collapses a neighbourhood by natural-neighbour interpolation,
   % which is exact for a linear field: f = 2x + 3y at (2.3, 3.7) = 15.7.

   [X, Y] = ndgrid(0:6, 0:6);
   loc = [2.3, 3.7];
   [start, count, collapse] = ...
      icemodel.forcing.helpers.gridLocation(X, Y, loc, "natural");

   rows = start(1):start(1) + count(1) - 1;
   cols = start(2):start(2) + count(2) - 1;
   xs = X(rows, cols);
   ys = Y(rows, cols);
   block = 2 * xs(:) + 3 * ys(:);   % planar field over the neighbourhood cells

   testCase.verifyEqual(collapse(block), 2 * 2.3 + 3 * 3.7, 'AbsTol', 1e-9);
end

function test_gridLocation_natural_excludes_masked_values(testCase)
   % Natural-neighbour sampling uses only valid local cells, so a masked
   % sentinel cannot contaminate an otherwise planar field.

   [X, Y] = ndgrid(0:6, 0:6);
   loc = [2.3, 3.7];
   validmask = true(size(X));
   validmask(3, 5) = false;
   [start, count, collapse] = ...
      icemodel.forcing.helpers.gridLocation(X, Y, loc, "natural", ...
      validmask=validmask);
   rows = start(1):start(1) + count(1) - 1;
   cols = start(2):start(2) + count(2) - 1;
   xs = X(rows, cols);
   ys = Y(rows, cols);
   block = 2 * xs(:) + 3 * ys(:);
   localmask = validmask(rows, cols);
   block(~localmask(:)) = 1e9;

   testCase.verifyEqual(collapse(block), 2 * 2.3 + 3 * 3.7, 'AbsTol', 1e-9);
end

%% conservative polygon remap (exactremap-guarded)

function test_remapPolygon_conservative_areaavg(testCase)
   % Conservative area-weighted remap of a linear field f = 10x + y over
   % the square [1,3]x[1,3]: the exact area-weighted mean of a linear
   % field equals its value at the centroid (2,2) = 22 (scaled x2 -> 44).
   % Inputs are ndgrid (the icemodel builder convention); remapPolygon
   % transposes to exactremap's meshgrid contract and reshapes the
   % cells-by-time block to the 3-D form exactremap weights correctly.

   [X, Y] = ndgrid(0:5, 0:5);
   f = 10 * X + Y;
   block = [f(:), 2 * f(:)];   % two timesteps scaled [1, 2]
   P = polyshape([1 3 3 1], [1 1 3 3]);

   s = icemodel.forcing.helpers.remapPolygon(X, Y, block, P);
   testCase.verifyEqual(s, [22; 44], 'AbsTol', 1e-6);
end

function test_gridLocation_conservative_vs_equal(testCase)
   % The polygon collapse: "equal" weights only centre-in-polygon cells
   % (here just one -> 22), "conservative" area-weights the overlap (the
   % analytic linear-field mean over [1,3]^2, also 22 at the centroid).

   [X, Y] = ndgrid(0:6, 0:6);
   P = polyshape([1 3 3 1], [1 1 3 3]);

   [s, c, collapse] = icemodel.forcing.helpers.gridLocation( ...
      X, Y, P, "nearest", remap="conservative");
   rows = s(1):s(1) + c(1) - 1;
   cols = s(2):s(2) + c(2) - 1;
   xs = X(rows, cols);
   ys = Y(rows, cols);
   block = 10 * xs(:) + ys(:);
   v = collapse(block);
   testCase.verifyEqual(v, 22, 'AbsTol', 1e-6);
end

%% validatemet

function test_validatemet_rejects_missing_variable(testCase)
   % Required met variables must be present.

   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met = removevars(met, 'psfc');

   testCase.verifyError(@() icemodel.forcing.helpers.validatemet(met), ...
      'icemodel:forcing:validatemet:missingVariables');
end

function test_validatemet_rejects_irregular_time(testCase)
   % The time axis must have a uniform timestep.

   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met.Time(3) = met.Time(3) + minutes(1);

   testCase.verifyError(@() icemodel.forcing.helpers.validatemet(met), ...
      'icemodel:forcing:validatemet:irregularTimeAxis');
end

function test_validatemet_accepts_all_nan_required_placeholder(testCase)
   % An all-NaN required channel is an explicit placeholder, not an absent var.

   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met.swd(:) = NaN;

   testCase.verifyWarningFree(@() icemodel.forcing.helpers.validatemet(met));
end

%% staged artifact identity

function test_artifactMetadata_preserves_userdata_and_fills_point(testCase)
   % The shared writer adapter keeps UserData, fills legacy custom-property
   % coordinates without overriding explicit identity, and handles bare structs.
   Data = makeSyntheticData(datetime(2016, 1, 1), 4);
   Data.Properties.UserData = struct('sample_method', 'nearest');
   metadata = icemodel.forcing.helpers.artifactMetadata(Data);
   testCase.verifyEqual(string(metadata.sample_method), "nearest");
   testCase.verifyEqual(metadata.lat_wgs84, 67.0);
   testCase.verifyEqual(metadata.lon_wgs84, -49.5);
   testCase.verifyEqual(metadata.artifact_cadence_seconds, 3600);

   Data.Properties.UserData.lat_wgs84 = 66.0;
   metadata = icemodel.forcing.helpers.artifactMetadata(Data);
   testCase.verifyEqual(metadata.lat_wgs84, 66.0);
   direct = struct('source_id', 'mar3.11');
   testCase.verifyEqual( ...
      icemodel.forcing.helpers.artifactMetadata(direct), direct);
   testCase.verifyEmpty( ...
      fieldnames(icemodel.forcing.helpers.artifactMetadata(1)));
end

function test_artifactCadenceMatches_current_legacy_and_invalid(testCase)
   % Current top-level cadence stays source-light; legacy files require a uniform
   % timetable, and missing, mislabeled, or irregular payloads are rejected.
   Data = makeSyntheticData(datetime(2016, 1, 1), 8);
   filename = fullfile(testCase.TestData.outdir, 'cadence-current.mat');
   artifact_metadata = ...
      icemodel.forcing.helpers.artifactMetadata(Data);
   save(filename, 'Data', 'artifact_metadata');
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactCadenceMatches( ...
      filename, "Data", 3600));
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactCadenceMatches( ...
      filename, "Data", 1800));

   legacy_file = fullfile(testCase.TestData.outdir, 'cadence-legacy.mat');
   save(legacy_file, 'Data');
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactCadenceMatches( ...
      legacy_file, "Data", 3600));
   Data.Time(4) = Data.Time(4) + minutes(15);
   irregular_file = fullfile(testCase.TestData.outdir, 'cadence-irregular.mat');
   save(irregular_file, 'Data');
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactCadenceMatches( ...
      irregular_file, "Data", 3600));
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactCadenceMatches( ...
      fullfile(testCase.TestData.outdir, 'missing.mat'), "Data", 3600));

   bogus = 1;
   bogus_file = fullfile(testCase.TestData.outdir, 'cadence-bogus.mat');
   save(bogus_file, 'bogus');
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactCadenceMatches( ...
      bogus_file, "Data", 3600));
end

function test_artifactIdentityMatches_compares_concrete_metadata(testCase)
   % Matching current metadata reuses safely; concrete method/product/schema and
   % point changes are rejected while unknown legacy metadata remains compatible.
   filename = fullfile(testCase.TestData.outdir, 'identity-current.mat');
   Data = makeSyntheticData(datetime(2016, 1, 1), 4);
   Data.Properties.UserData = struct('sample_method', 'nearest', ...
      'source_id', 'mar3.11', 'product_id', 'MARv3.11', ...
      'schema_version', '1');
   artifact_metadata = Data.Properties.UserData;
   save(filename, 'Data', 'artifact_metadata');

   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      filename, Data, "Data"));
   expected = Data.Properties.UserData;
   expected.sample_method = 'natural';
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      filename, expected, "Data"));
   expected = struct('product_id', 'MARv3.12');
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      filename, expected, "Data"));
   expected = struct('schema_version', '2');
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      filename, expected, "Data"));

   % Production builder aliases compare against repaired/manifest aliases.
   production_file = fullfile(testCase.TestData.outdir, ...
      'identity-production-alias.mat');
   artifact_metadata = struct('method', 'nearest', ...
      'source_id', 'mar3.11');
   save(production_file, 'artifact_metadata');
   expected = struct('sample_method', 'natural', ...
      'source_id', 'mar3.11');
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      production_file, expected));
   expected.sample_method = 'nearest';
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      production_file, expected));
   expected.source_id = 'mar3.12';
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      production_file, expected));

   artifact_metadata.sample_method = 'natural';
   save(production_file, 'artifact_metadata');
   expected.source_id = 'mar3.11';
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      production_file, expected));

   legacy_file = fullfile(testCase.TestData.outdir, 'identity-legacy.mat');
   legacy = timetable(Data.Time, Data.tair, 'VariableNames', {'tair'});
   save(legacy_file, 'legacy');
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      legacy_file, expected));
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      fullfile(testCase.TestData.outdir, 'missing.mat'), expected));
end

function test_artifactScalarIdentityMatches_normalizes_only_documented_alias(testCase)
   % The shared pure comparator keeps independent fields independent while
   % treating production method and repaired sample_method as one identity.
   existing = struct('method', 'nearest', 'source_id', 'mar3.11', ...
      'schema_version', 1);
   incoming = struct('sample_method', 'nearest', 'source_id', 'mar3.11', ...
      'schema_version', 1);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      existing, incoming));

   incoming.sample_method = 'natural';
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      existing, incoming));
   incoming.sample_method = 'nearest';
   incoming.source_id = 'mar3.12';
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      existing, incoming));

   % Contradictory values inside one alias group are invalid; absent legacy
   % identity and non-scalar values remain unknown-compatible.
   existing.sample_method = 'natural';
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      existing, struct('method', 'nearest')));
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      struct(), struct('schema_version', [1, 2])));

   % Native IMAU/RetMIP producer fields use these exact spellings. Equal facts
   % match, while each independently conflicting concrete fact rejects reuse.
   [fields, matching_values, conflicting_values] = ...
      nativeProducerIdentityCases();
   for k = 1:numel(fields)
      name = char(fields(k));
      existing = struct();
      incoming = struct();
      existing.(name) = matching_values(k);
      incoming.(name) = matching_values(k);
      testCase.verifyTrue( ...
         icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
         existing, incoming));
      incoming.(name) = conflicting_values(k);
      testCase.verifyFalse( ...
         icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
         existing, incoming));
   end
end

function test_artifactIdentityMatches_uses_table_and_rejects_bad_point(testCase)
   % Table CustomProperties provide legacy point identity, and a partial saved
   % coordinate pair is malformed when the requested point is concrete.
   filename = fullfile(testCase.TestData.outdir, 'identity-table.mat');
   Data = makeSyntheticData(datetime(2016, 1, 1), 4);
   save(filename, 'Data');
   expected = struct('lat_wgs84', 67.0, 'lon_wgs84', -49.5);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      filename, expected, "Data"));
   expected.lat_wgs84 = 68.0;
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      filename, expected, "Data"));

   partial_file = fullfile(testCase.TestData.outdir, 'identity-partial.mat');
   artifact_metadata = struct('lat_wgs84', 67.0);
   save(partial_file, 'artifact_metadata');
   expected = struct('lat_wgs84', 67.0, 'lon_wgs84', -49.5);
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      partial_file, expected));
end

%% writemet

function test_writemet_window_writes_and_loads(testCase)
   % Window naming produces one 15-minute model-met file by default.

   met = makeSyntheticMet(datetime(2016, 1, 1), 48);

   filenames = icemodel.forcing.helpers.writemet(met, "tst", "src", ...
      outdir=testCase.TestData.outdir);

   testCase.verifyEqual(numel(filenames), 1);
   [~, name, ext] = fileparts(filenames);
   testCase.verifyEqual(name + ext, "met_tst_src_20160101_20160102_15m.mat");

   loaded = load(filenames, 'met', 'artifact_metadata');
   testCase.verifyEqual(seconds(median(diff(loaded.met.Time))), 900);
   testCase.verifyEqual(loaded.met.Time(1), met.Time(1));
   testCase.verifyEqual(loaded.met.Time(end), met.Time(end) + minutes(45));
   testCase.verifyEqual(string(loaded.met.Properties.VariableUnits( ...
      string(loaded.met.Properties.VariableNames) == "ppt")), "m s-1");
   testCase.verifyEqual(loaded.met.Properties.UserData, ...
      loaded.artifact_metadata);
   testCase.verifyEqual(loaded.artifact_metadata.artifact_cadence_seconds, 900);
end

function test_writemet_resample_preserves_nan_outages_and_endpoints(testCase)
   % The 15-minute default holds finite intervals and preserves NaN support.
   t0 = datetime(2016, 1, 1, 0, 0, 0, TimeZone="UTC");
   met = makeSyntheticMet(t0, 8);
   met.swd(3:4) = NaN;
   met.tair([1, end]) = NaN;
   met.ppt(6) = NaN;

   filename = icemodel.forcing.helpers.writemet(met, "gap", "src", ...
      outdir=testCase.TestData.outdir);
   loaded = load(filename, 'met', 'artifact_metadata');

   testCase.verifyEqual(loaded.met.Time.TimeZone, 'UTC');
   testCase.verifyTrue(all(isnan(loaded.met.swd( ...
      loaded.met.Time >= met.Time(3) & loaded.met.Time < met.Time(5)))));
   testCase.verifyTrue(isfinite(loaded.met.swd(loaded.met.Time == met.Time(2))));
   testCase.verifyTrue(isfinite(loaded.met.swd(loaded.met.Time == met.Time(5))));
   testCase.verifyTrue(isnan(loaded.met.tair(1)));
   testCase.verifyTrue(isnan(loaded.met.tair(end)));
   testCase.verifyTrue(all(isfinite(loaded.met.lwd)));
   testCase.verifyEqual(string(loaded.artifact_metadata.met_resample_policy), ...
      "interval_start_zero_order_hold");
   testCase.verifyGreaterThan( ...
      loaded.artifact_metadata.met_resample_expected_missing_counts.swd, ...
      nnz(isnan(met.swd)));
end

function test_writemet_resample_preserves_omitted_time_gap(testCase)
   % Missing source rows become a line-breaking outage, not a long interpolation.
   met = makeSyntheticMet( ...
      datetime(2016, 1, 1, 0, 0, 0, TimeZone="UTC"), 12);
   met(5:8, :) = [];

   filename = icemodel.forcing.helpers.writemet(met, "omitted", "src", ...
      outdir=testCase.TestData.outdir);
   loaded = load(filename, 'met', 'artifact_metadata');
   interior = loaded.met.Time >= met.Time(4) + hours(1) ...
      & loaded.met.Time < met.Time(5);

   testCase.verifyTrue(all(all(isnan(loaded.met{interior, :}))));
   testCase.verifyTrue(all(isfinite(loaded.met{loaded.met.Time == met.Time(4), :})));
   testCase.verifyTrue(all(isfinite(loaded.met{loaded.met.Time == met.Time(5), :})));
   testCase.verifyEqual( ...
      loaded.artifact_metadata.met_resample_source_time_gap_count, 1);
end

function test_writemet_blank_dt_preserves_native_cadence(testCase)
   % Blank window output preserves native rows and truthful guarded provenance.
   met = makeSyntheticMet(datetime(2016, 1, 1), 48);

   filenames = icemodel.forcing.helpers.writemet(met, "tst", "src", ...
      outdir=testCase.TestData.outdir, dt_out="");

   testCase.verifyTrue(endsWith(filenames, "_1hr.mat"));
   loaded = load(filenames, 'met');
   expected = icemodel.forcing.helpers.stampMetadata(met);
   expected.Properties.UserData = ...
      icemodel.forcing.helpers.artifactMetadata(expected);
   testCase.verifyEqual(loaded.met, expected);

   guarded = icemodel.forcing.helpers.resampleMetTimestep(met, "15m");
   guarded.Properties.UserData = rmfield(guarded.Properties.UserData, ...
      {'met_resample_support_start_inclusive', ...
      'met_resample_yearly_summaries'});
   guarded_file = icemodel.forcing.helpers.writemet( ...
      guarded, "guarded_native", "src", ...
      outdir=testCase.TestData.outdir, dt_out="");
   guarded_loaded = load(guarded_file, 'met', 'artifact_metadata');
   expected = icemodel.forcing.helpers.stampMetadata(guarded);
   expected.Properties.UserData = ...
      icemodel.forcing.helpers.artifactMetadata(expected);
   testCase.verifyEqual(guarded_loaded.met, ...
      expected);
   testCase.verifyEqual(guarded_loaded.artifact_metadata, ...
      icemodel.forcing.helpers.artifactMetadata(guarded_loaded.met));
end

function test_writemet_yearly_splits_years(testCase)
   % Raw yearly slices independently record local source/support provenance.

   met = makeSyntheticMet(datetime(2015, 12, 31), 48);   % spans 2015/2016

   filenames = icemodel.forcing.helpers.writemet(met, "tst", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");

   testCase.verifyEqual(numel(filenames), 2);
   testCase.verifyTrue(endsWith(filenames(1), "met_tst_src_2015_15m.mat"));
   testCase.verifyTrue(endsWith(filenames(2), "met_tst_src_2016_15m.mat"));

   first = load(filenames(1), 'met', 'artifact_metadata');
   second = load(filenames(2), 'met', 'artifact_metadata');
   testCase.verifyEqual(unique(year(first.met.Time)), 2015);
   testCase.verifyEqual(unique(year(second.met.Time)), 2016);
   testCase.verifyEqual(seconds(median(diff(second.met.Time))), 900);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_row_count, 24);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_source_row_count, 24);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 1));
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 2));
end

function test_writemet_yearly_localizes_guarded_source(testCase)
   % Public build outputs may already be guarded 15-minute data. Exact compact
   % summaries localize lineage without treating derived gap blocks as source.
   source = makeSyntheticMet(datetime(2015, 12, 31), 48);
   source.swd(5) = NaN;
   guarded = icemodel.forcing.helpers.resampleMetTimestep(source, "15m");

   filenames = icemodel.forcing.helpers.writemet(guarded, "guarded", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");

   first = load(filenames(1), 'met', 'artifact_metadata');
   second = load(filenames(2), 'met', 'artifact_metadata');
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_row_count, 24);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_source_row_count, 24);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_missing_counts.swd, 1);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_expected_missing_counts.swd, 4);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_expected_missing_counts.swd, 0);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 1));
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 2));
end

function test_writemet_yearly_handles_two_row_cross_year_source(testCase)
   % A two-row hourly parent split at midnight leaves one row in each year;
   % full-source derivation defines each row's exact one-hour support first.
   source = makeSyntheticMet(datetime(2015, 12, 31, 23, 0, 0), 2);

   filenames = icemodel.forcing.helpers.writemet(source, "boundary", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");

   first = load(filenames(1), 'met', 'artifact_metadata');
   second = load(filenames(2), 'met', 'artifact_metadata');
   testCase.verifyEqual(height(first.met), 4);
   testCase.verifyEqual(height(second.met), 4);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_row_count, 1);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 1));
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 1, 1, 0, 0));
end

function test_writemet_yearly_handles_two_row_three_hour_slice(testCase)
   % Full-source derivation proves 3 h before the Dec31 output is sliced.
   hourly = makeSyntheticMet(datetime(2015, 12, 31, 18, 0, 0), 13);
   source = hourly(1:3:end, :);

   filenames = icemodel.forcing.helpers.writemet(source, "three_hour", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");

   first = load(filenames(1), 'met', 'artifact_metadata');
   second = load(filenames(2), 'met', 'artifact_metadata');
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_row_count, 2);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_source_row_count, 3);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_cadence_seconds, 10800);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 1));
end

function test_writemet_yearly_rejects_unprovable_guarded_source(testCase)
   % Missing summary/cadence metadata or a changed derived row is unprovable.
   source = makeSyntheticMet(datetime(2015, 12, 31), 48);
   guarded = icemodel.forcing.helpers.resampleMetTimestep(source, "15m");
   missing_cadence = guarded;
   missing_cadence.Properties.UserData = rmfield( ...
      missing_cadence.Properties.UserData, ...
      'met_resample_source_cadence_seconds');
   invalid_root = fullfile(testCase.TestData.outdir, 'invalid-guarded-root');

   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      missing_cadence, "missing", "src", ...
      outdir=invalid_root, naming="yearly"), ...
      'icemodel:forcing:writemet:guardedYearlyProvenanceMissing');
   testCase.verifyFalse(isfolder(invalid_root));

   missing_summary = guarded;
   missing_summary.Properties.UserData = rmfield( ...
      missing_summary.Properties.UserData, ...
      'met_resample_yearly_summaries');
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      missing_summary, "missing_summary", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly"), ...
      'icemodel:forcing:writemet:guardedYearlyProvenanceMissing');

   malformed = guarded;
   summaries = malformed.Properties.UserData.met_resample_yearly_summaries;
   summaries(1).support_end_exclusive = ...
      summaries(1).support_end_exclusive + minutes(15);
   malformed.Properties.UserData.met_resample_yearly_summaries = summaries;
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      malformed, "malformed", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly"), ...
      'icemodel:forcing:writemet:guardedYearlyProvenanceInvalid');

   deleted_row = guarded;
   deleted_row(2, :) = [];
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.resampleMetTimestep(deleted_row, "15m"), ...
      'icemodel:forcing:resampleMetTimestep:guardedInputInvalid');
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      deleted_row, "deleted", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly"), ...
      'icemodel:forcing:resampleMetTimestep:guardedInputInvalid');

   guarded.swd(2) = guarded.swd(2) + 1;
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      guarded, "changed", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly"), ...
      'icemodel:forcing:writemet:guardedSourceNotConstant');
end

function test_writemet_yearly_native_guarded_localizes_provenance(testCase)
   % Explicit native output keeps supplied guarded rows and replaces only the
   % parent aggregate facts with exact per-year lineage.
   source = makeSyntheticMet(datetime(2015, 12, 31), 48);
   guarded = icemodel.forcing.helpers.resampleMetTimestep(source, "15m");

   filenames = icemodel.forcing.helpers.writemet(guarded, "native", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly", dt_out="");

   first = load(filenames(1), 'met', 'artifact_metadata');
   second = load(filenames(2), 'met', 'artifact_metadata');
   testCase.verifyEqual(string(first.artifact_metadata.met_resample_policy), ...
      "interval_start_zero_order_hold");
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_row_count, 24);
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_support_end_exclusive, ...
      datetime(2016, 1, 1));
   testCase.verifyEqual(numel( ...
      first.artifact_metadata.met_resample_yearly_summaries), 1);
   testCase.verifyEqual(unique(year(first.met.Time)), 2015);
   testCase.verifyEqual(unique(year(second.met.Time)), 2016);
end

function test_writemet_yearly_preserves_cross_year_outage(testCase)
   % A Dec31 23:00 -> Jan1 02:00 outage stays missing in the new-year file;
   % all-NaN derived blocks are not reclassified as native source rows.
   source = makeSyntheticMet(datetime(2015, 12, 31, 22, 0, 0), 6);
   source(3:4, :) = [];
   source{end, :} = nan(1, width(source));

   filenames = icemodel.forcing.helpers.writemet(source, "outage", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");
   first = load(filenames(1), 'met', 'artifact_metadata');
   second = load(filenames(2), 'met', 'artifact_metadata');
   outage = second.met.Time < datetime(2016, 1, 1, 2, 0, 0);

   testCase.verifyTrue(all(all(isnan(second.met{outage, :}))));
   testCase.verifyEqual( ...
      first.artifact_metadata.met_resample_source_time_gap_count, 0);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_source_time_gap_count, 1);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_source_row_count, 2);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_source_missing_counts.tair, 1);
   testCase.verifyEqual( ...
      second.artifact_metadata.met_resample_expected_missing_counts.tair, 12);

   % A guarded native-cadence request preserves identical rows while replacing
   % only the full-window aggregates with the same exact yearly summaries.
   guarded = icemodel.forcing.helpers.resampleMetTimestep(source, "15m");
   guarded_files = icemodel.forcing.helpers.writemet( ...
      guarded, "outage_guarded", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly", dt_out="");
   guarded_second = load(guarded_files(2), 'met', 'artifact_metadata');
   testCase.verifyEqual(guarded_second.met, second.met);
   testCase.verifyEqual(guarded_second.artifact_metadata, ...
      second.artifact_metadata);
end

function test_writemet_yearly_rejects_fractional_cross_year_step(testCase)
   % A 75-minute Jan1 boundary step is invalid even when each year is regular.
   source = makeSyntheticMet(datetime(2015, 12, 31, 21, 0, 0), 5);
   source.Time(4:5) = source.Time(4:5) + minutes(15);

   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      source, "fractional", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly"), ...
      'icemodel:forcing:resampleMetTimestep:irregularSourceTime');
end

function test_writemet_repeat_is_noop_until_overwrite(testCase)
   % Ordinary repeated writes preserve bytes; overwrite=true replaces them.
   met = makeSyntheticMet(datetime(2016, 1, 1), 48);
   filename = icemodel.forcing.helpers.writemet(met, "tst", "src", ...
      outdir=testCase.TestData.outdir);
   original = fileBytes(filename);

   changed = met;
   changed.tair(:) = changed.tair + 10;
   narrower = changed(5:20, :);
   repeated = icemodel.forcing.helpers.writemet(narrower, "tst", "src", ...
      outdir=testCase.TestData.outdir);

   testCase.verifyEqual(repeated, filename);
   testCase.verifyEqual(fileBytes(filename), original);
   testCase.verifyNumElements(dir(fullfile(testCase.TestData.outdir, ...
      'src', '*.mat')), 1);
   lastwarn("");
   evalc(['icemodel.forcing.helpers.writemet(changed, "tst", "src", ' ...
      'outdir=testCase.TestData.outdir, overwrite=true);']);
   [~, warning_id] = lastwarn;
   testCase.verifyEqual(string(warning_id), ...
      "icemodel:forcing:writemet:overwrite");
   testCase.verifyFalse(isequal(fileBytes(filename), original));
end

function test_writemet_rejects_identity_collision_and_scopes_pruning(testCase)
   % Same-name provenance conflicts require explicit overwrite, and a wider
   % artifact sampled by another method cannot prune the retained shorter file.
   met = makeSyntheticMet(datetime(2016, 1, 1), 96);
   met.Properties.UserData = struct('method', 'nearest', ...
      'lat_wgs84', 67, 'lon_wgs84', -49.5);
   runtime_root = fullfile(testCase.TestData.outdir, 'identity_runtime');
   met_root = fullfile(runtime_root, 'met');
   nearest_short = icemodel.forcing.helpers.writemet( ...
      met(25:48, :), "identity", "src", outdir=met_root);

   natural_short = met(25:48, :);
   natural_short.Properties.UserData.method = 'natural';
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      natural_short, "identity", "src", outdir=met_root), ...
      'icemodel:forcing:writemet:identityConflict');

   natural_wide = met;
   natural_wide.Properties.UserData.method = 'natural';
   natural_wide_file = icemodel.forcing.helpers.writemet( ...
      natural_wide, "identity", "src", outdir=met_root);
   testCase.verifyTrue(isfile(nearest_short));
   testCase.verifyTrue(isfile(natural_wide_file));

   % Runtime resolves the exact narrow name before the compatible broad file.
   % The writer must reject that exact identity collision rather than return the
   % broad artifact and leave runtime pointed at stale provenance.
   runtime_opts = struct('sitename', 'identity', 'forcings', 'src', ...
      'simyears', 2016, 'dt', 900, ...
      'startdate', natural_short.Time(1), ...
      'enddate', natural_short.Time(end), 'pathinput', runtime_root);
   runtime_names = icemodel.createMetFileNames(runtime_opts);
   [~, exact_stem, exact_ext] = fileparts(nearest_short);
   testCase.verifyEqual(string(runtime_names{1}), exact_stem + exact_ext);
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      natural_short, "identity", "src", outdir=met_root), ...
      'icemodel:forcing:writemet:identityConflict');

   % Explicit overwrite repairs the exact runtime target and returns that same
   % path, preserving the established destructive boundary.
   repaired = icemodel.forcing.helpers.writemet( ...
      natural_short, "identity", "src", outdir=met_root, overwrite=true);
   testCase.verifyEqual(repaired, nearest_short);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      repaired, natural_short, "met"));

   % The yearly branch applies the same exact-target identity guard.
   icemodel.forcing.helpers.writemet(met(1:24, :), "identity_year", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      natural_wide(1:24, :), "identity_year", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly"), ...
      'icemodel:forcing:writemet:identityConflict');
end

function test_writemet_checks_native_producer_identity(testCase)
   % Exact-window reuse is safe only when every concrete native producer fact
   % agrees; a conflict must leave the original met artifact byte-stable.
   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   [fields, matching_values, conflicting_values] = ...
      nativeProducerIdentityCases();
   metadata = struct();
   for k = 1:numel(fields)
      metadata.(char(fields(k))) = matching_values(k);
   end
   met.Properties.UserData = metadata;

   for k = 1:numel(fields)
      % A distinct site name isolates each exact-field decision while keeping
      % the artifact window, source, cadence, and payload identical.
      site = "producer_met_" + string(k);
      filename = icemodel.forcing.helpers.writemet( ...
         met, site, "src", outdir=testCase.TestData.outdir);
      original = fileBytes(filename);
      repeated = icemodel.forcing.helpers.writemet( ...
         met, site, "src", outdir=testCase.TestData.outdir);
      testCase.verifyEqual(repeated, filename);
      testCase.verifyEqual(fileBytes(filename), original);

      conflicting = met;
      name = char(fields(k));
      conflicting.Properties.UserData.(name) = conflicting_values(k);
      testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
         conflicting, site, "src", outdir=testCase.TestData.outdir), ...
         'icemodel:forcing:writemet:identityConflict');
      testCase.verifyEqual(fileBytes(filename), original);
   end
end

function test_writemet_rejects_mislabeled_cadence_and_scopes_pruning(testCase)
   % Filename suffixes are not cadence proof: forged hourly payloads named 15m
   % must not be reused, accepted at an exact target, or pruned by a wider write.
   source = makeSyntheticMet( ...
      datetime(2016, 1, 1, TimeZone="UTC"), 96);
   forged_short = icemodel.forcing.helpers.writemet( ...
      source(25:48, :), "forged", "src", ...
      outdir=testCase.TestData.outdir);
   met = source(25:48, :);
   artifact_metadata = ...
      icemodel.forcing.helpers.artifactMetadata(met);
   save(forged_short, 'met', 'artifact_metadata');

   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      source(25:48, :), "forged", "src", ...
      outdir=testCase.TestData.outdir), ...
      'icemodel:forcing:writemet:cadenceConflict');
   wide_file = icemodel.forcing.helpers.writemet( ...
      source, "forged", "src", outdir=testCase.TestData.outdir);
   testCase.verifyTrue(isfile(forged_short));
   testCase.verifyTrue(isfile(wide_file));

   % Yearly exact targets apply the same content-cadence proof.
   forged_year = icemodel.forcing.helpers.writemet( ...
      source(1:24, :), "forged_year", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");
   met = source(1:24, :);
   artifact_metadata = ...
      icemodel.forcing.helpers.artifactMetadata(met);
   save(forged_year, 'met', 'artifact_metadata');
   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      source(1:24, :), "forged_year", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly"), ...
      'icemodel:forcing:writemet:cadenceConflict');

   % A correctly sampled legacy artifact without top-level metadata remains
   % reusable after its actual timetable axis proves the requested cadence.
   legacy_file = icemodel.forcing.helpers.writemet( ...
      source(1:48, :), "legacy_met", "src", ...
      outdir=testCase.TestData.outdir);
   saved = load(legacy_file, 'met');
   met = saved.met;
   save(legacy_file, 'met');
   selected = icemodel.forcing.helpers.writemet( ...
      source(13:24, :), "legacy_met", "src", ...
      outdir=testCase.TestData.outdir);
   testCase.verifyEqual(selected, legacy_file);
end

function test_writemet_wider_window_prunes_only_contained_same_cadence(testCase)
   % A wider successful write removes only the shorter same-class artifact.
   broad = makeSyntheticMet(datetime(2016, 1, 1), 96);
   narrow = broad(25:48, :);

   narrow_file = icemodel.forcing.helpers.writemet( ...
      narrow, "tst", "src", outdir=testCase.TestData.outdir);
   unrelated_file = icemodel.forcing.helpers.writemet( ...
      narrow, "other", "src", outdir=testCase.TestData.outdir);
   native_file = icemodel.forcing.helpers.writemet( ...
      narrow, "tst", "src", outdir=testCase.TestData.outdir, dt_out="");
   unrelated_bytes = fileBytes(unrelated_file);
   native_bytes = fileBytes(native_file);

   % The warning makes the deliberate stale-window cleanup visible to callers.
   lastwarn("");
   evalc(['wide_file = icemodel.forcing.helpers.writemet(' ...
      'broad, "tst", "src", outdir=testCase.TestData.outdir);']);
   [~, warning_id] = lastwarn;

   testCase.verifyEqual(string(warning_id), ...
      "icemodel:forcing:pruneSupersededWindowFiles:removed");
   testCase.verifyTrue(isfile(wide_file));
   testCase.verifyFalse(isfile(narrow_file));
   testCase.verifyTrue(isfile(unrelated_file));
   testCase.verifyTrue(isfile(native_file));
   testCase.verifyEqual(fileBytes(unrelated_file), unrelated_bytes);
   testCase.verifyEqual(fileBytes(native_file), native_bytes);
   testCase.verifyNumElements(dir(fullfile(testCase.TestData.outdir, ...
      'src', 'met_tst_src_*_15m.mat')), 1);
end

function test_writemet_allows_extra_diagnostic_columns(testCase)
   % Metadata stamping should annotate known met channels without rejecting
   % caller-specific diagnostics that validatemet otherwise permits.
   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met.custom_diagnostic = (1:height(met))';

   filenames = icemodel.forcing.helpers.writemet(met, "tst", "src", ...
      outdir=testCase.TestData.outdir);

   loaded = load(filenames, 'met');
   testCase.verifyTrue(ismember("custom_diagnostic", ...
      string(loaded.met.Properties.VariableNames)));
   testCase.verifyEqual(string(loaded.met.Properties.VariableUnits( ...
      string(loaded.met.Properties.VariableNames) == "custom_diagnostic")), "");
   testCase.verifyEqual(string(loaded.met.Properties.VariableUnits( ...
      string(loaded.met.Properties.VariableNames) == "tair")), "K");
end

function test_writemet_rejects_noncanonical_ppt_units(testCase)
   % Existing ppt units must be validated before metadata stamping fills the
   % canonical unit labels or any output/source directory is created.
   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met.Properties.VariableUnits = repmat({''}, 1, width(met));
   met.Properties.VariableUnits{string(met.Properties.VariableNames) == "ppt"} = ...
      'mm hr-1';
   invalid_root = fullfile(testCase.TestData.outdir, 'invalid-met-root');

   testCase.verifyError(@() icemodel.forcing.helpers.writemet( ...
      met, "tst", "invalid", outdir=invalid_root), ...
      'icemodel:forcing:validatemet:pptUnit');
   testCase.verifyFalse(isfolder(invalid_root));
end

%% writeuserdata

function test_writeuserdata_per_year_with_metadata(testCase)
   % Default naming="yearly": Data files split per year, keep the Data variable
   % name, and carry the location CustomProperties.

   Data = makeSyntheticData(datetime(2015, 12, 31), 48);

   filenames = icemodel.forcing.helpers.writeuserdata(Data, "tst", "src", ...
      outdir=testCase.TestData.outdir);

   testCase.verifyEqual(numel(filenames), 2);
   testCase.verifyTrue(endsWith(filenames(1), "tst_src_2015.mat"));
   testCase.verifyTrue(endsWith(filenames(2), "tst_src_2016.mat"));

   loaded = load(filenames(1), 'Data', 'artifact_metadata');
   testCase.verifyEqual(unique(year(loaded.Data.Time)), 2015);
   testCase.verifyEqual(loaded.Data.Properties.CustomProperties.Lat, 67.0);
   testCase.verifyEqual(string(loaded.Data.Properties.VariableUnits( ...
      string(loaded.Data.Properties.VariableNames) == "tair")), "K");
   testCase.verifyEqual(loaded.Data.Properties.UserData, ...
      loaded.artifact_metadata);
   testCase.verifyEqual(loaded.artifact_metadata.artifact_cadence_seconds, 3600);
end

function test_writeuserdata_window_is_full_period(testCase)
   % naming="window": ONE full-period file named
   % <site>_<source>_<YYYYMMDD>_<YYYYMMDD>.mat that preserves the entire span
   % (not split per year) and keeps the location CustomProperties.

   Data = makeSyntheticData(datetime(2015, 12, 31), 48);

   filenames = icemodel.forcing.helpers.writeuserdata(Data, "tst", "src", ...
      outdir=testCase.TestData.outdir, naming="window");

   testCase.verifyEqual(numel(filenames), 1);
   expected = sprintf('tst_src_%s_%s.mat', ...
      char(min(Data.Time), 'yyyyMMdd'), char(max(Data.Time), 'yyyyMMdd'));
   testCase.verifyTrue(endsWith(filenames(1), expected));

   loaded = load(filenames(1), 'Data');
   % Full span retained in the single file (both calendar years present).
   testCase.verifyEqual(unique(year(loaded.Data.Time))', [2015, 2016]);
   testCase.verifyEqual(height(loaded.Data), height(Data));
   testCase.verifyEqual(loaded.Data.Properties.CustomProperties.Lat, 67.0);
end

function test_writeuserdata_repeat_preserves_native_data_until_overwrite(testCase)
   % Hourly userdata remains unchanged and ordinary repeats preserve bytes.
   Data = makeSyntheticData(datetime(2016, 1, 1), 48);
   filename = icemodel.forcing.helpers.writeuserdata(Data, "tst", "src", ...
      outdir=testCase.TestData.outdir, naming="window");
   original = fileBytes(filename);

   changed = Data;
   changed.tair(:) = changed.tair + 10;
   narrower = changed(5:20, :);
   repeated = icemodel.forcing.helpers.writeuserdata( ...
      narrower, "tst", "src", outdir=testCase.TestData.outdir, ...
      naming="window");

   testCase.verifyEqual(repeated, filename);
   testCase.verifyEqual(fileBytes(filename), original);
   testCase.verifyNumElements(dir(fullfile(testCase.TestData.outdir, ...
      'src', '*.mat')), 1);
   loaded = load(filename, 'Data');
   testCase.verifyEqual(seconds(median(diff(loaded.Data.Time))), 3600);
   lastwarn("");
   evalc(['icemodel.forcing.helpers.writeuserdata(changed, "tst", "src", ' ...
      'outdir=testCase.TestData.outdir, naming="window", overwrite=true);']);
   [~, warning_id] = lastwarn;
   testCase.verifyEqual(string(warning_id), ...
      "icemodel:forcing:writeuserdata:overwrite");
   testCase.verifyFalse(isequal(fileBytes(filename), original));
end

function test_writeuserdata_rejects_identity_collision_and_scopes_pruning(testCase)
   % Userdata reuse/pruning is limited to artifacts with compatible concrete
   % sampling and point identity, for both window and yearly naming.
   Data = makeSyntheticData(datetime(2016, 1, 1), 96);
   Data.Properties.UserData = struct('method', 'nearest', ...
      'source_id', 'mar3.11');
   nearest_short = icemodel.forcing.helpers.writeuserdata( ...
      Data(25:48, :), "identity_data", "src", ...
      outdir=testCase.TestData.outdir, naming="window");

   natural_short = Data(25:48, :);
   natural_short.Properties.UserData.method = 'natural';
   testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
      natural_short, "identity_data", "src", ...
      outdir=testCase.TestData.outdir, naming="window"), ...
      'icemodel:forcing:writeuserdata:identityConflict');

   natural_wide = Data;
   natural_wide.Properties.UserData.method = 'natural';
   natural_wide_file = icemodel.forcing.helpers.writeuserdata( ...
      natural_wide, "identity_data", "src", ...
      outdir=testCase.TestData.outdir, naming="window");
   testCase.verifyTrue(isfile(nearest_short));
   testCase.verifyTrue(isfile(natural_wide_file));

   % A compatible broad file must not bypass validation of the conflicting
   % exact target, even though normal userdata runtime selection prefers broad
   % enclosing support.
   testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
      natural_short, "identity_data", "src", ...
      outdir=testCase.TestData.outdir, naming="window"), ...
      'icemodel:forcing:writeuserdata:identityConflict');

   % Point conflicts are concrete even when the sampling method agrees.
   moved = Data(25:48, :);
   moved.Properties.CustomProperties.Lat = 68;
   testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
      moved, "identity_data", "src", outdir=testCase.TestData.outdir, ...
      naming="window"), ...
      'icemodel:forcing:writeuserdata:identityConflict');

   icemodel.forcing.helpers.writeuserdata(Data(1:24, :), ...
      "identity_data_year", "src", outdir=testCase.TestData.outdir);
   testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
      natural_wide(1:24, :), "identity_data_year", "src", ...
      outdir=testCase.TestData.outdir), ...
      'icemodel:forcing:writeuserdata:identityConflict');
end

function test_writeuserdata_checks_native_producer_identity(testCase)
   % Exact-window userdata reuse follows the same native producer identity
   % boundary as met and never relabels an existing payload on conflict.
   Data = makeSyntheticData(datetime(2016, 1, 1), 24);
   [fields, matching_values, conflicting_values] = ...
      nativeProducerIdentityCases();
   metadata = struct();
   for k = 1:numel(fields)
      metadata.(char(fields(k))) = matching_values(k);
   end
   Data.Properties.UserData = metadata;

   for k = 1:numel(fields)
      % Keep all non-producer identity and exact-window naming inputs stable so
      % the tested field is the only possible reuse boundary.
      site = "producer_data_" + string(k);
      filename = icemodel.forcing.helpers.writeuserdata( ...
         Data, site, "src", outdir=testCase.TestData.outdir, ...
         naming="window");
      original = fileBytes(filename);
      repeated = icemodel.forcing.helpers.writeuserdata( ...
         Data, site, "src", outdir=testCase.TestData.outdir, ...
         naming="window");
      testCase.verifyEqual(repeated, filename);
      testCase.verifyEqual(fileBytes(filename), original);

      conflicting = Data;
      name = char(fields(k));
      conflicting.Properties.UserData.(name) = conflicting_values(k);
      testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
         conflicting, site, "src", outdir=testCase.TestData.outdir, ...
         naming="window"), ...
         'icemodel:forcing:writeuserdata:identityConflict');
      testCase.verifyEqual(fileBytes(filename), original);
   end
end

function test_writeuserdata_defaults_hourly_and_supports_native_override(testCase)
   % Finer source data become hourly only at the shared writer boundary.
   source = makeSyntheticData(datetime(2016, 1, 1), 4);
   half_hour = (source.Time(1):minutes(30):source.Time(end))';
   Data = retime(source, half_hour, 'linear');
   direction = repmat([350; 10], ceil(height(Data) / 2), 1);
   Data.wdir = direction(1:height(Data));

   hourly_file = icemodel.forcing.helpers.writeuserdata( ...
      Data, "hourly", "src", outdir=testCase.TestData.outdir, ...
      naming="window");
   hourly = load(hourly_file, 'Data');
   testCase.verifyEqual(seconds(median(diff(hourly.Data.Time))), 3600);
   testCase.verifyEqual(string(hourly.Data.Properties.UserData. ...
      userdata_resample_policy), "hourly_mean");
   testCase.verifyLessThan(min(abs([hourly.Data.wdir(1), ...
      hourly.Data.wdir(1) - 360])), 1e-10);

   native_file = icemodel.forcing.helpers.writeuserdata( ...
      Data, "native", "src", outdir=testCase.TestData.outdir, ...
      naming="window", dt_out="");
   native = load(native_file, 'Data');
   testCase.verifyEqual(seconds(median(diff(native.Data.Time))), 1800);
   testCase.verifyEqual(string(native.Data.Properties.UserData. ...
      userdata_resample_policy), "explicit_native");
end

function test_writeuserdata_hourly_enclosing_file_is_not_reused_for_native(testCase)
   % A broader hourly artifact cannot satisfy an explicit 30-minute request.
   hourly = makeSyntheticData(datetime(2016, 1, 1), 97);
   native_time = (hourly.Time(1):minutes(30):hourly.Time(end))';
   native = retime(hourly, native_time, 'linear');
   hourly_file = icemodel.forcing.helpers.writeuserdata( ...
      native, "h2n", "src", outdir=testCase.TestData.outdir, ...
      naming="window");

   % The public native override has a cadence-qualified name, so the exact same
   % site/source/window can coexist without reuse or destructive replacement.
   native_file = icemodel.forcing.helpers.writeuserdata( ...
      native, "h2n", "src", outdir=testCase.TestData.outdir, ...
      naming="window", dt_out="");
   testCase.verifyNotEqual(native_file, hourly_file);
   testCase.verifyTrue(endsWith(native_file, "_30m.mat"));

   narrow = native(native.Time >= datetime(2016, 1, 2) ...
      & native.Time < datetime(2016, 1, 3), :);
   selected_native = icemodel.forcing.helpers.writeuserdata( ...
      narrow, "h2n", "src", outdir=testCase.TestData.outdir, ...
      naming="window", dt_out="");
   saved_hourly = load(hourly_file, 'Data');
   saved_native = load(native_file, 'Data');

   testCase.verifyEqual(selected_native, native_file);
   testCase.verifyTrue(isfile(hourly_file));
   testCase.verifyEqual(seconds(median(diff(saved_hourly.Data.Time))), 3600);
   testCase.verifyEqual(seconds(median(diff(saved_native.Data.Time))), 1800);

   % Same-cadence containment remains byte-stable and reuses the native file.
   before = fileBytes(native_file);
   repeated = icemodel.forcing.helpers.writeuserdata( ...
      narrow(5:20, :), "h2n", "src", outdir=testCase.TestData.outdir, ...
      naming="window", dt_out="");
   testCase.verifyEqual(repeated, native_file);
   testCase.verifyEqual(fileBytes(native_file), before);

   % The yearly naming form uses the same cadence identity and cannot collide.
   hourly_year = icemodel.forcing.helpers.writeuserdata( ...
      native, "h2n_year", "src", outdir=testCase.TestData.outdir);
   native_year = icemodel.forcing.helpers.writeuserdata( ...
      native, "h2n_year", "src", outdir=testCase.TestData.outdir, dt_out="");
   testCase.verifyTrue(all(isfile(hourly_year)));
   testCase.verifyTrue(all(isfile(native_year)));
   testCase.verifyTrue(all(endsWith(native_year, "_30m.mat")));
   testCase.verifyEmpty(intersect(hourly_year, native_year));
end

function test_writeuserdata_native_enclosing_file_is_not_reused_for_hourly(testCase)
   % A broader native 30-minute artifact cannot satisfy an hourly request.
   hourly = makeSyntheticData(datetime(2016, 1, 1), 97);
   native_time = (hourly.Time(1):minutes(30):hourly.Time(end))';
   native = retime(hourly, native_time, 'linear');
   native_file = icemodel.forcing.helpers.writeuserdata( ...
      native, "n2h", "src", outdir=testCase.TestData.outdir, ...
      naming="window", dt_out="");
   narrow = native(native.Time >= datetime(2016, 1, 2) ...
      & native.Time < datetime(2016, 1, 3), :);
   hourly_file = icemodel.forcing.helpers.writeuserdata( ...
      narrow, "n2h", "src", outdir=testCase.TestData.outdir, ...
      naming="window");
   saved_native = load(native_file, 'Data');
   saved_hourly = load(hourly_file, 'Data');

   testCase.verifyNotEqual(hourly_file, native_file);
   testCase.verifyEqual(seconds(median(diff(saved_native.Data.Time))), 1800);
   testCase.verifyEqual(seconds(median(diff(saved_hourly.Data.Time))), 3600);

   % Wider same-name cleanup is cadence-scoped: a new hourly window must retain
   % a contained native file even though both use the suffix-less userdata class.
   contained_native = icemodel.forcing.helpers.writeuserdata( ...
      narrow, "prune_cadence", "src", outdir=testCase.TestData.outdir, ...
      naming="window", dt_out="");
   wider_hourly = icemodel.forcing.helpers.writeuserdata( ...
      native, "prune_cadence", "src", outdir=testCase.TestData.outdir, ...
      naming="window");
   testCase.verifyTrue(isfile(contained_native));
   testCase.verifyTrue(isfile(wider_hourly));
end

function test_writeuserdata_legacy_enclosing_file_uses_proven_time_axis(testCase)
   % Legacy files without cadence metadata remain reusable only when their saved
   % timetable proves the requested uniform cadence.
   Data = makeSyntheticData(datetime(2016, 1, 1), 72);
   filename = icemodel.forcing.helpers.writeuserdata( ...
      Data, "legacy", "src", outdir=testCase.TestData.outdir, ...
      naming="window");
   saved = load(filename, 'Data');
   Data = saved.Data;
   metadata = Data.Properties.UserData;
   fields = string(fieldnames(metadata));
   fields = fields(startsWith(fields, "userdata_resample_"));
   metadata = rmfield(metadata, cellstr(fields));
   Data.Properties.UserData = metadata;
   save(filename, 'Data');
   before = fileBytes(filename);

   returned = icemodel.forcing.helpers.writeuserdata( ...
      Data(13:36, :), "legacy", "src", outdir=testCase.TestData.outdir, ...
      naming="window");

   testCase.verifyEqual(returned, filename);
   testCase.verifyEqual(fileBytes(filename), before);
end

function test_writeuserdata_interpolates_coarser_source_to_hourly(testCase)
   % The unconditional hourly default also handles a coarser native product.
   Data = makeSyntheticData(datetime(2016, 1, 1), 7);
   Data = Data(1:2:end, :);
   Data.wdir = [350; 10; 30; 50];

   filename = icemodel.forcing.helpers.writeuserdata( ...
      Data, "coarse", "src", outdir=testCase.TestData.outdir, ...
      naming="window");
   saved = load(filename, 'Data');
   testCase.verifyEqual(seconds(median(diff(saved.Data.Time))), 3600);
   testCase.verifyEqual(height(saved.Data), 7);
   testCase.verifyEqual(string(saved.Data.Properties.UserData. ...
      userdata_resample_policy), "hourly_linear");
   testCase.verifyLessThan(min(abs([saved.Data.wdir(2), ...
      saved.Data.wdir(2) - 360])), 1e-10);
end

function test_writeuserdata_wider_refresh_prunes_only_contained_window(testCase)
   % A successful wider refresh removes its shorter predecessor, not siblings.
   Data = makeSyntheticData(datetime(2016, 1, 1), 96);
   short = Data(25:72, :);
   short_file = icemodel.forcing.helpers.writeuserdata( ...
      short, "wide", "src", outdir=testCase.TestData.outdir, ...
      naming="window");
   sibling_file = icemodel.forcing.helpers.writeuserdata( ...
      short, "sibling", "src", outdir=testCase.TestData.outdir, ...
      naming="window");

   testCase.verifyWarning(@() ...
      icemodel.forcing.helpers.writeuserdata( ...
      Data, "wide", "src", outdir=testCase.TestData.outdir, ...
      naming="window"), ...
      'icemodel:forcing:pruneSupersededWindowFiles:removed');
   wide_file = fullfile(fileparts(short_file), ...
      'wide_src_20160101_20160104.mat');
   testCase.verifyFalse(isfile(short_file));
   testCase.verifyTrue(isfile(wide_file));
   testCase.verifyTrue(isfile(sibling_file));

   % A later narrow request reuses the widest artifact without rewriting it.
   before = fileBytes(wide_file);
   returned = icemodel.forcing.helpers.writeuserdata( ...
      short, "wide", "src", outdir=testCase.TestData.outdir, ...
      naming="window");
   testCase.verifyEqual(returned, string(wide_file));
   testCase.verifyEqual(fileBytes(wide_file), before);
end

function test_findEnclosingWindowFile_prefers_widest_then_latest(testCase)
   % Selection is independent of filesystem directory order.
   folder = fullfile(testCase.TestData.outdir, 'windows');
   mkdir(folder)
   names = ["x_20140101_20160101.mat", ...
      "x_20150101_20170101.mat", "x_20130101_20180101.mat"];
   for name = names
      fclose(fopen(fullfile(folder, name), 'w'));
   end
   q1 = datetime(2015, 6, 1, 'TimeZone', 'UTC');
   q2 = datetime(2015, 7, 1, 'TimeZone', 'UTC');
   selected = icemodel.forcing.helpers.findEnclosingWindowFile( ...
      folder, "x", ".mat", q1, q2);
   testCase.verifyEqual(selected, "x_20130101_20180101.mat");

   % Equal-width containing windows prefer the one with the latest end date.
   delete(fullfile(folder, "x_20130101_20180101.mat"));
   selected = icemodel.forcing.helpers.findEnclosingWindowFile( ...
      folder, "x", ".mat", q1, q2);
   testCase.verifyEqual(selected, "x_20150101_20170101.mat");
end

function test_pruneSupersededWindowFiles_keeps_noncontained_files(testCase)
   % Cleanup is restricted to strictly contained files in one naming class.
   folder = fullfile(testCase.TestData.outdir, 'prune');
   mkdir(folder)
   names = ["x_20150101_20151231.mat", ...
      "x_20150301_20150401.mat", "x_20140101_20150301.mat", ...
      "x_20140101_20160101.mat", "x_bad.mat"];
   for name = names
      fclose(fopen(fullfile(folder, name), 'w'));
   end
   new_file = fullfile(folder, names(1));

   testCase.verifyWarning(@() ...
      icemodel.forcing.helpers.pruneSupersededWindowFiles( ...
      new_file, "x", ".mat"), ...
      'icemodel:forcing:pruneSupersededWindowFiles:removed');
   testCase.verifyFalse(isfile(fullfile(folder, names(2))));
   testCase.verifyTrue(isfile(fullfile(folder, names(3))));
   testCase.verifyTrue(isfile(fullfile(folder, names(4))));
   testCase.verifyTrue(isfile(fullfile(folder, names(5))));
   testCase.verifyEmpty( ...
      icemodel.forcing.helpers.pruneSupersededWindowFiles( ...
      new_file, "x", ".mat"));
   testCase.verifyEmpty( ...
      icemodel.forcing.helpers.pruneSupersededWindowFiles( ...
      fullfile(folder, "invalid.mat"), "x", ".mat"));
   testCase.verifyEmpty( ...
      icemodel.forcing.helpers.pruneSupersededWindowFiles( ...
      fullfile(folder, 'missing', names(1)), "x", ".mat"));
end

function test_writeuserdata_rejects_missing_metadata(testCase)
   % The CustomProperties location-metadata contract is enforced.

   Data = makeSyntheticData(datetime(2016, 1, 1), 24);
   Data = removeCustomProperties(Data);

   testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
      Data, "tst", "src", outdir=testCase.TestData.outdir), ...
      'icemodel:forcing:writeuserdata:missingMetadata');

   Data = makeSyntheticData(datetime(2016, 1, 1), 24);
   Data.Properties.UserData = "invalid";
   testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
      Data, "tst", "src", outdir=testCase.TestData.outdir), ...
      'icemodel:forcing:writeuserdata:badMetadata');
end

%% psnProjection

function test_psnProjection_matches_legacy_projsipsn(testCase)
   % psnProjection (EPSG:3413) must reproduce the legacy projsipsn.mat
   % forward projection. Skips when matfunclib (which carries
   % projsipsn.mat) is not on the path.

   legacy_file = which('projsipsn.mat');
   testCase.assumeTrue(~isempty(legacy_file), ...
      'projsipsn.mat not on path (matfunclib not loaded)');

   legacy = load(legacy_file).('projsipsn');
   lat = [61.1; 67.0; 72.5];
   lon = [-48.3; -49.5; -40.0];

   [x_old, y_old] = projfwd(legacy, lat, lon);
   [x_new, y_new] = projfwd(icemodel.forcing.helpers.psnProjection(), lat, lon);

   testCase.verifyEqual(x_new, x_old, 'AbsTol', 0.01);
   testCase.verifyEqual(y_new, y_old, 'AbsTol', 0.01);
end

%% sourceAlbedo

function test_sourceAlbedo_rejects_low_insolation_ratios(testCase)
   % Low-sun and nonpositive-reflection ratios are sensor noise, while the
   % first supported sample at 10 W m-2 remains observational.
   [albedo, counts] = icemodel.forcing.helpers.sourceAlbedo( ...
      [0; 5; 10; 20; 20; NaN], [0; 4; 8; 0; 1; 1]);

   testCase.verifyEqual(albedo, [NaN; NaN; 0.8; NaN; 0.05; NaN]);
   testCase.verifyEqual(counts.low_light, 1);
   testCase.verifyEqual(counts.nonpositive_swu, 1);
   testCase.verifyEqual(counts.below_minimum, 0);
   testCase.verifyEqual(counts.total, 2);
end

function test_sourceAlbedo_supports_source_specific_minimum(testCase)
   % Accumulation-site adapters may reject isolated physically implausible
   % ratios without changing the shared 10 W m-2 low-light threshold.
   [albedo, counts] = icemodel.forcing.helpers.sourceAlbedo( ...
      [100; 100; 100], [20; 30; 50], minimum=0.3);

   testCase.verifyEqual(albedo, [NaN; 0.3; 0.5], 'AbsTol', 1e-12);
   testCase.verifyEqual(counts.below_minimum, 1);
   testCase.verifyEqual(counts.total, 1);
end

function test_sourceAlbedo_rejects_low_solar_elevation(testCase)
   % Location-aware radiometer QC rejects a finite nighttime ratio while
   % preserving the same ratio at solar noon; callers may override the angle.
   Time = datetime(2021, 6, 21, [0; 12], 0, 0, 'TimeZone', 'UTC');
   swd = [100; 100];
   swu = [50; 50];

   [albedo, counts] = icemodel.forcing.helpers.sourceAlbedo( ...
      swd, swu, Time=Time, latitude=0, longitude=0);
   unfiltered = icemodel.forcing.helpers.sourceAlbedo( ...
      swd, swu, Time=Time, latitude=0, longitude=0, ...
      minimum_solar_elevation=-90);

   testCase.verifyEqual(albedo, [NaN; 0.5], 'AbsTol', 1e-12);
   testCase.verifyEqual(unfiltered, [0.5; 0.5], 'AbsTol', 1e-12);
   testCase.verifyEqual(counts.low_solar_elevation, 1);
   testCase.verifyEqual(counts.total, 1);
end

function test_sourceAlbedo_requires_complete_matching_solar_geometry(testCase)
   % Partial geometry and mismatched timestamp arrays must fail explicitly.
   Time = datetime(2021, 6, 21, 'TimeZone', 'UTC');

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.sourceAlbedo(100, 50, Time=Time), ...
      'icemodel:forcing:helpers:sourceAlbedo:incompleteSolarGeometry');
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.sourceAlbedo([100; 100], [50; 50], ...
      Time=Time, latitude=0, longitude=0), ...
      'icemodel:forcing:helpers:sourceAlbedo:timeSizeMismatch');
end

%% dailyAlbedoAnomalyFlags

function test_dailyAlbedoAnomalyFlags_expands_only_seeded_episode(testCase)
   % One conservative collapse seeds the neighboring weaker transient days,
   % while raw radiation and the caller's row-vector shape remain unchanged.
   [Time, swd, swu, dates] = transientAlbedoFixture();
   swd_before = swd;
   swu_before = swu;

   [flags, report] = ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      Time.', swd.', swu.');

   expected_days = dates([44 46 48]);
   flagged_days = unique(dateshift(Time(flags(:)), 'start', 'day'));
   testCase.verifyEqual(flagged_days, expected_days);
   testCase.verifySize(flags, [1 numel(Time)]);
   testCase.verifyEqual(report.seed_day_count, 1);
   testCase.verifyEqual(report.episode_day_count, 3);
   testCase.verifyEqual(report.flagged_row_count, 72);
   testCase.verifyEqual(report.seed_dates, ...
      string(dates(46), 'yyyy-MM-dd'));
   testCase.verifyEqual(report.episode_dates, ...
      string(expected_days, 'yyyy-MM-dd'));
   testCase.verifyEqual(swd, swd_before);
   testCase.verifyEqual(swu, swu_before);
   testCase.verifyTrue(contains(report.policy, "two-sided drop"));
end

function test_dailyAlbedoAnomalyFlags_rejects_broken_native_grids(testCase)
   % Duplicate-plus-missing, off-grid, missing, and reversed timestamp support
   % must fail closed instead of passing a row-count-only completeness test.
   [Time, swd, swu, dates] = transientAlbedoFixture();
   seed_rows = find(dateshift(Time, 'start', 'day') == dates(46));

   duplicate_time = Time;
   duplicate_time(seed_rows(9)) = duplicate_time(seed_rows(8));
   duplicate = icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      duplicate_time, swd, swu);
   testCase.verifyFalse(any(duplicate));

   off_grid_time = Time;
   off_grid_time(seed_rows(9)) = off_grid_time(seed_rows(9)) + minutes(5);
   off_grid = icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      off_grid_time, swd, swu);
   testCase.verifyFalse(any(off_grid));

   keep = true(size(Time));
   keep(seed_rows(9)) = false;
   missing = icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      Time(keep), swd(keep), swu(keep));
   testCase.verifyFalse(any(missing));

   reversed = icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      flipud(Time), flipud(swd), flipud(swu));
   testCase.verifyFalse(any(reversed));
end

function test_dailyAlbedoAnomalyFlags_rejects_persistent_dark_surface(testCase)
   % A sustained dark-ice season is not a recovered sensor-collapse episode,
   % even though its albedo is lower than the conservative absolute cap.
   [Time, swd, ~, dates] = transientAlbedoFixture();
   alpha = 0.85 * ones(size(swd));
   dark_days = dates(20:75);
   alpha(ismember(dateshift(Time, 'start', 'day'), dark_days)) = 0.4;
   swu = alpha .* swd;

   [flags, report] = ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags(Time, swd, swu);

   testCase.verifyFalse(any(flags));
   testCase.verifyEqual(report.seed_day_count, 0);
   testCase.verifyEqual(report.episode_day_count, 0);
end

function test_dailyAlbedoAnomalyFlags_enforces_support_gates(testCase)
   % Low incident energy, an incomplete day, and insufficient reflected-energy
   % coverage each prevent an otherwise abrupt ratio collapse from seeding.
   [Time, swd, swu, dates] = transientAlbedoFixture();
   anomaly_rows = dateshift(Time, 'start', 'day') == dates(46);

   low_energy_swd = swd;
   low_energy_swd(anomaly_rows) = 50;
   low_energy_swu = swu;
   low_energy_swu(anomaly_rows) = 0.4 .* low_energy_swd(anomaly_rows);
   low_energy = icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      Time, low_energy_swd, low_energy_swu);
   testCase.verifyFalse(any(low_energy));

   incomplete = ~findFirst(anomaly_rows);
   incomplete_flags = ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      Time(incomplete), swd(incomplete), swu(incomplete));
   testCase.verifyFalse(any(incomplete_flags));

   low_coverage_swu = swu;
   anomaly_indices = find(anomaly_rows);
   low_coverage_swu(anomaly_indices(1:13)) = NaN;
   low_coverage = icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      Time, swd, low_coverage_swu);
   testCase.verifyFalse(any(low_coverage));
end

function test_dailyAlbedoAnomalyFlags_handles_short_and_mismatched_inputs( ...
      testCase)
   % A one-row record has no daily context, while mismatched channels are a
   % caller error rather than an implicit truncation policy.
   Time = datetime(2017, 1, 1, 'TimeZone', 'UTC');
   [flags, report] = ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags(Time, 100, 80);

   testCase.verifyFalse(flags);
   testCase.verifyEqual(report.episode_day_count, 0);
   testCase.verifyEmpty(report.episode_dates);
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      Time, [100; 100], 80), ...
      'icemodel:forcing:helpers:dailyAlbedoAnomalyFlags:sizeMismatch');
end

%% modisAlbedoChannel

function test_modisAlbedoChannel_leaves_missing_years_nan(testCase)
   % MODIS is optional for MAR/MERRA/RACMO staging, so missing coverage for a
   % requested year leaves that year's channel missing instead of failing the
   % parent RCM source leg.

   empty_dir = string(tempname);
   mkdir(empty_dir);
   testCase.addTeardown(@() rmdir(empty_dir, 's'));

   Time = (datetime(2009, 1, 1):hours(1):datetime(2009, 1, 1, 5, 0, 0))';
   modis = icemodel.forcing.helpers.modisAlbedoChannel( ...
      empty_dir, 2009, [67.0 -49.5], "nearest", "conservative", Time);

   testCase.verifySize(modis, [numel(Time), 1]);
   testCase.verifyTrue(all(isnan(modis)));
end

function test_modisAlbedoChannel_errors_on_duplicate_files(testCase)
   % Duplicate files for the same MODIS year are a bad source layout, not a
   % coverage gap, so the helper must refuse to choose one silently.

   modis_dir = string(tempname);
   mkdir(modis_dir);
   testCase.addTeardown(@() rmdir(modis_dir, 's'));
   fclose(fopen(fullfile(modis_dir, 'Greenland_Reflectivity_2009_a.nc'), 'w'));
   fclose(fopen(fullfile(modis_dir, 'Greenland_Reflectivity_2009_b.nc'), 'w'));

   Time = (datetime(2009, 1, 1):hours(1):datetime(2009, 1, 1, 5, 0, 0))';
   testCase.verifyError(@() icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, 2009, [67.0 -49.5], "nearest", "conservative", Time), ...
      'icemodel:forcing:modisAlbedoChannel:ambiguousFile');
end

function test_modisAlbedoChannel_accepts_column_years(testCase)
   % Staging passes unique(year(Time)) as a column vector. The helper must
   % still iterate one year at a time instead of building one impossible
   % multi-year filename pattern.

   % Resolve the scoped root installed by the test bootstrap so provisioned
   % optional forcing data cannot be skipped behind the former data/test path.
   config = icemodel.config('getenv', true);
   modis_dir = string(fullfile(config.ICEMODEL_DATA_PATH, 'forcing', ...
      'geus', 'albedo', 'gris'));
   testCase.assumeTrue(isfolder(modis_dir), ...
      'GEUS MODIS albedo source not available');

   Time = (datetime(2012, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'))';
   modis = icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, 2012, [67.1556, -49.9226], "nearest", ...
      "conservative", Time);
   modis_col = icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, [2011; 2012], [67.1556, -49.9226], "nearest", ...
      "conservative", Time);

   testCase.verifyEqual(modis_col, modis);
   testCase.verifyGreaterThan(nnz(isfinite(modis_col)), 0);
end

%% Local fixture helpers

function [Time, swd, swu, dates] = transientAlbedoFixture()
   %TRANSIENTALBEDOFIXTURE Make a seeded three-day reflected-SW collapse.
   dates = (datetime(2017, 1, 1, 'TimeZone', 'UTC') ...
      + caldays(0:89))';
   Time = (dates(1):hours(1):(dates(end) + hours(23)))';
   swd = 100 * ones(size(Time));
   alpha = 0.85 * ones(size(Time));

   % Two weaker days separated by one calendar day surround one conservative
   % seed, exercising episode expansion without creating a persistent season.
   day_key = dateshift(Time, 'start', 'day');
   alpha(day_key == dates(44)) = 0.72;
   alpha(day_key == dates(46)) = 0.40;
   alpha(day_key == dates(48)) = 0.72;
   swu = alpha .* swd;
end

function selected = findFirst(mask)
   %FINDFIRST Select only the first true element of a logical vector.
   selected = false(size(mask));
   first = find(mask, 1);
   if ~isempty(first)
      selected(first) = true;
   end
end

function met = makeSyntheticMet(t0, nsteps)
   %MAKESYNTHETICMET Hourly timetable satisfying the met contract.
   Time = (t0:hours(1):(t0 + hours(nsteps - 1)))';
   n = numel(Time);
   tair = 265 + 5*sin(2*pi*(1:n)'/24);
   swd = max(0, 400*sin(2*pi*(1:n)'/24));
   lwd = 250*ones(n, 1);
   albedo = 0.8*ones(n, 1);
   wspd = 4 + sin(2*pi*(1:n)'/24);
   rh = 70*ones(n, 1);
   psfc = 1e5*ones(n, 1);
   ppt = zeros(n, 1);
   met = timetable(Time, tair, swd, lwd, albedo, wspd, rh, psfc, ppt);
end

function Data = makeSyntheticData(t0, nsteps)
   %MAKESYNTHETICDATA Data timetable with location CustomProperties.
   Data = makeSyntheticMet(t0, nsteps);
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = -2.1e5;
   Data.Properties.CustomProperties.Y = -2.5e6;
   Data.Properties.CustomProperties.Lat = 67.0;
   Data.Properties.CustomProperties.Lon = -49.5;
   Data.Properties.CustomProperties.Elev = 1270;
   Data.Properties.CustomProperties.Slope = 0.01;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];
end

function Data = removeCustomProperties(Data)
   %REMOVECUSTOMPROPERTIES Rebuild DATA without CustomProperties.
   Data = timetable(Data.Time, Data.tair, 'VariableNames', {'tair'});
end

function [fields, matching, conflicting] = nativeProducerIdentityCases()
   %NATIVEPRODUCERIDENTITYCASES Return exact production metadata spellings.
   fields = ["source_family", "station", "doi", "bundle_doi"];
   matching = ["imau", "S21", "10.1594/PANGAEA.969585", ...
      "10.1594/PANGAEA.971647"];
   conflicting = ["gcnet_vandecrux", "S22", ...
      "10.1594/PANGAEA.969629", "10.1594/PANGAEA.970127"];
end

function bytes = fileBytes(filename)
   %FILEBYTES Read a binary artifact for byte-stability assertions.
   fid = fopen(filename, 'r');
   cleanup = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
end

function cleanupMarkerTree(filename, mar_dir, had_mar_dir, forcing_dir, ...
      had_forcing_dir)
   %CLEANUPMARKERTREE Delete temporary marker files and directories.
   if isfile(filename)
      delete(filename)
   end
   if ~had_mar_dir && isfolder(mar_dir)
      rmdir(mar_dir)
   end
   if ~had_forcing_dir && isfolder(forcing_dir)
      rmdir(forcing_dir)
   end
end
