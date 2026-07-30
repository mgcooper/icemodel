function tests = test_modis_staging
   %TEST_MODIS_STAGING Verify MODIS albedo staging and met-cadence attachment.
   %
   % Covers icemodel.forcing.stageModisAlbedo (staging against the repo-local
   % 2012 GEUS fixture; no dependence on the external source volume), the
   % multi-point extension of icemodel.forcing.readGeusModis, and the SSOT
   % attachment helper icemodel.forcing.modisToMetCadence on synthetic daily
   % series (exactness, gap refusal, support mask, bounds fail-closed).
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install project dependencies and allocate one isolated artifact tree.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = string(tempname);
   mkdir(testCase.TestData.tmp)
end

function teardown(testCase)
   % Remove generated met and artifact fixtures after each test.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

%% modisToMetCadence

function test_modisToMetCadence_matches_daily_samples_exactly(testCase)
   % Met samples landing exactly on finite daily samples reproduce the daily
   % values without interpolation error.
   time_daily = dailyAxis(10);
   albedo_daily = linspace(0.3, 0.8, 10)';
   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      albedo_daily, time_daily, time_daily);
   expected = albedo_daily;
   testCase.verifyEqual(returned, expected);
   testCase.verifyTrue(all(support));
end

function test_modisToMetCadence_linear_between_neighbors(testCase)
   % A midday met sample is the exact mean of the bracketing daily samples
   % under linear interpolation.
   time_daily = dailyAxis(3);
   albedo_daily = [0.2; 0.6; 0.4];
   time_met = time_daily(1) + hours(12) + days(0:1)';
   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      albedo_daily, time_daily, time_met);
   expected = [mean(albedo_daily(1:2)); mean(albedo_daily(2:3))];
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(support));
end

function test_modisToMetCadence_refuses_wide_gaps_by_default(testCase)
   % An 8-day hole exceeds the documented 5-day default, so met samples
   % strictly inside the hole are unsupported while the finite edges and the
   % bridged regions stay supported.
   time_daily = dailyAxis(15);
   albedo_daily = 0.5 * ones(15, 1);
   albedo_daily(4:10) = NaN;
   time_met = (time_daily(1):hours(6):time_daily(end))';
   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      albedo_daily, time_daily, time_met);

   % Samples strictly between the finite neighbors of the hole (day 3 and
   % day 11) are refused; samples exactly on those days remain supported.
   inside = time_met > time_daily(3) & time_met < time_daily(11);
   testCase.verifyTrue(all(~support(inside)));
   testCase.verifyTrue(all(isnan(returned(inside))));
   testCase.verifyTrue(all(support(~inside)));
   testCase.verifyEqual(returned(~inside), 0.5 * ones(sum(~inside), 1));
end

function test_modisToMetCadence_max_gap_parameter_bridges_wide_gap(testCase)
   % Raising max_gap above the hole width restores the linear bridge, so the
   % parameter, not a hidden constant, governs refusal.
   time_daily = dailyAxis(15);
   albedo_daily = 0.5 * ones(15, 1);
   albedo_daily(4:10) = NaN;
   time_met = (time_daily(1):hours(6):time_daily(end))';
   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      albedo_daily, time_daily, time_met, max_gap=days(9));
   testCase.verifyTrue(all(support));
   testCase.verifyEqual(returned, 0.5 * ones(numel(time_met), 1), ...
      'AbsTol', 1e-12);
end

function test_modisToMetCadence_no_extrapolation_outside_support(testCase)
   % Met samples before the first and after the last finite daily sample are
   % unsupported: attachment never invents albedo beyond the observed span.
   time_daily = dailyAxis(5);
   albedo_daily = [NaN; 0.4; 0.5; 0.6; NaN];
   time_met = (time_daily(1) - days(1):hours(12):time_daily(end) + days(1))';
   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      albedo_daily, time_daily, time_met);
   outside = time_met < time_daily(2) | time_met > time_daily(4);
   testCase.verifyTrue(all(~support(outside)));
   testCase.verifyTrue(all(isnan(returned(outside))));
   testCase.verifyTrue(all(support(~outside)));
end

function test_modisToMetCadence_single_and_empty_support(testCase)
   % One finite daily sample supports only exact met matches; an all-NaN
   % series supports nothing.
   time_daily = dailyAxis(3);
   albedo_daily = [NaN; 0.4; NaN];
   time_met = [time_daily(2) - hours(1); time_daily(2); ...
      time_daily(2) + hours(1)];
   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      albedo_daily, time_daily, time_met);
   testCase.verifyEqual(returned, [NaN; 0.4; NaN]);
   testCase.verifyEqual(support, [false; true; false]);

   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      nan(3, 1), time_daily, time_met);
   testCase.verifyTrue(all(isnan(returned)));
   testCase.verifyFalse(any(support));
end

function test_modisToMetCadence_rejects_bad_inputs(testCase)
   % Unsorted axes, misaligned sizes, non-positive gaps, and out-of-bounds
   % daily values all fail loudly instead of producing silent output.
   time_daily = dailyAxis(3);
   time_met = dailyAxis(2);
   testCase.verifyError(@() icemodel.forcing.modisToMetCadence( ...
      [0.4; 0.5; 0.6], time_daily([2 1 3]), time_met), ...
      'icemodel:forcing:modisToMetCadence:unsortedDailyAxis');
   testCase.verifyError(@() icemodel.forcing.modisToMetCadence( ...
      [0.4; 0.5], time_daily, time_met), ...
      'icemodel:forcing:modisToMetCadence:sizeMismatch');
   testCase.verifyError(@() icemodel.forcing.modisToMetCadence( ...
      [0.4; 0.5; 0.6], time_daily, time_met, max_gap=days(0)), ...
      'icemodel:forcing:modisToMetCadence:invalidMaxGap');

   % The reconstruction bounds are enforced fail-closed through the canonical
   % interpolation helper, not restated locally.
   bounds = icemodel.forcing.reconstruct.physicalBounds("albedo");
   testCase.verifyError(@() icemodel.forcing.modisToMetCadence( ...
      [0.4; bounds(2) + 0.01; 0.6], time_daily, time_met), ...
      'icemodel:forcing:dailyToHourly:sourceOutOfBounds');
end

function test_modisToMetCadence_aligns_time_zones(testCase)
   % Staged artifacts are UTC while some met axes are unzoned; the met axis
   % owns the zone (loadmet swap parity), so mixed callers still convert.
   time_daily = dailyAxis(3);
   albedo_daily = [0.4; 0.5; 0.6];
   time_met = datetime(2012, 1, 2, 0, 0, 0) + hours(0:12)';
   [returned, support] = icemodel.forcing.modisToMetCadence( ...
      albedo_daily, time_daily, time_met);
   testCase.verifyEqual(returned(1), 0.5);
   testCase.verifyTrue(all(support));
end

%% readGeusModis multi-point extension

function test_readGeusModis_multipoint_matches_single(testCase)
   % A [lat lon] row list returns one column per point, each identical to
   % the corresponding single-point read, plus per-column cell provenance.
   filename = fixtureModisFile();
   points = [67.067 -48.836; 67.096 -49.951];   % KAN_M, KAN_L
   [multi, time_multi, selection] = icemodel.forcing.readGeusModis( ...
      filename, points, "nearest");
   [single_1, time_single] = icemodel.forcing.readGeusModis( ...
      filename, points(1, :), "nearest");
   single_2 = icemodel.forcing.readGeusModis(filename, points(2, :));

   testCase.verifyEqual(size(multi, 2), 2);
   testCase.verifyEqual(multi(:, 1), single_1);
   testCase.verifyEqual(multi(:, 2), single_2);
   testCase.verifyEqual(time_multi, time_single);

   % Nearest-cell selections pin one cell each; the query offset can never
   % exceed the half-diagonal of a 5 km cell and the pinned cell centre must
   % sit beside the query point.
   testCase.verifyEqual(numel(selection), 2);
   for q = 1:2
      testCase.verifyEqual(selection(q).count, [1 1]);
      testCase.verifyEqual(strlength(string(selection(q).grid_sha256)), 64);
      testCase.verifyLessThan(selection(q).distance_m, 5000 / sqrt(2));
      testCase.verifyEqual(selection(q).cell_lat, points(q, 1), ...
         'AbsTol', 0.1);
      testCase.verifyEqual(selection(q).cell_lon, points(q, 2), ...
         'AbsTol', 0.25);
   end
end

%% stageModisAlbedo

function test_stageModisAlbedo_stages_fixture_artifact(testCase)
   % Staging one site from the 2012 fixture writes a window-stamped daily
   % artifact through the canonical writer, with values matching a direct
   % bounds-masked reader extraction and full byte-pinned provenance.
   met_dir = fullfile(testCase.TestData.tmp, "met");
   outdir = fullfile(testCase.TestData.tmp, "userdata");
   latlon = [67.067 -48.836];
   met_file = writeSyntheticMetArtifact(met_dir, "kanm", "promice", latlon);

   info = icemodel.forcing.stageModisAlbedo("kanm", ...
      modis_dir=fixtureModisDir(), met_dir=met_dir, outdir=outdir, ...
      years=2012);

   % Canonical window-stamped native-daily name inside the per-source folder.
   expected = string(fullfile(outdir, "modis", ...
      "kanm_modis_20120101_20121231_86400s.mat"));
   testCase.verifyEqual(info.filename, expected);
   testCase.verifyTrue(isfile(expected));

   % The artifact holds a daily timetable named Data spanning all of 2012
   % with the daily cadence stamped in the top-level metadata.
   loaded = load(expected);
   testCase.verifyTrue(istimetable(loaded.Data));
   testCase.verifyEqual(height(loaded.Data), 366);
   testCase.verifyEqual( ...
      loaded.artifact_metadata.artifact_cadence_seconds, 86400);
   testCase.verifyEqual( ...
      string(loaded.Data.Properties.VariableNames), "albedo");

   % Values equal a direct reader extraction with the reconstruction bounds
   % applied; axis days beyond the truncated fixture stay missing.
   [expected_albedo, time_daily] = icemodel.forcing.readGeusModis( ...
      fixtureModisFile(), latlon, "nearest");
   bounds = icemodel.forcing.reconstruct.physicalBounds("albedo");
   expected_albedo(expected_albedo < bounds(1) ...
      | expected_albedo > bounds(2)) = NaN;
   returned = loaded.Data.albedo(1:numel(time_daily));
   testCase.verifyEqual(returned, expected_albedo);
   testCase.verifyTrue(all(isnan(loaded.Data.albedo( ...
      numel(time_daily) + 1:end))));

   % Extraction and byte-pinning provenance ride in the artifact metadata.
   metadata = loaded.Data.Properties.UserData;
   testCase.verifyEqual(string(metadata.extraction_method), "nearest");
   testCase.verifyEqual(metadata.extraction_cell_count, [1 1]);
   testCase.verifyLessThan(metadata.extraction_cell_distance_m, ...
      5000 / sqrt(2));
   testCase.verifyEqual(metadata.albedo_physical_bounds, bounds);
   source_info = dir(fixtureModisFile());
   testCase.verifyEqual(metadata.source_files(1).size_bytes, ...
      source_info.bytes);
   testCase.verifyEqual(string(metadata.source_files(1).sha256), ...
      icemodel.verification.setup.fileSha256(fixtureModisFile()));
   testCase.verifyEqual(metadata.source_files(1).year, 2012);
   testCase.verifyEqual(metadata.source_grid_size, selectionGridSize( ...
      fixtureModisFile()));
   testCase.verifyEqual(strlength(string(metadata.source_grid_sha256)), 64);
   testCase.verifyEqual(string(metadata.attachment_helper), ...
      "icemodel.forcing.modisToMetCadence");
   testCase.verifyEqual(string(metadata.met_identity_file), ...
      string(met_file));

   % The location contract is copied verbatim from the met artifact.
   custom = loaded.Data.Properties.CustomProperties;
   testCase.verifyEqual(custom.Lat, latlon(1));
   testCase.verifyEqual(custom.Lon, latlon(2));
   testCase.verifyEqual(numel(custom.ScalarUnits), 6);
end

function test_stageModisAlbedo_clips_fixed_366_day_records(testCase)
   % The C6 product stores a fixed 366-day record per file, so a non-leap
   % file carries one trailing day of the next year: the next file's native
   % day must win the interior overlap and days beyond the staged axis must
   % be dropped, keeping the artifact axis strictly one row per date.
   modis_dir = fullfile(testCase.TestData.tmp, "gris");
   met_dir = fullfile(testCase.TestData.tmp, "met");
   outdir = fullfile(testCase.TestData.tmp, "userdata");

   % 2013 ends with a distinctive overhang value that 2014's native first
   % day must replace; 2014's own overhang value must be dropped entirely.
   values_2013 = 0.5 * ones(1, 366);
   values_2013(366) = 0.9;
   values_2014 = 0.4 * ones(1, 366);
   values_2014(1) = 0.2;
   values_2014(366) = 0.7;
   [~, latlon] = writeTinyReflectivityFile(modis_dir, 2013, values_2013);
   writeTinyReflectivityFile(modis_dir, 2014, values_2014);
   writeSyntheticMetArtifact(met_dir, "kanm", "promice", latlon);

   info = icemodel.forcing.stageModisAlbedo("kanm", ...
      modis_dir=modis_dir, met_dir=met_dir, outdir=outdir);

   % One row per date across 2013-2014 (730 days, no leap year involved).
   loaded = load(info.filename);
   testCase.verifyEqual(height(loaded.Data), 730);
   returned = loaded.Data.albedo;
   testCase.verifyEqual(returned(365), 0.5);   % native 2013 Dec 31
   testCase.verifyEqual(returned(366), 0.2);   % native 2014 Jan 1 wins
   testCase.verifyEqual(returned(367), 0.4);   % native 2014 Jan 2
   testCase.verifyEqual(returned(730), 0.4);   % native 2014 Dec 31
   testCase.verifyFalse(any(returned == 0.9)); % 2013 overhang replaced
   testCase.verifyFalse(any(returned == 0.7)); % 2014 overhang dropped
end

function test_stageModisAlbedo_rejects_shifted_yearly_grid(testCase)
   % Matching selected-cell metadata is insufficient identity: an interior
   % coordinate can move while the selected cell and reconstructed endpoints
   % stay fixed, so the complete raw coordinate grid must be fingerprinted.
   modis_dir = fullfile(testCase.TestData.tmp, "gris");
   met_dir = fullfile(testCase.TestData.tmp, "met");
   outdir = fullfile(testCase.TestData.tmp, "userdata");
   values = 0.5 * ones(1, 366);
   [~, latlon] = writeTinyReflectivityFile(modis_dir, 2013, values);
   writeTinyReflectivityFile(modis_dir, 2014, values, 100);
   writeSyntheticMetArtifact(met_dir, "kanm", "promice", latlon);

   testCase.verifyError(@() icemodel.forcing.stageModisAlbedo("kanm", ...
      modis_dir=modis_dir, met_dir=met_dir, outdir=outdir), ...
      'icemodel:forcing:stageModisAlbedo:gridMismatch');
end

function test_stageModisAlbedo_rejects_unsafe_site_before_paths(testCase)
   % Station tokens reach met globs and output names, so traversal text must
   % fail before an unavailable source directory can be inspected.
   testCase.verifyError(@() icemodel.forcing.stageModisAlbedo("../x", ...
      modis_dir=fullfile(testCase.TestData.tmp, "missing")), ...
      'icemodel:reconstruct:mustBeStationToken:invalidToken');
end

function test_stageModisAlbedo_resolves_metadata_location(testCase)
   % Gap-filled PROMICE met artifacts save their location as top-level
   % artifact_metadata lat/lon/elev instead of CustomProperties; staging
   % resolves that route and projects the missing EPSG:3413 identity.
   met_dir = fullfile(testCase.TestData.tmp, "met");
   outdir = fullfile(testCase.TestData.tmp, "userdata");
   latlon = [67.067 -48.836];
   writeSyntheticMetArtifact(met_dir, "kanm", "promice", latlon, "metadata");

   info = icemodel.forcing.stageModisAlbedo("kanm", ...
      modis_dir=fixtureModisDir(), met_dir=met_dir, outdir=outdir, ...
      years=2012);

   loaded = load(info.filename);
   custom = loaded.Data.Properties.CustomProperties;
   testCase.verifyEqual(custom.Lat, latlon(1));
   testCase.verifyEqual(custom.Lon, latlon(2));
   testCase.verifyEqual(custom.Elev, 1270);
   testCase.verifyTrue(isfinite(custom.X) && isfinite(custom.Y));
end

function test_stageModisAlbedo_bare_met_location_errors(testCase)
   % A met artifact with neither metadata nor CustomProperties location
   % cannot supply a colocation identity and must fail loudly.
   met_dir = fullfile(testCase.TestData.tmp, "met");
   writeSyntheticMetArtifact(met_dir, "kanm", "promice", ...
      [67.067 -48.836], "bare");
   testCase.verifyError(@() icemodel.forcing.stageModisAlbedo("kanm", ...
      modis_dir=fixtureModisDir(), met_dir=met_dir, ...
      outdir=fullfile(testCase.TestData.tmp, "userdata"), years=2012), ...
      'icemodel:forcing:stageModisAlbedo:badMetLocation');
end

function test_stageModisAlbedo_missing_met_errors(testCase)
   % A site with no staged met artifact must fail before any writing: the
   % MODIS artifact cannot claim a colocation identity it does not have.
   met_dir = fullfile(testCase.TestData.tmp, "met");
   mkdir(met_dir)
   testCase.verifyError(@() icemodel.forcing.stageModisAlbedo("kanm", ...
      modis_dir=fixtureModisDir(), met_dir=met_dir, ...
      outdir=fullfile(testCase.TestData.tmp, "userdata"), years=2012), ...
      'icemodel:forcing:stageModisAlbedo:metArtifactNotFound');
end

function test_stageModisAlbedo_missing_year_errors(testCase)
   % Explicitly requested years without a source file fail loudly so a mount
   % hiccup cannot silently stage a partial window.
   testCase.verifyError(@() icemodel.forcing.stageModisAlbedo("kanm", ...
      modis_dir=fixtureModisDir(), ...
      met_dir=fullfile(testCase.TestData.tmp, "met"), ...
      outdir=fullfile(testCase.TestData.tmp, "userdata"), years=1999), ...
      'icemodel:forcing:stageModisAlbedo:sourceYearMissing');
end

%% Local fixtures

function time_daily = dailyAxis(ndays)
   %DAILYAXIS UTC daily datetimes matching the staged-artifact axis contract.
   time_daily = datetime(2012, 1, 1, 'TimeZone', 'UTC') + days(0:ndays - 1)';
end

function gris = fixtureModisDir()
   %FIXTUREMODISDIR Repo-local GEUS 2012 fixture directory (no volume dep).
   gris = string(fullfile(icemodel.internal.fullpath('test'), ...
      'data', 'forcing', 'geus', 'albedo', 'gris'));
end

function filename = fixtureModisFile()
   %FIXTUREMODISFILE The single 2012 reflectivity fixture NetCDF.
   filename = string(fullfile(fixtureModisDir(), ...
      'Greenland_Reflectivity_2012_5km_C6.nc'));
end

function [filename, latlon] = writeTinyReflectivityFile(modis_dir, ...
      yyyy, values, shift_m)
   %WRITETINYREFLECTIVITYFILE Attribute-free two-cell GEUS-style source file
   % holding one uniform daily albedo record, so overlap/clip behavior can
   % be asserted against exactly known values. LATLON is the first cell
   % centre, a valid nearest-cell staging target.
   if nargin < 4
      shift_m = 0;
   end
   if ~isfolder(modis_dir)
      mkdir(modis_dir)
   end
   filename = fullfile(modis_dir, sprintf( ...
      'Greenland_Reflectivity_%d_5km_C6.nc', yyyy));
   [X, Y] = ndgrid([0; 5000; 10000], [0, 5000, 10000]);
   X(2, 2) = X(2, 2) + shift_m;
   [lat, lon] = projinv( ...
      icemodel.forcing.helpers.geusModisProjection(), X, Y);
   latlon = [lat(1, 1), lon(1, 1)];
   nccreate(filename, 'lat', 'Dimensions', {'x', 3, 'y', 3});
   nccreate(filename, 'lon', 'Dimensions', {'x', 3, 'y', 3});
   nccreate(filename, 'albedo', ...
      'Dimensions', {'x', 3, 'y', 3, 'time', numel(values)});
   ncwrite(filename, 'lat', lat);
   ncwrite(filename, 'lon', lon);
   ncwrite(filename, 'albedo', repmat( ...
      reshape(values, 1, 1, []), 3, 3, 1));
end

function grid_size = selectionGridSize(filename)
   %SELECTIONGRIDSIZE Return the source grid dimensions asserted in metadata.
   grid_size = size(ncread(filename, 'lat'));
end

function met_file = writeSyntheticMetArtifact(met_dir, site, source, ...
      latlon, form)
   %WRITESYNTHETICMETARTIFACT Minimal met artifact carrying a location
   % identity; only that identity, not the met channels, is consumed by the
   % staging lookup. FORM selects the identity route under test:
   % "customprops" attaches the full location CustomProperties contract,
   % "metadata" saves the gap-filled PROMICE style top-level
   % artifact_metadata (lat/lon/elev) beside a bare met timetable, and
   % "bare" writes neither so the identity lookup must fail.
   arguments
      met_dir (1, 1) string
      site (1, 1) string
      source (1, 1) string
      latlon (1, 2) double
      form (1, 1) string {mustBeMember(form, ...
         ["customprops", "metadata", "bare"])} = "customprops"
   end
   mkdir(met_dir)
   Time = datetime(2012, 1, 1, 'TimeZone', 'UTC') + minutes(0:15:45)';
   met = timetable(Time, zeros(4, 1), 'VariableNames', {'tair'});
   met_file = fullfile(met_dir, sprintf( ...
      'met_%s_%s_20120101_20121231_15m.mat', site, source));
   switch form
      case "customprops"
         location = struct('lat_wgs84', latlon(1), 'lon_wgs84', latlon(2), ...
            'elev_m', 1270, 'slope', NaN);
         met = icemodel.forcing.helpers.attachLocationMetadata( ...
            met, location);
         save(met_file, 'met')
      case "metadata"
         artifact_metadata = struct('site', upper(site), ...
            'lat', latlon(1), 'lon', latlon(2), 'elev', 1270);
         save(met_file, 'met', 'artifact_metadata')
      case "bare"
         save(met_file, 'met')
   end
end
