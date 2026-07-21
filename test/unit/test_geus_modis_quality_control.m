function tests = test_geus_modis_quality_control
   %TEST_GEUS_MODIS_QUALITY_CONTROL Verify source masking and bounded repair.
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
   % Remove generated source and artifact fixtures after each test.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

function test_normalizer_masks_only_nonphysical_values(testCase)
   % GEUS 999 and other nonphysical/nonfinite samples become missing while the
   % complete closed physical albedo interval remains byte-for-byte numeric.
   raw = [NaN, -Inf, -0.01, 0, 0.5, 1, 1.01, 999, Inf];
   actual = icemodel.forcing.helpers.normalizeGeusModisAlbedo(raw);
   expected = [NaN, NaN, NaN, 0, 0.5, 1, NaN, NaN, NaN];

   testCase.verifyEqual(actual, expected);
   testCase.verifyClass(actual, 'double');
   testCase.verifyTrue(isnan( ...
      icemodel.forcing.helpers.normalizeGeusModisAlbedo(0.5 + 1i)));
end

function test_reader_masks_undocumented_999_sentinel(testCase)
   % A source file with no fill/range attributes must still expose its finite 999
   % sentinel as NaN before point or polygon collapse. The misleading year in
   % the parent path must not replace the product year encoded in the basename.
   modis_dir = fullfile(testCase.TestData.tmp, ...
      "scratch_2825_parent", "reader-source");
   [filename, location] = writeTinyGeusModis( ...
      modis_dir, 2012, [0.4, 999, 0.8]);
   info = ncinfo(filename, 'albedo');

   [albedo, Time] = icemodel.forcing.readGeusModis( ...
      filename, location, "nearest");

   testCase.verifyEmpty(info.Attributes);
   testCase.verifyEqual(albedo, [0.4; NaN; 0.8]);
   expected_time = (datetime(2012, 1, 1, TimeZone="UTC") + days(0:2))';
   testCase.verifyEqual(posixtime(Time), posixtime(expected_time));
   testCase.verifyEqual(string(Time.TimeZone), "UTC");
end

function test_reader_renormalizes_partial_polygon_coverage(testCase)
   % A sentinel in one selected cell must not contaminate the valid-cell mean;
   % collapsing validity weights separately yields the mean of physical samples.
   modis_dir = fullfile(testCase.TestData.tmp, "polygon-source");
   [filename, ~] = writeTinyGeusModis(modis_dir, 2012, 0);
   values = reshape([0.2, 0.4, 0.6, 999], 2, 2, 1);
   ncwrite(filename, 'albedo', values);
   lat = ncread(filename, 'lat');
   lon = ncread(filename, 'lon');
   [x, y] = projfwd(projcrs(3413), lat, lon);
   margin = 1000;
   polygon = polyshape( ...
      [min(x(:)) - margin, max(x(:)) + margin, ...
      max(x(:)) + margin, min(x(:)) - margin], ...
      [min(y(:)) - margin, min(y(:)) - margin, ...
      max(y(:)) + margin, max(y(:)) + margin]);

   albedo = icemodel.forcing.readGeusModis( ...
      filename, polygon, remap="equal");

   testCase.verifyEqual(albedo, 0.4, AbsTol=1e-12);
end

function test_channel_preserves_missing_year_and_rejects_duplicates(testCase)
   % Missing optional years remain NaN, while duplicate files for one requested
   % year are rejected instead of choosing an arbitrary source. The second
   % output classifies the exact full, partial, and absent coverage states.
   modis_dir = fullfile(testCase.TestData.tmp, "channel-source");
   [filename, location] = writeTinyGeusModis( ...
      modis_dir, 2012, [0.25, 0.75]);
   Time = [datetime(2012, 1, 1, TimeZone="UTC"); ...
      datetime(2012, 1, 2, TimeZone="UTC"); ...
      datetime(2013, 1, 1, TimeZone="UTC")];

   [modis, partial] = icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, [2012, 2013], location, "nearest", ...
      "conservative", Time);
   testCase.verifyEqual(modis, [0.25; 0.75; NaN]);
   testCase.verifyEqual(string(partial.modis_product), ...
      "GEUS Greenland Reflectivity 5km C6");
   testCase.verifyEqual(string(partial.modis_status), ...
      "partial_source_coverage");
   testCase.verifyEqual(double(partial.modis_coverage_years(:))', 2012);

   [~, complete] = icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, 2012, location, "nearest", "conservative", Time(1:2));
   testCase.verifyEqual(string(complete.modis_status), "source_coverage");

   % An existing cache with no matching source year is explicit absence, not a
   % fabricated all-missing observation claimed as coverage.
   [missing, absent] = icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, 2013, location, "nearest", "conservative", Time(3));
   testCase.verifyTrue(all(isnan(missing)));
   testCase.verifyEqual(string(absent.modis_status), "no_source_coverage");
   testCase.verifyEmpty(absent.modis_coverage_years);

   % A second matching name makes the year ambiguous even when bytes agree.
   copyfile(filename, fullfile(modis_dir, "duplicate_2012_copy.nc"));
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, 2012, location, "nearest", "conservative", Time(1:2)), ...
       'icemodel:forcing:modisAlbedoChannel:ambiguousFile');
end

function test_coverage_contract_normalizes_and_rejects_false_coverage(testCase)
   % Shape/duplicates cannot change canonical bytes, coverage outside the target
   % axis fails, and a present file with no physical target values is not covered.
   metadata = icemodel.forcing.helpers.geusModisCoverageMetadata( ...
      [2013; 2012; 2013], [2013; 2012; 2012]);
   testCase.verifyEqual(string(metadata.modis_status), "source_coverage");
   testCase.verifyEqual(double(metadata.modis_coverage_years(:))', [2012, 2013]);
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.geusModisCoverageMetadata(2012, 2013), ...
      'icemodel:forcing:geusModisCoverageMetadata:outsideRequest');

   modis_dir = fullfile(testCase.TestData.tmp, "invalid-source");
   [~, location] = writeTinyGeusModis(modis_dir, 2014, [999, 999]);
   Time = datetime(2014, 1, 1:2, TimeZone="UTC")';
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.modisAlbedoChannel( ...
      modis_dir, 2014, location, "nearest", "conservative", Time), ...
      'icemodel:forcing:modisAlbedoChannel:noPhysicalValues');
end

function [filename, location] = writeTinyGeusModis(modis_dir, yyyy, values)
   %WRITETINYGEUSMODIS Create an attribute-free two-cell GEUS-style source.
   if ~isfolder(modis_dir)
      mkdir(modis_dir)
   end
   filename = fullfile(modis_dir, sprintf( ...
      'Greenland_Reflectivity_%d_5km_C6.nc', yyyy));
   [X, Y] = ndgrid([0; 5000], [0, 5000]);
   [lat, lon] = projinv( ...
      icemodel.forcing.helpers.geusModisProjection(), X, Y);
   location = [lat(1, 1), lon(1, 1)];
   nccreate(filename, 'lat', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'lon', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'albedo', ...
      'Dimensions', {'x', 2, 'y', 2, 'time', numel(values)});
   ncwrite(filename, 'lat', lat);
   ncwrite(filename, 'lon', lon);
   ncwrite(filename, 'albedo', repmat(reshape(values, 1, 1, []), 2, 2, 1));
end
