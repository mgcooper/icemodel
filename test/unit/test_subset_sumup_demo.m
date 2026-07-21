function tests = test_subset_sumup_demo
   %TEST_SUBSET_SUMUP_DEMO Verify deterministic committed-fixture derivation.
   tests = functiontests(localfunctions);
end

function test_selects_nearest_profile_and_caps_each_channel(testCase)
   % Each supported channel keeps one nearest profile in its native row order.
   point = [67, -48];
   keys = ["near"; "near"; "near"; "far"; "far"];
   names = ["Near"; "Near"; "Near"; "Far"; "Far"];
   lat = [67; 67.01; 66.99; 75; 75.01];
   lon = [-48; -48.01; -47.99; -40; -40.01];
   dates = datetime(2000, 1, 1, TimeZone="UTC") + days((0:4)');

   density = table(keys, names, lat, lon, dates, (1:5)', ...
      VariableNames={'name_key', 'name', 'latitude', 'longitude', ...
      'datetime', 'density'});
   subsurface_temperature = timetable(dates, keys, names, lat, lon, ...
      (11:15)', VariableNames={'name_key', 'name', 'latitude', ...
      'longitude', 'temperature'});
   smb = table(keys, names, lat, lon, dates, dates + days(1), (21:25)', ...
      VariableNames={'name_key', 'name', 'latitude', 'longitude', ...
      'start_date', 'end_date', 'smb'});
   observations = struct('density', density, ...
      'subsurface_temperature', subsurface_temperature, 'smb', smb);

   [actual, provenance] = ...
      icemodel.verification.setup.subsetSumupDemo( ...
      observations, point, max_rows=2);

   testCase.verifyEqual(actual.density.name_key, ["near"; "near"]);
   testCase.verifyEqual(actual.subsurface_temperature.name_key, ...
      ["near"; "near"]);
   testCase.verifyEqual(actual.smb.name_key, ["near"; "near"]);
   testCase.verifyEqual(actual.density.density, [1; 2]);
   testCase.verifyEqual(actual.subsurface_temperature.temperature, [11; 12]);
   testCase.verifyEqual(actual.smb.smb, [21; 22]);
   testCase.verifyEqual(string({provenance.channels.channel}), ...
      ["density", "subsurface_temperature", "smb"]);
   testCase.verifyEqual([provenance.channels.n_rows], [2, 2, 2]);
   testCase.verifyEqual(string({provenance.channels.name_key}), ...
      ["near", "near", "near"]);
   testCase.verifyEqual(provenance.max_rows, 2);
end

function test_skips_absent_or_empty_channels(testCase)
   % An intentionally narrow fixture may contain only one available channel.
   point = [67, -48];
   density = table("near", 67, -48, 300, ...
      VariableNames={'name_key', 'latitude', 'longitude', 'density'});
   observations = struct('density', density, 'smb', table());

   [actual, provenance] = ...
      icemodel.verification.setup.subsetSumupDemo(observations, point);

   testCase.verifyEqual(height(actual.density), 1);
   testCase.verifyEqual(height(actual.smb), 0);
   testCase.verifyEqual(numel(provenance.channels), 1);
   testCase.verifyEqual(string(provenance.channels.channel), "density");
   testCase.verifyEqual(string(provenance.channels.name), "");
end

function test_requires_integer_row_cap(testCase)
   % Fractional caps cannot define deterministic table indexing.
   observations = struct();
   testCase.verifyError(@() ...
      icemodel.verification.setup.subsetSumupDemo( ...
      observations, [67, -48], max_rows=1.5), ...
      'MATLAB:validators:mustBeInteger');
end
