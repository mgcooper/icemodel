function tests = test_retmip_parsers
   %TEST_RETMIP_PARSERS Verify RetMIP protocol table and output inventory helpers.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the canonical test path and allocate a temporary source folder.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = tempname;
   mkdir(testCase.TestData.tmp);
end

function teardown(testCase)
   % Remove temporary parser fixtures.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

function test_protocol_table_parses_three_hour_cadence(testCase)
   % RetMIP surface protocol tables should parse as 3-hourly canonical series.
   filename = fullfile(testCase.TestData.tmp, 'surface.tab');
   writeText(filename, sprintf([ ...
      'time\tTsurf_K\tmelt_mmweq\tacc_subl_mmweq\n' ...
      '2012-01-01 00:00:00\t260\t0\t1\n' ...
      '2012-01-01 03:00:00\t261\t0.1\t2\n']));

   [returned, meta] = ...
      icemodel.verification.setup.readRetmipProtocolTable(filename);

   testCase.verifyTrue(isdatetime(returned.time));
   expected = 3;
   testCase.verifyEqual(meta.timestep_hours, expected);
   testCase.verifyTrue(ismember("Tsurf_K", meta.raw_variables));
   testCase.verifyTrue(ismember("tsfc", meta.variables));
   testCase.verifyTrue(ismember("melt", meta.variables));
   testCase.verifyTrue(ismember("snowf_subl", meta.variables));
   testCase.verifyTrue(all(ismember(["ppt", "rainf", "snowf"], ...
      meta.variables)));
   testCase.verifyEqual(returned.tsfc(1), 260);
   testCase.verifyEqual(returned.melt(2), 0.1e-3 / expected, ...
      'AbsTol', 1e-15);
   testCase.verifyEqual(returned.snowf_subl(2), 2e-3 / expected, ...
      'AbsTol', 1e-15);
   testCase.verifyEqual(meta.mass_flux_conversion_factor, 1e-3 / expected, ...
      'AbsTol', 1e-15);
   testCase.verifyTrue(contains(meta.mass_flux_policy, "mWE/h"));
   testCase.verifyTrue(isnan(returned.rainf(1)));
end

function test_protocol_table_converts_surface_temperature_celsius(testCase)
   % Celsius protocol aliases should stage as Kelvin tsfc.
   filename = fullfile(testCase.TestData.tmp, 'surface_c.tab');
   writeText(filename, sprintf([ ...
      'time\tTsurf_degC\tnet_acc_mmweq\n' ...
      '2012-01-01 00:00:00\t-13.15\t1\n' ...
      '2012-01-01 03:00:00\t-12.15\t2\n']));

   returned = icemodel.verification.setup.readRetmipProtocolTable(filename);

   testCase.verifyEqual(returned.tsfc(1), 260, 'AbsTol', 1e-12);
   testCase.verifyEqual(returned.snowf_subl(2), 2e-3 / 3, ...
      'AbsTol', 1e-15);
end

function test_protocol_table_rejects_bad_cadence(testCase)
   % A non-3-hour RetMIP protocol series should fail before staging.
   filename = fullfile(testCase.TestData.tmp, 'bad_surface.tab');
   writeText(filename, sprintf([ ...
      'time\tTsurf\n' ...
      '2012-01-01 00:00:00\t-12\n' ...
      '2012-01-01 02:00:00\t-11\n']));

   testCase.verifyError(@() ...
      icemodel.verification.setup.readRetmipProtocolTable(filename), ...
      'icemodel:verification:readRetmipProtocolTable:badCadence');
end

function test_profile_table_requires_depth_and_state(testCase)
   % Initial profile tables must include depth plus a firn-state variable.
   filename = fullfile(testCase.TestData.tmp, 'profile.tab');
   writeText(filename, sprintf([ ...
      'depth\trho\ttemp\tlwc\n' ...
      '0\t350\t263\t0\n']));

   [returned, meta] = ...
      icemodel.verification.setup.readRetmipProfileTable(filename);

   expected = 1;
   testCase.verifyEqual(height(returned), expected);
   testCase.verifyTrue(ismember("depth", meta.variables));
   testCase.verifyTrue(ismember("density", meta.variables));
   testCase.verifyTrue(ismember("subsurface_temperature", meta.variables));
   testCase.verifyTrue(ismember("lwc", meta.variables));
end

function test_output_inventory_reads_netcdf_header(testCase)
   % RetMIP model outputs are indexed by NetCDF header without loading arrays.
   filename = fullfile(testCase.TestData.tmp, 'model_output.nc');
   nccreate(filename, 'temp', 'Dimensions', {'time', 2});
   nccreate(filename, 'rho', 'Dimensions', {'time', 2});
   nccreate(filename, 'lwc', 'Dimensions', {'time', 2});

   returned = icemodel.verification.setup.retmipOutputInventory(filename);

   testCase.verifyTrue(ismember("temp", returned.variables));
   testCase.verifyTrue(ismember("rho", returned.variables));
   testCase.verifyTrue(ismember("lwc", returned.variables));
end

function writeText(filename, text)
   %WRITETEXT Write one small text fixture.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', text);
   clear cleaner
end
