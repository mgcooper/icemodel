function tests = test_imau_parsers
   %TEST_IMAU_PARSERS Verify IMAU hourly PANGAEA parser contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the verification path and allocate a parser fixture folder.
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

function test_hourly_parser_maps_corrected_native_channels(testCase)
   % Corrected IMAU channels should map to icemodel-native variable names and
   % units while preserving source metadata.
   filename = fullfile(testCase.TestData.tmp, 'VanTiggelen-etal_2024_S21.tab');
   writeImauFixture(filename, "S21", false);

   [returned, meta] = ...
      icemodel.forcing.helpers.readImauHourlyTable(filename);

   expected = ["height_1", "height_2", "tair", "rh", "wspd", "wdir", ...
      "psfc", "swd", "swu", "lwd", "lwu", "surface_elevation_raw", ...
      "surface_elevation", "surface_height_change_sonic", ...
      "surface_height_change", "boom_height", "lat_wgs84", ...
      "lon_wgs84", "albedo", "tsfc"];
   testCase.verifyEqual(string(returned.Properties.VariableNames), expected);
   testCase.verifyEqual(string(returned.Properties.DimensionNames(1)), ...
      "time");
   testCase.verifyEqual(returned.tair(1), 273.15);
   testCase.verifyEqual(returned.tsfc(1), 268.15);
   testCase.verifyEqual(returned.psfc(1), 90000);
   testCase.verifyEqual(returned.rh(1), 80);
   testCase.verifyEqual(returned.wspd(1), 5);
   testCase.verifyEqual(returned.time(1), ...
      datetime(2014, 4, 12, 0, 0, 0, 'TimeZone', 'UTC'));
   testCase.verifyEqual(string(meta.site_id), "S21");
   testCase.verifyEqual(string(meta.doi), "10.1594/PANGAEA.969585");
   testCase.verifyEqual(string(meta.bundle_doi), "10.1594/PANGAEA.971647");
   testCase.verifyEqual(meta.event.lat_wgs84, 66.36);
   testCase.verifyEqual(meta.row_summary.lat_median, 66.181304);
end

function test_hourly_parser_accepts_optional_wire_length(testCase)
   % S22/S23 include a draw-wire column while S21 does not.
   filename = fullfile(testCase.TestData.tmp, 'VanTiggelen-etal_2024_S22.tab');
   writeImauFixture(filename, "S22", true);

   returned = icemodel.forcing.helpers.readImauHourlyTable(filename);

   testCase.verifyTrue(ismember("wire_length", ...
      string(returned.Properties.VariableNames)));
   testCase.verifyEqual(returned.wire_length(1), 6.7);
end

function test_hourly_parser_preserves_blank_cells_as_nan(testCase)
   % Blank numeric PANGAEA cells are missing values, not zeros.
   filename = fullfile(testCase.TestData.tmp, 'VanTiggelen-etal_2024_S23.tab');
   writeImauFixture(filename, "S23", true, blank_corrected=true);

   returned = icemodel.forcing.helpers.readImauHourlyTable(filename);

   testCase.verifyTrue(isnan(returned.tair(1)));
   testCase.verifyTrue(isnan(returned.rh(1)));
   testCase.verifyTrue(isnan(returned.wspd(1)));
   testCase.verifyTrue(isnan(returned.albedo(1)));
end

function test_hourly_parser_normalizes_surface_temperature_fill(testCase)
   % The S22/S23 derived-temperature fill must become missing before the
   % source Celsius values are converted to canonical kelvin.
   filename = fullfile(testCase.TestData.tmp, 'VanTiggelen-etal_2024_S22.tab');
   writeImauFixture(filename, "S22", true, surface_temp="-1273.05");

   returned = icemodel.forcing.helpers.readImauHourlyTable(filename);

   testCase.verifyTrue(isnan(returned.tsfc(1)));
end

function test_hourly_parser_rejects_missing_tabular_section(testCase)
   % A cache artifact without the Date/Time table should fail clearly.
   filename = fullfile(testCase.TestData.tmp, 'bad.tab');
   writeText(filename, sprintf(['/* DATA DESCRIPTION:\n' ...
      'Citation:\tmissing table\n*/\nnot a data table\n']));

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readImauHourlyTable(filename), ...
      'icemodel:forcing:readImauHourlyTable:missingHeader');
end

function writeImauFixture(filename, site_id, include_wire, kwargs)
   %WRITEIMAUFIXTURE Write a tiny PANGAEA-shaped hourly table.
   arguments
      filename (1, 1) string
      site_id (1, 1) string
      include_wire (1, 1) logical
      kwargs.blank_corrected (1, 1) logical = false
      kwargs.surface_temp (1, 1) string = "-5"
   end

   station_doi = stationDoi(site_id);
   metadata = sprintf(['/* DATA DESCRIPTION:\n' ...
      'Citation:\tVan Tiggelen et al. (2024): Hourly data at %s ' ...
      '[dataset]. PANGAEA, https://doi.org/%s,\n' ...
      '\tIn: bundled publication]. PANGAEA, ' ...
      'https://doi.org/10.1594/PANGAEA.971647\n' ...
      'Event(s):\tGRL_%s (%s) * LATITUDE: 66.360000 * ' ...
      'LONGITUDE: -39.310000 * DATE/TIME START: ' ...
      '2014-04-12T00:00:00 * DATE/TIME END: 2014-04-12T01:00:00 ' ...
      '* ELEVATION: 1615.0 m\n' ...
      '*/\n'], site_id, station_doi, site_id, site_id);

   headers = imauHeaders(include_wire);
   row = ["2014-04-12T00:00:00", "2", "10", "-1", "0", ...
      "1.1", "1.2", "85", "80", "4", "5", "180", "900", ...
      "100", "50", "200", "250", "3.4", "3.3", "0.1", ...
      "0.2", "3.5", "66.181304", "-39.042994", "0.7", ...
      kwargs.surface_temp];
   if kwargs.blank_corrected
      row([5 9 11 25]) = "";
   end
   if include_wire
      row(end + 1) = "6.7";
   end

   text = metadata + strjoin(headers, sprintf('\t')) + newline ...
      + strjoin(row, sprintf('\t')) + newline;
   writeText(filename, text);
end

function headers = imauHeaders(include_wire)
   %IMAUHEADERS Return the PANGAEA hourly header columns used by parser tests.
   headers = ["Date/Time", "Height [m] (height 1)", ...
      "Height [m] (height 2)", "TTT [degC]", ...
      "TTT [degC] (corrected at 2m height)", ...
      "Humidity spec [g/kg]", ...
      "Humidity spec [g/kg] (corrected at 2m height)", ...
      "RH [%]", "RH [%] (corrected at 2m height)", "ff [m/s]", ...
      "FF10 [m/s] (corrected at 10m height)", "dd [deg]", ...
      "PPPP [hPa]", "SWD [W/m**2]", "SWU [W/m**2]", ...
      "LWD [W/m**2]", "LWU [W/m**2]", "Surf elev [m] raw", ...
      "Surf elev [m] filtered", "Surf elev change [m] sonic", ...
      "Surf elev change [m] sonic and ablation", "Height [m] boom", ...
      "Latitude", "Longitude", "Alb frac", "Surf temp [degC]"];
   if include_wire
      headers(end + 1) = "Wire length [m]";
   end
end

function doi = stationDoi(site_id)
   %STATIONDOI Return the per-station PANGAEA DOI used in fixture metadata.
   switch string(site_id)
      case "S21"
         doi = "10.1594/PANGAEA.969585";
      case "S22"
         doi = "10.1594/PANGAEA.969629";
      case "S23"
         doi = "10.1594/PANGAEA.969631";
   end
end

function writeText(filename, text)
   %WRITETEXT Write one small parser fixture.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', text);
   clear cleaner
end
