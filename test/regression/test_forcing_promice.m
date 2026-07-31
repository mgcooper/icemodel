function tests = test_forcing_promice
   %TEST_FORCING_PROMICE Verify the PROMICE L3 forcing/evaluation builders.
   %
   % Reads the staged pypromice Level-3 product from the repo-root
   % data/verification/promice tree (hour/<STATION>_hour.nc) and skips
   % cleanly when it is not available.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve the staged L3 product; skip the whole file when absent. Require
   % readable station NetCDF files, not just the folder.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;

   source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice'));
   hasdata = ~isempty(dir(fullfile(source_dir, "hour", "*_hour.nc")));
   testCase.assumeTrue(hasdata, ...
      'PROMICE L3 product not available (data/verification/promice missing)');
   testCase.TestData.source_dir = source_dir;
end

function test_readPromiceAws_standard_channels(testCase)
   % The reader returns icemodel-named channels on a UTC hourly axis with
   % station metadata, with the L3 unit conversions applied.

   [aws, metadata] = icemodel.forcing.readPromiceAws("KAN_L", ...
      source_dir=testCase.TestData.source_dir);

   needed = ["tair", "swd", "lwd", "albedo", "wspd", "rh", "psfc", ...
      "tsfc", "snow_height", "z_ice_surf", "tice1"];
   testCase.verifyTrue(all(ismember(needed, ...
      string(aws.Properties.VariableNames))));
   testCase.verifyEqual(aws.Time.TimeZone, 'UTC');
   testCase.verifyEqual(minute(aws.Time(1)), 0);   % snapped to hour start

   % Unit conversions: tair/tsfc kelvin, psfc pascal.
   testCase.verifyGreaterThan(median(aws.tair, 'omitnan'), 230);   % kelvin
   testCase.verifyLessThan(median(aws.tair, 'omitnan'), 290);
   testCase.verifyGreaterThan(median(aws.psfc, 'omitnan'), 5e4);   % pascal
   testCase.verifyGreaterThan(median(aws.rh, 'omitnan'), 50);      % percent

   testCase.verifyEqual(metadata.lat, 67.097, 'AbsTol', 1e-2);
   testCase.verifyEqual(metadata.n_booms, 1);   % KAN_L is one-boom
   testCase.verifyTrue(metadata.swd_source_file_present);
   testCase.verifyTrue(metadata.swd_corrected_source_file_present);
   testCase.verifyTrue(metadata.swd_source_file_observations_present);
   testCase.verifyTrue(metadata.swu_source_file_observations_present);
end

function test_readPromiceAws_resolves_canonical_station_id(testCase)
   % "KAN_M" and the compact alias "kanm" resolve to the same file.

   [~, m1] = icemodel.forcing.readPromiceAws("KAN_M", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 6, 1), enddate=datetime(2015, 6, 2));
   [~, m2] = icemodel.forcing.readPromiceAws("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 6, 1), enddate=datetime(2015, 6, 2));
   testCase.verifyEqual(m1.source_file, m2.source_file);
end

function test_cloud_fraction_scale_and_corrected_shortwave_persist(testCase)
   % The reader preserves native fractional cloud cover, accepts legacy
   % percent-form input once, and makes invalid values missing. The public
   % native Data table and its saved userdata retain corrected-first SWD/SWU.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   filename = writeSyntheticPromice(root, [-0.1; 0; 0.25; 1]);

   aws = icemodel.forcing.readPromiceAws("KAN_M", source_dir=root);
   testCase.verifyEqual(aws.cfrac, [NaN; 0; 0.25; 1]);

   % A legacy percent-form series is scaled as one product; its sub-one-percent
   % samples must divide too, while values above 100 are explicit missing data.
   ncwrite(filename, 'cc', [0.5; 25; 100; 150]);
   aws = icemodel.forcing.readPromiceAws("KAN_M", source_dir=root);
   testCase.verifyEqual(aws.cfrac, [0.005; 0.25; 1; NaN]);

   Data = icemodel.forcing.buildPromiceData("KAN_M", source_dir=root);
   testCase.verifyEqual(Data.cfrac, [0.005; 0.25; 1; NaN]);
   testCase.verifyEqual(Data.swd, [11; 20; 33; 0]);
   testCase.verifyEqual(Data.swu, [1.1; 2; 3.3; 0], 'AbsTol', 1e-12);

   % writeuserdata is the native runtime persistence boundary used by the
   % importer, so reload the saved table rather than trusting its input.
   outdir = fullfile(root, "userdata");
   filenames = icemodel.forcing.helpers.writeuserdata(Data, ...
      "kanm", "promice", outdir=outdir);
   saved = load(filenames(1), 'Data');
   for name = ["cfrac", "swd", "swu"]
      testCase.verifyEqual(saved.Data.(name), Data.(name));
   end
end

function test_promice_builders_enforce_paired_utc_window(testCase)
   % Both public products share one paired-window boundary through the reader.
   source_dir = testCase.TestData.source_dir;
   error_id = 'icemodel:internal:pairedWindow:invalidWindow';
   missing_source = fullfile(source_dir, 'missing-window-precedence');
   testCase.verifyError(@() icemodel.forcing.buildPromiceMet("KAN_M", ...
      source_dir=missing_source, startdate="2015-06-01"), error_id);
   testCase.verifyError(@() icemodel.forcing.buildPromiceData("KAN_M", ...
      source_dir=missing_source, enddate="2015-06-02"), error_id);
   testCase.verifyFalse(isfolder(missing_source));

   builders = { ...
      @(a, b) icemodel.forcing.buildPromiceMet("KAN_M", ...
      source_dir=source_dir, startdate=a, enddate=b)
      @(a, b) icemodel.forcing.buildPromiceData("KAN_M", ...
      source_dir=source_dir, startdate=a, enddate=b)};
   for k = 1:numel(builders)
      build = builders{k};
      testCase.verifyError(@() build("2015-06-01", ""), error_id);
      testCase.verifyError(@() build("", "2015-06-02"), error_id);
      testCase.verifyError(@() build("2015-06-02", "2015-06-01"), error_id);
   end

   % A non-UTC representation must retain its instant and return a UTC axis.
   expected_start = datetime(2015, 6, 1, 0, 0, 0, 'TimeZone', 'UTC');
   expected_end = expected_start + hours(2);
   window_start = expected_start;
   window_end = expected_end;
   window_start.TimeZone = 'America/New_York';
   window_end.TimeZone = 'America/New_York';
   for k = 1:numel(builders)
      built = builders{k}(window_start, window_end);
      testCase.verifyEqual(built.Time.TimeZone, 'UTC');
      testCase.verifyEqual(built.Time(1), expected_start);
      testCase.verifyEqual(built.Time(end), expected_end);
   end
end

function test_negative_derived_net_shortwave_is_missing(testCase)
   % Physically inconsistent swu>swd samples remain visible in their source
   % components but cannot enter evaluation through a negative net shortwave.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   filename = writeSyntheticPromice(root, [0; 0; 0; 0]);
   ncwrite(filename, 'usr_cor', [12; NaN; 3.3; -0.4]);

   [Data, metadata] = icemodel.forcing.buildPromiceData( ...
      "KAN_M", source_dir=root);

   testCase.verifyEqual(Data.swd(1), 11);
   testCase.verifyEqual(Data.swu(1), 12);
   testCase.verifyTrue(isnan(Data.swn(1)));
   testCase.verifyEqual(Data.swn(2), 18);
   testCase.verifyEqual(metadata.swn_negative_invalid_count, 1);
   testCase.verifyTrue(contains(metadata.swn_policy, ...
      "physically inconsistent"));
end

function test_public_shortwave_prefers_corrected_product_at_zac_sites(testCase)
   % Raw ZAC pyranometer data contain source-valid negative offsets, whereas
   % pypromice's corrected product is nonnegative. Met, Data, and therefore the
   % observations.mat source product use one corrected-first public series.
   sites = ["ZAC_L", "ZAC_U"];
   starts = [datetime(2008, 3, 29, 15, 0, 0), ...
      datetime(2012, 8, 17, 0, 0, 0)];
   ends = [datetime(2008, 4, 1, 23, 0, 0), ...
      datetime(2012, 8, 20, 23, 0, 0)];
   for k = 1:numel(sites)
      aws = icemodel.forcing.readPromiceAws(sites(k), ...
         source_dir=testCase.TestData.source_dir, ...
         startdate=starts(k), enddate=ends(k));
      [met, met_metadata] = icemodel.forcing.buildPromiceMet(sites(k), ...
         source_dir=testCase.TestData.source_dir, ...
         startdate=starts(k), enddate=ends(k));
      [Data, data_metadata] = icemodel.forcing.buildPromiceData(sites(k), ...
         source_dir=testCase.TestData.source_dir, ...
         startdate=starts(k), enddate=ends(k));

      % Construct the source-backed rule independently from the public builder:
      % corrected finite values, then raw fallback, then nonnegative invariant.
      expected = aws.swd_cor;
      fallback = ~isfinite(expected) & isfinite(aws.swd);
      expected(fallback) = aws.swd(fallback);
      expected(isfinite(expected) & expected < 0) = 0;

      testCase.verifyGreaterThan(nnz(aws.swd < 0), 0);
      testCase.verifyGreaterThan(nnz(isfinite(aws.swd_cor)), 0);
      present = isfinite(expected);
      testCase.verifyEqual(Data.swd(present), expected(present));
      testCase.verifyEqual(met.swd(present), expected(present));
      darkness_filled = ~present & isfinite(Data.swd);
      testCase.verifyTrue(all(Data.swd(darkness_filled) == 0));
      testCase.verifyEqual(met.swd, Data.swd);
      testCase.verifyEqual(nnz(darkness_filled), ...
         data_metadata.swd_darkness_zero_filled_count);
      testCase.verifyEqual(data_metadata.swd_corrected_used_count, ...
         nnz(isfinite(aws.swd_cor)));
      testCase.verifyTrue(all(Data.swd(isfinite(Data.swd)) >= 0));
      testCase.verifyTrue(all(met.swd(isfinite(met.swd)) >= 0));
      % Match the public derived-flux policy: retain the source components but
      % exclude a physically inconsistent negative net shortwave sample.
      expected_swn = Data.swd - Data.swu;
      expected_swn(isfinite(expected_swn) & expected_swn < 0) = NaN;
      testCase.verifyEqual(Data.swn, expected_swn);
      testCase.verifyGreaterThan(data_metadata.swd_raw_negative_count, 0);
      testCase.verifyGreaterThan(data_metadata.swd_corrected_used_count, 0);
      testCase.verifyEqual(met_metadata.swd_policy, ...
         data_metadata.swd_policy);
      testCase.verifyEqual(Data.Properties.UserData.swd_policy, ...
         data_metadata.swd_policy);
   end
end

function test_missing_promice_radiation_has_explicit_placeholder_metadata(testCase)
   % LYN_L omits dsr/dsr_cor and albedo. The met contract retains all-NaN
   % placeholders with channel-specific policies, while observational Data
   % omits those absent source channels and records the same provenance.
   [met, met_metadata] = icemodel.forcing.buildPromiceMet("LYN_L", ...
      source_dir=testCase.TestData.source_dir);
   [Data, data_metadata] = icemodel.forcing.buildPromiceData("LYN_L", ...
      source_dir=testCase.TestData.source_dir);

   testCase.verifyTrue(all(isnan(met.swd)));
   testCase.verifyTrue(all(isnan(met.albedo)));
   testCase.verifyTrue(all(contains(string(met_metadata.swd_policy), ...
      ["swd", "placeholder"])));
   testCase.verifyTrue(all(contains(string(met_metadata.albedo_policy), ...
      ["albedo", "placeholder"])));
   testCase.verifyFalse(met_metadata.swd_source_present);
   testCase.verifyFalse(met_metadata.swd_observations_present);
   testCase.verifyFalse(met_metadata.swd_source_file_present);
   testCase.verifyFalse(met_metadata.swd_corrected_source_file_present);
   testCase.verifyFalse(met_metadata.swd_source_file_observations_present);
   testCase.verifyFalse(met_metadata.albedo_source_present);
   testCase.verifyFalse(met_metadata.albedo_observations_present);
   testCase.verifyFalse(ismember("swd", ...
      string(Data.Properties.VariableNames)));
   testCase.verifyFalse(ismember("albedo", ...
      string(Data.Properties.VariableNames)));
   testCase.verifyTrue(contains(string(data_metadata.swd_policy), ...
      "placeholder"));
   testCase.verifyTrue(contains(string(data_metadata.albedo_policy), ...
      "placeholder"));
end

function test_readPromiceAws_clamps_ice_temperature_string(testCase)
   % The tice string is clamped to the dictionary physical range (-80..1 C),
   % so no out-of-range thermistor spikes reach the output.

   aws = icemodel.forcing.readPromiceAws("KAN_L", ...
      source_dir=testCase.TestData.source_dir);
   Tf = icemodel.physicalConstant('Tf');
   names = string(aws.Properties.VariableNames);
   tice = names(~cellfun('isempty', ...
      regexp(cellstr(names), '^tice\d+$', 'once')));
   testCase.verifyNotEmpty(tice);
   for tv = tice
      degc = aws.(tv) - Tf;
      degc = degc(isfinite(degc));
      testCase.verifyTrue(all(degc >= -80 - 1e-9 & degc <= 1 + 1e-9), ...
         sprintf('%s outside the dictionary clamp range', tv));
   end
end

function test_buildPromiceMet_satisfies_met_contract(testCase)
   % A full-year KAN_M met build passes the met contract with a complete
   % hourly axis and the canonical m s-1 ppt unit.

   [met, metadata] = icemodel.forcing.buildPromiceMet("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 1, 1), ...
      enddate=datetime(2015, 12, 31, 23, 0, 0));

   icemodel.forcing.helpers.validatemet(met)
   testCase.verifyEqual(height(met), 8760);
   testCase.verifyTrue(all(isnan(met.rainf)));
   testCase.verifyTrue(all(isnan(met.snowf)));
   testCase.verifyTrue(all(isnan(met.ppt)));
   testCase.verifyTrue(all(isfinite(met.tair)));
   testCase.verifyTrue(all(met.albedo >= 0.05 & met.albedo <= 0.98));
   source_info = dir(string(metadata.source_file));
   testCase.verifyEqual(double(metadata.source_size_bytes), ...
      source_info.bytes);
   testCase.verifyEqual(string(metadata.source_sha256), ...
      icemodel.verification.setup.fileSha256(metadata.source_file));
   aws = icemodel.forcing.readPromiceAws("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 1, 1), ...
      enddate=datetime(2015, 12, 31, 23, 0, 0));
   testCase.verifyTrue(ismember("boom_height", ...
      string(met.Properties.VariableNames)));
   [~, expected_swu] = icemodel.forcing.helpers.promiceShortwave(aws, ...
      fill_darkness=true, latitude=metadata.lat, ...
      longitude=metadata.lon, ...
      swd_source_file_observations_present= ...
      metadata.swd_source_file_observations_present, ...
      swu_source_file_observations_present= ...
      metadata.swu_source_file_observations_present);
   testCase.verifyEqual(met.swu, expected_swu);
   testCase.verifyEqual(met.boom_height, aws.boom_height);
   testCase.verifyTrue(isfield(metadata, 'station_transition_times'));
   testCase.verifyTrue(isfield(metadata, 'station_transition_record'));

   [~, ~, pptunit] = icemodel.forcing.helpers.metvariables();
   units = string(met.Properties.VariableUnits);
   names = string(met.Properties.VariableNames);
   testCase.verifyEqual(units(names == "rainf"), pptunit);
   testCase.verifyEqual(units(names == "snowf"), pptunit);
   testCase.verifyEqual(units(names == "ppt"), pptunit);
   testCase.verifyEqual(string(met.Properties.UserData.station_name), ...
      "KAN_M");
   testCase.verifyEqual(string(met.Properties.UserData.site), "kanm");
   testCase.verifyEqual(met.Properties.UserData.met_variables, ...
      string(met.Properties.VariableNames(:)));
   testCase.verifyFalse(met.Properties.UserData.rainf_source_present);
   testCase.verifyFalse(met.Properties.UserData.rainf_observations_present);
   testCase.verifyTrue(contains( ...
      met.Properties.UserData.rainf_policy, "source channel absent"));
end

function test_buildPromiceMet_uses_selected_station_inventory(testCase)
   % A custom source root owns both its NetCDF and composing-station inventory;
   % repository metadata must not add handovers from another product snapshot.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   writeSyntheticPromice(root, [0; 0.25; 0.5; 0.75]);

   site_header = "site_id,project,location_type,stations," + ...
      "date_installation,latitude_installation,longitude_installation," + ...
      "altitude_installation,date_last_valid,latitude_last_valid," + ...
      "longitude_last_valid,altitude_last_valid";
   writelines([site_header; ...
      "KAN_M,PROMICE,ice sheet,CUSTOM_A CUSTOM_B,2020-01-01," + ...
      "67,-48,1200,2020-01-02,67,-48,1200"], ...
      fullfile(root, 'AWS_sites_metadata.csv'));
   station_header = "station_id,project,number_of_booms,location_type," + ...
      "date_installation,latitude_installation,longitude_installation," + ...
      "altitude_installation,date_last_valid,latitude_last_valid," + ...
      "longitude_last_valid,altitude_last_valid";
   writelines([station_header; ...
      "CUSTOM_A,PROMICE,1,ice sheet,2020-01-01,67,-48,1200," + ...
      "2020-01-02,67,-48,1200"; ...
      "CUSTOM_B,PROMICE,1,ice sheet,2020-01-02,67,-48,1200," + ...
      "2020-01-03,67,-48,1200"], ...
      fullfile(root, 'AWS_stations_metadata.csv'));

   [~, metadata] = icemodel.forcing.buildPromiceMet("KAN_M", ...
      source_dir=root);

   testCase.verifyEqual(string(metadata.composing_stations), ...
      ["CUSTOM_A"; "CUSTOM_B"]);
   testCase.verifyEqual(string({metadata.station_transition_record.station}), ...
      ["CUSTOM_A", "CUSTOM_B"]);
end

function test_buildPromiceMet_does_not_escape_custom_source_root(testCase)
   % A custom product root without local inventories must not borrow
   % composing-station or transition metadata from its unrelated parent.
   outer = string(tempname);
   selected = fullfile(outer, "selected");
   mkdir(fullfile(selected, "hour"))
   cleanup = onCleanup(@() rmdir(outer, 's'));
   writeSyntheticPromice(selected, [0; 0.25; 0.5; 0.75]);

   site_header = "site_id,project,location_type,stations," + ...
      "date_installation,latitude_installation,longitude_installation," + ...
      "altitude_installation,date_last_valid,latitude_last_valid," + ...
      "longitude_last_valid,altitude_last_valid";
   writelines([site_header; ...
      "KAN_M,OTHER,ice sheet,WRONG_A WRONG_B,2020-01-01," + ...
      "67,-48,1200,2020-01-02,67,-48,1200"], ...
      fullfile(outer, 'AWS_sites_metadata.csv'));
   station_header = "station_id,project,number_of_booms,location_type," + ...
      "date_installation,latitude_installation,longitude_installation," + ...
      "altitude_installation,date_last_valid,latitude_last_valid," + ...
      "longitude_last_valid,altitude_last_valid";
   writelines([station_header; ...
      "WRONG_A,OTHER,1,ice sheet,2020-01-01,67,-48,1200," + ...
      "2020-01-02,67,-48,1200"; ...
      "WRONG_B,OTHER,1,ice sheet,2020-01-02,67,-48,1200," + ...
      "2020-01-03,67,-48,1200"], ...
      fullfile(outer, 'AWS_stations_metadata.csv'));

   [~, metadata] = icemodel.forcing.buildPromiceMet("KAN_M", ...
      source_dir=selected);

   testCase.verifyEqual(string(metadata.composing_stations), "KAN_M");
   testCase.verifyEmpty(metadata.station_transition_record);
end

function test_buildPromiceMet_preserves_long_station_outage(testCase)
   % CEN has a documented multi-year swd/tair outage beginning at this hour;
   % neither builder QA nor the shared 15-minute writer may bridge it.
   startdate = datetime(2007, 5, 3, 17, 0, 0, TimeZone="UTC");
   enddate = datetime(2007, 5, 4, 20, 0, 0, TimeZone="UTC");
   [met, metadata] = icemodel.forcing.buildPromiceMet("CEN", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=startdate, enddate=enddate);

   testCase.verifyTrue(any(isfinite(met.swd)));
   testCase.verifyTrue(any(isnan(met.swd)));
   testCase.verifyTrue(any(isnan(met.tair)));
   testCase.verifyTrue(contains(string(metadata.gap_policy), ...
      "other source gaps preserved"));

   outdir = string(tempname);
   mkdir(outdir)
   testCase.addTeardown(@() rmdir(outdir, 's'));
   filename = icemodel.forcing.helpers.writemet(met, "cen", "promice", ...
      outdir=outdir);
   loaded = load(filename, 'met', 'artifact_metadata');
   outage = loaded.met.Time >= datetime(2007, 5, 3, 19, 0, 0, ...
      TimeZone="UTC");

   testCase.verifyTrue(all(isnan(loaded.met.swd(outage))));
   testCase.verifyTrue(all(isnan(loaded.met.tair(outage))));
   testCase.verifyEqual(string( ...
      loaded.artifact_metadata.met_resample_policy), ...
      "interval_start_zero_order_hold");
end

function test_promice_builders_persist_deep_night_zero_and_outage_nan(testCase)
   % A KAN_L outage spanning night and day receives derived zeros only in deep
   % civil night. An outage-only surgical build must match a broader build at
   % the same timestamps, and both persistence boundaries retain exact counts.
   broad_start = datetime(2008, 11, 8, 0, 0, 0, TimeZone="UTC");
   outage_start = datetime(2008, 11, 9, 0, 0, 0, TimeZone="UTC");
   outage_end = datetime(2008, 11, 9, 23, 0, 0, TimeZone="UTC");
   [broad_data, broad_data_metadata] = ...
      icemodel.forcing.buildPromiceData("KAN_L", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=broad_start, enddate=outage_end);
   [broad_met, broad_met_metadata] = ...
      icemodel.forcing.buildPromiceMet("KAN_L", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=broad_start, enddate=outage_end);
   [Data, data_metadata] = icemodel.forcing.buildPromiceData("KAN_L", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=outage_start, enddate=outage_end);
   [met, met_metadata] = icemodel.forcing.buildPromiceMet("KAN_L", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=outage_start, enddate=outage_end);

   % Whole-source support makes output values independent of the requested
   % window even though source-selection counts remain local to each artifact.
   broad_rows = broad_data.Time >= outage_start;
   testCase.verifyEqual(Data.Time, broad_data.Time(broad_rows));
   testCase.verifyEqual(Data.swd, broad_data.swd(broad_rows));
   testCase.verifyEqual(Data.swu, broad_data.swu(broad_rows));
   testCase.verifyEqual(met.Time, broad_met.Time(broad_rows));
   testCase.verifyEqual(met.swd, broad_met.swd(broad_rows));
   testCase.verifyEqual(broad_data_metadata.swd_darkness_zero_filled_count, 14);
   testCase.verifyEqual(broad_met_metadata.swd_darkness_zero_filled_count, 14);

   % The complete 9 November source-radiation outage contains 14 conservative
   % deep-night hours; its ten twilight/daylight hours remain unavailable.
   deep_night = ismember(hour(Data.Time), [0:9, 20:23]);
   other_outage = ~deep_night;
   testCase.verifyEqual(nnz(deep_night), 14);
   testCase.verifyTrue(all(Data.swd(deep_night) == 0));
   testCase.verifyTrue(all(Data.swu(deep_night) == 0));
   testCase.verifyTrue(all(isnan(Data.swd(other_outage))));
   testCase.verifyTrue(all(isnan(Data.swu(other_outage))));
   testCase.verifyEqual(met.swd, Data.swd);
   testCase.verifyEqual(data_metadata.swd_darkness_zero_filled_count, 14);
   testCase.verifyEqual(data_metadata.swu_darkness_zero_filled_count, 14);
   testCase.verifyEqual(met_metadata.swd_darkness_zero_filled_count, 14);
   testCase.verifyEqual(data_metadata.swd_remaining_missing_count, 10);
   testCase.verifyFalse(data_metadata.swd_observations_present);
   testCase.verifyTrue(data_metadata.swd_source_file_observations_present);
   testCase.verifyTrue(met_metadata.swd_source_file_observations_present);

   % Exercise the actual persistence boundaries used by PROMICE imports.
   outdir = string(tempname);
   mkdir(outdir)
   testCase.addTeardown(@() rmdir(outdir, 's'));
   userdata_file = icemodel.forcing.helpers.writeuserdata(Data, ...
      "kanl", "promice", outdir=outdir);
   met_file = icemodel.forcing.helpers.writemet(met, ...
      "kanl", "promice", outdir=outdir);
   saved_data = load(userdata_file(1), 'Data', 'artifact_metadata');
   saved_met = load(met_file, 'met', 'artifact_metadata');
   testCase.verifyEqual(saved_data.Data.swd, Data.swd);
   testCase.verifyEqual(saved_data.Data.swu, Data.swu);

   % Source-selection, clamp, darkness, and remaining-gap counts are source-hour
   % quantities. Writers preserve them exactly in table and top-level metadata.
   count_suffixes = ["_corrected_used_count", "_raw_fallback_count", ...
      "_raw_negative_count", "_corrected_negative_count", ...
      "_raw_fallback_negative_count", "_negative_clamped_count", ...
      "_darkness_zero_filled_count", "_remaining_missing_count"];
   for channel = ["swd", "swu"]
      for suffix = count_suffixes
         field = char(channel + suffix);
         expected = data_metadata.(field);
         testCase.verifyEqual(met_metadata.(field), expected);
         testCase.verifyEqual( ...
            saved_data.Data.Properties.UserData.(field), expected);
         testCase.verifyEqual(saved_data.artifact_metadata.(field), expected);
         testCase.verifyEqual( ...
            saved_met.met.Properties.UserData.(field), expected);
         testCase.verifyEqual(saved_met.artifact_metadata.(field), expected);
      end
   end

   % The default writer emits four 15-minute samples for each source hour.
   dark_interval = saved_met.met.Time >= datetime(2008, 11, 9, 0, 0, 0, ...
      TimeZone="UTC") & saved_met.met.Time < datetime(2008, 11, 9, 1, 0, 0, ...
      TimeZone="UTC");
   day_interval = saved_met.met.Time >= datetime(2008, 11, 9, 12, 0, 0, ...
      TimeZone="UTC") & saved_met.met.Time < datetime(2008, 11, 9, 13, 0, 0, ...
      TimeZone="UTC");
   testCase.verifyEqual(nnz(dark_interval), 4);
   testCase.verifyEqual(nnz(day_interval), 4);
   testCase.verifyTrue(all(saved_met.met.swd(dark_interval) == 0));
   testCase.verifyTrue(all(isnan(saved_met.met.swd(day_interval))));
   testCase.verifyEqual(string( ...
      saved_met.artifact_metadata.met_resample_policy), ...
      "interval_start_zero_order_hold");
end

function test_buildPromiceMet_preserves_liquid_rain_only(testCase)
   % NUK_B is a short representative record with rainfall_cor_u. The builder
   % converts the corrected hourly LIQUID amount to the canonical rate without
   % inventing solid or total precipitation from the tipping-bucket channel.

   startdate = datetime(2023, 10, 6, 9, 0, 0);
   enddate = datetime(2023, 10, 6, 12, 0, 0);
   aws = icemodel.forcing.readPromiceAws("NUK_B", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=startdate, enddate=enddate);
   [met, metadata] = icemodel.forcing.buildPromiceMet("NUK_B", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=startdate, enddate=enddate);

   observed = isfinite(aws.rainf);
   testCase.verifyTrue(any(aws.rainf > 0));
   testCase.verifyEqual(met.rainf(observed), ...
      aws.rainf(observed) * 1e-3 / 3600, 'AbsTol', 1e-15);
   testCase.verifyTrue(all(isnan(met.snowf)));
   testCase.verifyTrue(all(isnan(met.ppt)));
   testCase.verifyTrue(metadata.rainf_source_present);
   testCase.verifyTrue(metadata.rainf_observations_present);
   testCase.verifyTrue(contains(metadata.rainf_policy, "converted"));
   testCase.verifyTrue(contains(metadata.precip_policy, "liquid"));
   testCase.verifyTrue(contains(metadata.precip_policy, ...
      "not reliable solid precipitation"));
end

function test_buildPromiceMet_uses_observed_lwd_by_default(testCase)
   % A station with a real longwave channel (KAN_M) reports lwd as observed:
   % lwd_estimated is false and the empirical fallback is not engaged.

   [met, metadata] = icemodel.forcing.buildPromiceMet("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 6, 1), ...
      enddate=datetime(2015, 6, 2));

   testCase.verifyFalse(metadata.lwd_estimated);
   testCase.verifyFalse(metadata.lwd_placeholder);
   testCase.verifyTrue(contains(metadata.lwd_policy, "observed"));
   testCase.verifyTrue(any(isfinite(met.lwd)));
end

function test_buildPromiceMet_fill_lwd_does_not_override_observed(testCase)
   % fill_lwd is an OPT-IN fallback for stations with NO longwave sensor; it
   % must never override a real observed channel. With fill_lwd=true at a
   % station that has lwd (KAN_M), lwd stays observed (lwd_estimated false).

   [~, metadata] = icemodel.forcing.buildPromiceMet("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 6, 1), ...
      enddate=datetime(2015, 6, 2), fill_lwd=true);

   testCase.verifyFalse(metadata.lwd_estimated);
   testCase.verifyFalse(metadata.lwd_placeholder);
end

function test_buildPromiceMet_defaults_missing_lwd_to_placeholder(testCase)
   % Missing required channels are explicit placeholders by default. LYN_L has
   % no observed longwave channel, but the builder should still produce a valid
   % met-contract timetable whose lwd can be filled by a later source swap.
   lyn_file = fullfile(testCase.TestData.source_dir, 'hour', 'LYN_L_hour.nc');
   testCase.assumeTrue(isfile(lyn_file), ...
      'LYN_L PROMICE source not available; skipping missing-lwd test.');

   [met, metadata] = icemodel.forcing.buildPromiceMet("LYN_L", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2020, 1, 1), ...
      enddate=datetime(2020, 1, 2));

   icemodel.forcing.helpers.validatemet(met)
   testCase.verifyTrue(all(isnan(met.lwd)));
   testCase.verifyTrue(all(isnan(met.rainf)));
   testCase.verifyTrue(all(isnan(met.snowf)));
   testCase.verifyTrue(all(isnan(met.ppt)));
   testCase.verifyTrue(metadata.lwd_placeholder);
   testCase.verifyTrue(contains(metadata.lwd_policy, "placeholder"));
   testCase.verifyFalse(metadata.rainf_source_present);
   testCase.verifyFalse(metadata.rainf_observations_present);
   testCase.verifyTrue(contains(metadata.rainf_policy, ...
      "source channel absent"));
end

function test_buildPromiceMet_strict_mode_rejects_missing_lwd(testCase)
   % fillwithmissing=false remains the explicit strict mode for callers that
   % want absent required channels to abort instead of becoming placeholders.
   lyn_file = fullfile(testCase.TestData.source_dir, 'hour', 'LYN_L_hour.nc');
   testCase.assumeTrue(isfile(lyn_file), ...
      'LYN_L PROMICE source not available; skipping strict-mode test.');

   f = @() icemodel.forcing.buildPromiceMet("LYN_L", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2020, 1, 1), ...
      enddate=datetime(2020, 1, 2), fillwithmissing=false);

   testCase.verifyError(f, ...
      'icemodel:forcing:buildPromiceMet:missingForcing');
end

function test_buildPromiceMet_strict_recognizes_derivable_shortwave(testCase)
   % Corrected downwelling shortwave satisfies strict source presence without
   % raw dsr, while absent upward shortwave remains derivable downstream.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   writeSyntheticPromice(root, [0; 0.25; 0.5; 0.75], false, false);

   caught = [];
   try
      icemodel.forcing.buildPromiceMet("KAN_M", source_dir=root, ...
         fillwithmissing=false);
   catch exception
      caught = exception;
   end

   testCase.assertNotEmpty(caught);
   testCase.verifyEqual(caught.identifier, ...
      'icemodel:forcing:buildPromiceMet:missingForcing');
   testCase.verifyTrue(contains(caught.message, 'ppt'));
   testCase.verifyFalse(contains(caught.message, 'swd'));
   testCase.verifyFalse(contains(caught.message, 'swu'));
end

function test_empirical_lwd_estimate_is_physical(testCase)
   % The empirical fallback the fill_lwd path calls returns finite, physically
   % plausible downwelling longwave for typical polar air temperature + vapor
   % pressure when a caller chooses estimated lwd instead of placeholders.

   tair = (250:5:270)';             % K
   ea = 100 + zeros(size(tair));    % Pa, low-humidity polar air
   lwd = icemodel.surface.empirical_incoming_longwave_radiation(tair, ea);

   testCase.verifyTrue(all(isfinite(lwd)));
   testCase.verifyTrue(all(lwd > 50 & lwd < 350));   % plausible polar range
   testCase.verifyTrue(issorted(lwd));               % warmer sky -> more lwd
end

function test_buildPromiceData_reads_l3_evaluation_channels(testCase)
   % The Data builder reads the QC'd L3 surface channels (snow depth from
   % snow_height, ablation from z_ice_surf) rather than deriving them, with
   % non-negative snow depth and physically plausible monotone ablation over
   % the bare-ice ablation site KAN_L (2009-2022).

   [Data, metadata] = icemodel.forcing.buildPromiceData("KAN_L", ...
      source_dir=testCase.TestData.source_dir, frequency="daily", ...
      startdate=datetime(2009, 1, 1), ...
      enddate=datetime(2022, 12, 31, 23, 0, 0));

   testCase.verifyEqual(string(metadata.site_surface_type), "ablation");
   testCase.verifyEqual(string(Data.Properties.UserData.station_name), ...
      "KAN_L");
   testCase.verifyEqual(size(Data.Properties.UserData.composing_stations, 2), ...
      1);
   testCase.verifyEqual(string(metadata.surface_channel), "ablation");
   testCase.verifyEqual(string(metadata.surface_source), "L3 z_ice_surf");
   testCase.verifyEqual(string(metadata.ablation_source), "L3 z_ice_surf");
   testCase.verifyTrue(startsWith(string(metadata.snow_depth_source), ...
      "L3 snow_height"));

   % Snow depth: clamped strictly non-negative per the readme (bare-ice
   % ablation site stays near zero).
   sd = Data.snow_depth(isfinite(Data.snow_depth));
   testCase.verifyGreaterThan(numel(sd), 1000);
   testCase.verifyTrue(all(sd >= 0));

   % Gap flag present and binary on the surface-height-derived channel.
   testCase.verifyTrue(ismember("surface_height_flag", ...
      string(Data.Properties.VariableNames)));
   fl = Data.surface_height_flag(isfinite(Data.surface_height_flag));
   testCase.verifyTrue(all(fl == 0 | fl == 1));

   % Ablation: cumulative surface lowering, positive-down, monotone-ish, with
   % a physically plausible multi-year magnitude (KAN_L lowers ~57 m over the
   % window; the homegrown derivation inflated multi-year totals).
   ab = Data.ablation(isfinite(Data.ablation));
   testCase.verifyGreaterThan(numel(ab), 1000);
   testCase.verifyEqual(ab(1), 0, 'AbsTol', 1e-6);   % zeroed at window start
   testCase.verifyGreaterThan(max(ab), 30);
   testCase.verifyLessThan(max(ab), 80);
   testCase.verifyGreaterThan(mean(diff(ab) >= 0), 0.95);   % near-monotone

   % Ice-temperature string and evolving depths read from t_i_*/d_t_i_*;
   % source/canonical/flag diagnostics make the primary 10 m target auditable.
   dn = string(Data.Properties.VariableNames);
   testCase.verifyTrue(all(ismember( ...
      "tice" + string(1:8), dn)));
   testCase.verifyTrue(all(ismember( ...
      "dtice" + string(1:8), dn)));
   testCase.verifyTrue(all(ismember( ...
      ["tice10m", "tice10m_source", "tice10m_qc_flag"], dn)));
   flagged_tice = Data.tice10m_qc_flag > 0;
   testCase.verifyTrue(all(isnan(Data.tice10m(flagged_tice))));
   testCase.verifyEqual(metadata.tice10m_qc_flagged_sample_count, ...
      nnz(flagged_tice));

   % Derived net fluxes present.
   testCase.verifyTrue(all(ismember(["swn", "lwn", "netr", "thf"], dn)));

   % Userdata metadata contract.
   testCase.verifyEqual(Data.Properties.CustomProperties.Lat, 67.097, ...
      'AbsTol', 1e-2);
   testCase.verifyEqual(numel(Data.Properties.CustomProperties.ScalarUnits), 6);
end

function test_buildPromiceData_accumulation_site_uses_surface_height(testCase)
   % An accumulation-zone station with no z_ice_surf (KAN_U) emits a
   % surface_height channel from z_surf_combined (NET accumulation, positive
   % up), NOT an ablation or snow_depth channel. The gap flag rides along.

   [Data, metadata] = icemodel.forcing.buildPromiceData("KAN_U", ...
      source_dir=testCase.TestData.source_dir, frequency="daily", ...
      startdate=datetime(2009, 1, 1), ...
      enddate=datetime(2022, 12, 31, 23, 0, 0));

   testCase.verifyEqual(string(metadata.site_surface_type), "accumulation");
   testCase.verifyEqual(string(metadata.surface_channel), "surface_height");
   testCase.verifyEqual(string(metadata.surface_source), ...
      "L3 z_surf_combined");

   names = string(Data.Properties.VariableNames);
   testCase.verifyTrue(ismember("surface_height", names));
   testCase.verifyTrue(ismember("surface_height_flag", names));
   % No ablation / snow_depth fabricated for an accumulation site.
   testCase.verifyFalse(ismember("ablation", names));
   testCase.verifyFalse(ismember("snow_depth", names));

   sh = Data.surface_height(isfinite(Data.surface_height));
   testCase.verifyGreaterThan(max(sh), 0);         % surface rose (accumulation)
   testCase.verifyLessThan(max(abs(sh)), 15);      % plausible magnitude
end

function test_buildPromiceData_keeps_wholly_masked_tice10m_contract(testCase)
   % A surgical window inside KAN_U's unresolved 2022 thermistor epoch must keep
   % the all-NaN canonical channel beside its finite source and code-3 flag.
   [Data, metadata] = icemodel.forcing.buildPromiceData("KAN_U", ...
      source_dir=testCase.TestData.source_dir, frequency="hourly", ...
      startdate=datetime(2022, 9, 4), enddate=datetime(2022, 9, 6, 23, 0, 0));

   names = string(Data.Properties.VariableNames);
   testCase.verifyTrue(all(ismember( ...
      ["tice10m", "tice10m_source", "tice10m_qc_flag"], names)));
   testCase.verifyTrue(all(isnan(Data.tice10m)));
   testCase.verifyTrue(all(isfinite(Data.tice10m_source)));
   testCase.verifyTrue(all(Data.tice10m_qc_flag == 3));
   testCase.verifyEqual(metadata.tice10m_qc_persistent_sample_count, ...
      height(Data));
end

function test_buildPromiceData_gap_flag_from_sensors_not_just_z_nan(testCase)
   % The improved gap flag is sensor-derived: it must catch slope-bridged
   % samples (all surface sensors NaN but z finite) that the old z-NaN-only
   % heuristic missed. MIT has thousands of such samples, so the new
   % gap-flagged count must EXCEED the bare z-NaN count.

   [Data, metadata] = icemodel.forcing.buildPromiceData("MIT", ...
      source_dir=testCase.TestData.source_dir, frequency="hourly");

   % Reconstruct the bare z-NaN count from the raw reader for the comparison.
   aws = icemodel.forcing.readPromiceAws("MIT", ...
      source_dir=testCase.TestData.source_dir);
   z_nan_only = nnz(~isfinite(aws.z_ice_surf));

   names = string(Data.Properties.VariableNames);
   testCase.verifyTrue(ismember("surface_height_flag", names));
   % The sensor-derived gap count must strictly exceed the z-NaN-only count
   % (slope-bridged segments are now flagged too).
   testCase.verifyGreaterThan(metadata.gap_flagged_samples, z_nan_only);
end

function test_buildPromiceData_gap_flag_does_not_overlay_observations(testCase)
   % A gap-flagged sample must coincide with all underlying surface sensors
   % being NaN (no direct observation), so the red gap markers never overlay a
   % genuinely-observed sample.

   Data = icemodel.forcing.buildPromiceData("MIT", ...
      source_dir=testCase.TestData.source_dir, frequency="hourly");
   aws = icemodel.forcing.readPromiceAws("MIT", ...
      source_dir=testCase.TestData.source_dir);

   sensor_names = intersect(["transducer_depth", "boom_height", ...
      "stake_height"], string(aws.Properties.VariableNames), 'stable');
   sensors = cell2mat(arrayfun(@(v) aws.(v), sensor_names, ...
      'UniformOutput', false));
   all_sensors_nan = all(~isfinite(sensors), 2);

   % Interior gap samples (z finite) must have all sensors NaN: no observed
   % sensor reading hides under a gap flag.
   gap = Data.surface_height_flag == 1;
   z_finite = isfinite(aws.z_ice_surf);
   interior_gap = gap & z_finite;
   testCase.verifyTrue(all(all_sensors_nan(interior_gap)), ...
      'a gap flag overlays a sample with a finite surface sensor');
end

function test_buildPromiceData_stages_station_transition_flag(testCase)
   % The station-transition flag channel is staged on every site (0/1), with the
   % composing-station provenance (multistation) recorded for the consumer.

   [Data, metadata] = icemodel.forcing.buildPromiceData("KAN_L", ...
      source_dir=testCase.TestData.source_dir, frequency="daily");

   names = string(Data.Properties.VariableNames);
   testCase.verifyTrue(ismember("station_transition_flag", names));
   stf = Data.station_transition_flag(isfinite(Data.station_transition_flag));
   testCase.verifyTrue(all(stf == 0 | stf == 1));
   testCase.verifyTrue(isfield(metadata, 'composing_stations'));
   testCase.verifyTrue(isfield(metadata, 'is_multistation'));
   testCase.verifyTrue(isfield(metadata, 'station_transition_times'));
end

function test_buildPromiceData_populates_station_transition_flag(testCase)
   % For a multi-station site (CEN = CEN2/CEN1/GITS) the per-station install
   % dates in AWS_stations_metadata.csv populate station_transition_flag: each
   % WITHIN-RECORD handover (an install strictly after the record start) opens a
   % flag window, so the flag is NOT all-zero and a window surrounds each
   % handover date. The founding station's install (== record start) is NOT a
   % handover; legacy GC-Net names absent from the CSV contribute no date.

   [Data, metadata] = icemodel.forcing.buildPromiceData("CEN", ...
      source_dir=testCase.TestData.source_dir, frequency="daily");

   testCase.verifyTrue(metadata.is_multistation);

   % CEN has two within-record handovers: CEN1 (2017-07-24), CEN2 (2021-08-12).
   tt = metadata.station_transition_times;
   testCase.verifyEqual(numel(tt), 2);
   testCase.verifyTrue(any(abs(days(tt - ...
      datetime(2017, 7, 24, 'TimeZone', 'UTC'))) < 1));
   testCase.verifyTrue(any(abs(days(tt - ...
      datetime(2021, 8, 12, 'TimeZone', 'UTC'))) < 1));

   % The flag is populated (not all-zero) and binary, with a window around each
   % handover (default +-14 days -> at least a handful of flagged daily samples).
   % metadata.station_transition_samples is the count on the underlying hourly
   % series; the daily Data carries the corresponding daily (max-aggregated)
   % windows, so both are strictly positive but not numerically equal.
   stf = Data.station_transition_flag;
   testCase.verifyTrue(all(stf == 0 | stf == 1));
   testCase.verifyGreaterThan(metadata.station_transition_samples, 0);
   testCase.verifyGreaterThan(nnz(stf == 1), 0);

   % Each handover date sits inside a flagged daily sample.
   for h = tt(:)'
      [~, near] = min(abs(Data.Time - h));
      testCase.verifyEqual(stf(near), 1, ...
         sprintf('no flag window around handover %s', string(h)));
   end

   % Provenance record: a row per composing station; the founding CEN1/CEN2 v3
   % ids are in the CSV, the legacy GC-Net "GITS" is not.
   rec = metadata.station_transition_record;
   testCase.verifyEqual(numel(rec), numel(metadata.composing_stations));
   in_csv = [rec.in_csv];
   testCase.verifyTrue(any(in_csv));    % some composing stations resolved
   testCase.verifyTrue(any(~in_csv));   % GITS legacy name absent from CSV
end

function test_buildPromiceData_transition_evidence_wired_to_destep(testCase)
   % A station-transition handover supplies the 'transition' evidence line to the
   % de-stepping detector: re-running destepSurface on the staged surface channel
   % with the staged transition_times reproduces 'transition' evidence on at
   % least one detected step coincident with a handover, strengthening its
   % classification. (CEN is accumulation -> surface_height channel.)

   [Data, metadata] = icemodel.forcing.buildPromiceData("CEN", ...
      source_dir=testCase.TestData.source_dir, frequency="hourly");

   tt = metadata.station_transition_times;
   testCase.assumeNotEmpty(tt);

   [~, record] = icemodel.forcing.destepSurface(Data.Time, ...
      Data.surface_height, mode="detect", transition_times=tt, ...
      season="accumulation");

   % With handover times supplied, at least one detected step near a handover
   % carries the 'transition' evidence line.
   has_transition = false;
   for r = record
      if any(r.evidence == "transition")
         has_transition = true;
         break
      end
   end
   testCase.verifyTrue(has_transition, ...
      'no detected step gained the station-transition evidence line');
end

function test_buildPromiceData_stages_step_flags_unaltered(testCase)
   % The de-stepping DETECTION is staged (step_detected/correctable + signed
   % magnitude) but the staged ablation series itself is UNALTERED: it must
   % still equal the raw -(z - z(start)) lowering, proving correction is not
   % baked into the staged data.

   [Data, metadata] = icemodel.forcing.buildPromiceData("MIT", ...
      source_dir=testCase.TestData.source_dir, frequency="hourly");

   names = string(Data.Properties.VariableNames);
   testCase.verifyTrue(all(ismember(["step_detected_flag", ...
      "step_correctable_flag", "step_magnitude"], names)));

   % MIT carries unambiguous (correctable) steps (~5.9 m installation jumps).
   testCase.verifyGreaterThan(metadata.steps_correctable, 0);

   % Staged ablation is faithful: equals the raw lowering, NOT de-stepped.
   aws = icemodel.forcing.readPromiceAws("MIT", ...
      source_dir=testCase.TestData.source_dir);
   z = aws.z_ice_surf;
   raw_ablation = -(z - z(find(isfinite(z), 1)));
   testCase.verifyEqual(Data.ablation, raw_ablation, 'AbsTol', 1e-9);

   % The raw series still carries the ~11.9 m of bogus installation jump, so it
   % overshoots the de-stepped magnitude (the correction is NOT in the staged
   % data).
   ab = raw_ablation(isfinite(raw_ablation));
   testCase.verifyGreaterThan(max(ab), 40);   % raw includes the bogus jump
end

function test_buildPromiceData_units_from_shared_map(testCase)
   % The Data output carries the shared canonical units (ablation/snow_depth
   % in m, tice in K).

   Data = icemodel.forcing.buildPromiceData("kanm", ...
      source_dir=testCase.TestData.source_dir, frequency="daily", ...
      startdate=datetime(2015, 1, 1), enddate=datetime(2015, 12, 31));

   names = string(Data.Properties.VariableNames);
   units = string(Data.Properties.VariableUnits);
   testCase.verifyEqual(units(names == "ablation"), "m");
   testCase.verifyEqual(units(names == "snow_depth"), "m");
   testCase.verifyEqual(units(names == "surface_height_flag"), "1");
   testCase.verifyEqual(units(names == "tice1"), "K");
   if ismember("dtice1", names)
      testCase.verifyEqual(units(names == "dtice1"), "m");
   end
   if ismember("tice10m", names)
      testCase.verifyEqual(units(names == "tice10m"), "K");
      testCase.verifyEqual(units(names == "tice10m_source"), "K");
      testCase.verifyEqual(units(names == "tice10m_qc_flag"), "1");
   end
   testCase.verifyEqual(units(names == "tair"), "K");
end

function test_readPromiceAws_shum_kg_per_kg(testCase)
   % qh_u is g/kg in the L3 product (mislabeled kg/kg); the reader converts
   % to kg/kg so shum matches the MAR/MERRA vapor-kernel convention.

   aws = icemodel.forcing.readPromiceAws("KAN_M", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 6, 1), enddate=datetime(2015, 7, 1));
   s = aws.shum(isfinite(aws.shum));
   testCase.verifyNotEmpty(s);
   % kg/kg specific humidity over ice/firn is O(1e-3); g/kg would be O(1).
   testCase.verifyLessThan(median(s), 0.05);
   testCase.verifyGreaterThan(median(s), 1e-5);
end

function test_readPromiceAws_discards_surfaced_thermistors(testCase)
   % A surfaced thermistor (depth <= 0) is discarded: no tice sample survives
   % where its depth tag is at/above the surface.

   aws = icemodel.forcing.readPromiceAws("KAN_L", ...
      source_dir=testCase.TestData.source_dir);
   names = string(aws.Properties.VariableNames);
   tice = names(startsWith(names, "tice") & names ~= "tice10m");
   for tv = tice
      dv = "d" + tv;
      if ~ismember(dv, names); continue; end
      surfaced = aws.(dv) <= 0;
      testCase.verifyTrue(all(isnan(aws.(tv)(surfaced))), ...
         sprintf('%s retained a surfaced sample', tv));
   end
   % tice10m, the primary subsurface channel, is present.
   testCase.verifyTrue(ismember("tice10m", names));
end

function test_readPromiceAws_masks_reviewed_tice10m_discontinuity(testCase)
   % A physically impossible consecutive-hour target jump is masked at both
   % endpoints while the source-derived value and native-depth review survive.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   time = (0:7)';
   source = [-10; -9.8; -9.6; -3; -9.4; -9.3; -9.2; -9.1];
   temperatures = [ ...
      -5.0 -7.0 -9.5
      -5.1 -7.1 -9.4
      -5.2 -7.2 -9.3
      -5.2 -7.1 -3.1
      -5.3 -7.2 -9.2
      -5.4 -7.3 -9.1
      -5.5 -7.4 -9.0
      -5.6 -7.5 -8.9];
   depths = repmat([2, 7, 10.5], numel(time), 1);
   % A small depth-string adjustment at the sensor event proves that native
   % time-varying depths remain part of the review rather than being ignored.
   depths(4:end, :) = depths(4:end, :) + 0.2;
   writeSyntheticThermistors(root, time, source, temperatures, depths);

   [aws, metadata] = icemodel.forcing.readPromiceAws("KAN_U", ...
      source_dir=root);
   Tf = icemodel.physicalConstant('Tf');
   expected_flag = [0; 0; 1; 1; 1; 0; 0; 0];
   testCase.verifyEqual(aws.tice10m_source, source + Tf);
   testCase.verifyEqual(aws.tice10m_qc_flag, expected_flag);
   testCase.verifyTrue(all(isnan(aws.tice10m(expected_flag > 0))));
   testCase.verifyEqual(aws.tice10m(expected_flag == 0), ...
      aws.tice10m_source(expected_flag == 0));
   testCase.verifyEqual(metadata.tice10m_qc_failed_sample_count, 3);
   testCase.verifyEqual(metadata.tice10m_qc_unreviewed_sample_count, 0);
   testCase.verifyEqual(metadata.tice10m_qc_jump_threshold_K, 1);
end

function test_readPromiceAws_preserves_valid_changes_and_sensor_transitions(testCase)
   % Sub-threshold thermal evolution remains canonical even when the native
   % string changes support and begins a new thermistor-depth epoch.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   time = (0:5)';
   source = (-10:0.6:-7)';
   temperatures = [ ...
      -5.0 -7.0 -9.5
      -5.1 -7.1 -9.3
      -5.2 -7.2 NaN
      -5.3 -7.3 -8.9
      -5.4 -7.4 -8.7
      -5.5 -7.5 -8.5];
   depths = repmat([2, 7, 10.5], numel(time), 1);
   depths(4:end, 3) = 11.2;
   writeSyntheticThermistors(root, time, source, temperatures, depths);

   aws = icemodel.forcing.readPromiceAws("KAN_U", source_dir=root);
   testCase.verifyEqual(aws.tice10m_qc_flag, zeros(size(source)));
   testCase.verifyEqual(aws.tice10m, aws.tice10m_source);
   testCase.verifyEqual(aws.dtice3, depths(:, 3));
   testCase.verifyTrue(isnan(aws.tice3(3)));
end

function test_readPromiceAws_daily_skips_hourly_jump_rule(testCase)
   % The validated one-hour threshold must not be applied to daily products.
   root = string(tempname);
   mkdir(fullfile(root, "day"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   time = (0:2)';
   source = [-10; -3; -9];
   temperatures = [-9; -8.5; -8];
   writeSyntheticThermistors(root, time, source, temperatures, [], "day");

   aws = icemodel.forcing.readPromiceAws("KAN_U", source_dir=root, ...
      timescale="daily");
   testCase.verifyEqual(aws.tice10m_qc_flag, zeros(size(source)));
   testCase.verifyEqual(aws.tice10m, aws.tice10m_source);
end

function test_readPromiceAws_masks_isolated_out_of_range_target(testCase)
   % Physical-range QC must still run for a one-row source without a jump pair.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   writeSyntheticThermistors(root, 0, 2, -8, []);

   aws = icemodel.forcing.readPromiceAws("KAN_U", source_dir=root);
   testCase.verifyEqual(aws.tice10m_qc_flag, 1);
   testCase.verifyTrue(isnan(aws.tice10m));
   testCase.verifyEqual(aws.tice10m_source, ...
      icemodel.physicalConstant('Tf') + 2);
end

function test_readPromiceAws_masks_persistent_shift_to_depth_reset(testCase)
   % A large isolated sensor jump with no prompt recovery remains unreviewed
   % until the offending thermistor's next independent depth-epoch reset.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   time = (0:9)';
   source = [-10; -9.9; -9.8; -3; -3.1; -3.2; -3.3; -3.4; -3.5; -3.6];
   temperatures = [ ...
      -5.0 -7.0 -9.7
      -5.1 -7.1 -9.6
      -5.2 -7.2 -9.5
      -5.2 -7.2 -3.0
      -5.3 -7.3 -3.1
      -5.4 -7.4 -3.2
      -5.5 -7.5 -3.3
      -5.6 -7.6 -3.4
      -5.7 -7.7 -3.5
      -5.8 -7.8 -3.6];
   depths = repmat([2, 7, 10.5], numel(time), 1);
   depths(8:end, 3) = 11.2;
   writeSyntheticThermistors(root, time, source, temperatures, depths);

   [aws, metadata] = icemodel.forcing.readPromiceAws("KAN_U", ...
      source_dir=root);
   expected_flag = [0; 0; 3; 3; 3; 3; 3; 3; 0; 0];
   testCase.verifyEqual(aws.tice10m_qc_flag, expected_flag);
   testCase.verifyTrue(all(isnan(aws.tice10m(expected_flag > 0))));
   testCase.verifyEqual(metadata.tice10m_qc_persistent_sample_count, 6);
   testCase.verifyEqual(metadata.tice10m_qc_unreviewed_sample_count, 6);
end

function test_readPromiceAws_does_not_extend_shallow_sensor_jump(testCase)
   % A shallow native failure can invalidate its two target endpoints, but it
   % cannot justify masking a persistent 10 m epoch without target-depth support.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   time = (0:7)';
   source = [-10; -9.9; -3; -3.1; -3.2; -3.3; -3.4; -3.5];
   temperatures = [ ...
      -5.0 -7.0 -9.7
      -5.1 -7.1 -9.6
      -5.2 -1.0 -9.5
      -5.3 -1.1 -9.4
      -5.4 -1.2 -9.3
      -5.5 -1.3 -9.2
      -5.6 -1.4 -9.1
      -5.7 -1.5 -9.0];
   depths = repmat([2, 5, 10.5], numel(time), 1);
   writeSyntheticThermistors(root, time, source, temperatures, depths);

   [aws, metadata] = icemodel.forcing.readPromiceAws("KAN_U", ...
      source_dir=root);
   testCase.verifyEqual(aws.tice10m_qc_flag, [0; 1; 1; 0; 0; 0; 0; 0]);
   testCase.verifyEqual(metadata.tice10m_qc_persistent_sample_count, 0);
   testCase.verifyEqual(metadata.tice10m_qc_target_depth_tolerance_m, 2);
end

function test_readPromiceAws_marks_neighbor_insufficient_jump_unreviewed(testCase)
   % The target remains conservatively masked when native neighbors are absent,
   % but code 2 prevents the sparse event from masquerading as reviewed QC.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   time = (0:4)';
   source = [-10; -9.8; -4; -4.1; -4.2];
   temperatures = [-8; -8.1; -8.2; -8.3; -8.4];
   writeSyntheticThermistors(root, time, source, temperatures, []);

   [aws, metadata] = icemodel.forcing.readPromiceAws("KAN_U", ...
      source_dir=root);
   testCase.verifyEqual(aws.tice10m_qc_flag, [0; 2; 2; 0; 0]);
   testCase.verifyTrue(all(isnan(aws.tice10m(2:3))));
   testCase.verifyEqual(metadata.tice10m_qc_unreviewed_sample_count, 2);
   testCase.verifyEqual(metadata.tice10m_qc_failed_sample_count, 0);
end

function test_readPromiceAws_does_not_bridge_missing_thermistor_hours(testCase)
   % A large change across a source gap is insufficient temporal evidence for a
   % discontinuity. Sampling continuity, not row adjacency, governs the screen.
   root = string(tempname);
   mkdir(fullfile(root, "hour"))
   cleanup = onCleanup(@() rmdir(root, 's'));
   time = [0; 1; 12; 13];
   source = [-10; -9.9; -3; -3.1];
   temperatures = [-8; -8.1; -3.2; -3.3];
   writeSyntheticThermistors(root, time, source, temperatures, []);

   aws = icemodel.forcing.readPromiceAws("KAN_U", source_dir=root);
   testCase.verifyEqual(aws.tice10m_qc_flag, zeros(size(source)));
   testCase.verifyEqual(aws.tice10m, aws.tice10m_source);
end

function test_readPromiceAws_masks_known_kanu_thermistor_events(testCase)
   % The two source-level KAN_U failures that motivated the QC must remain
   % caught by the general native-neighbor rule, not a station/date exception.
   aws = icemodel.forcing.readPromiceAws("KAN_U", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2018, 8, 27), enddate=datetime(2022, 11, 3));
   events = [datetime(2018, 8, 28, 15, 0, 0, TimeZone="UTC"), ...
      datetime(2022, 9, 3, 7, 0, 0, TimeZone="UTC")];
   for event = events
      index = find(aws.Time == event, 1);
      testCase.assertNotEmpty(index);
      pair = index - 1:index;
      testCase.verifyGreaterThan(abs(diff(aws.tice10m_source(pair))), 5);
      expected_code = 1;
      if year(event) == 2022
         expected_code = 3;
      end
      testCase.verifyEqual(aws.tice10m_qc_flag(pair), ...
         expected_code * ones(2, 1));
      testCase.verifyTrue(all(isnan(aws.tice10m(pair))));
   end

   % The 2022 isolated tice8 shift remains unreviewed until the independent
   % >0.5 m native depth reset, then canonical comparisons resume.
   reset = datetime(2022, 11, 2, 5, 0, 0, TimeZone="UTC");
   unresolved = aws.Time >= events(2) - hours(1) & aws.Time <= reset;
   testCase.verifyTrue(all(aws.tice10m_qc_flag(unresolved) == 3));
   testCase.verifyTrue(all(isnan(aws.tice10m(unresolved))));
   after_reset = find(aws.Time > reset, 1);
   testCase.verifyEqual(aws.tice10m_qc_flag(after_reset), 0);
end

function test_buildPromiceData_writes_userdata(testCase)
   % A daily Data build is regularized to the hourly userdata default.

   Data = icemodel.forcing.buildPromiceData("kanm", ...
      source_dir=testCase.TestData.source_dir, frequency="daily", ...
      startdate=datetime(2016, 1, 1), enddate=datetime(2016, 12, 31));

   outdir = string(tempname);
   mkdir(outdir)
   cleanup = onCleanup(@() rmdir(outdir, 's'));

   filenames = icemodel.forcing.helpers.writeuserdata(Data, ...
      "kanm", "promice", outdir=outdir);

   testCase.verifyEqual(numel(filenames), 1);
   testCase.verifyTrue(endsWith(filenames(1), "kanm_promice_2016.mat"));
   loaded = load(filenames(1), 'Data', 'artifact_metadata');
   testCase.verifyEqual(height(loaded.Data), 24 * (height(Data) - 1) + 1);
   testCase.verifyEqual(seconds(median(diff(loaded.Data.Time))), 3600);
   testCase.verifyEqual(string(loaded.artifact_metadata.station_name), ...
      "KAN_M");
end

function test_readPromiceAws_two_boom_station(testCase)
   % A two-boom station (KAN_U) resolves through the upper-boom (_u) channels
   % and reports two booms.

   [aws, metadata] = icemodel.forcing.readPromiceAws("KAN_U", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 6, 1), enddate=datetime(2015, 6, 30));

   testCase.verifyEqual(metadata.n_booms, 2);
   testCase.verifyGreaterThan(height(aws), 24*28);
   testCase.verifyTrue(all(isfinite(aws.tair) | ~isfinite(aws.tair)));
   testCase.verifyTrue(ismember("tair", string(aws.Properties.VariableNames)));
end

function filename = writeSyntheticPromice(root, cfrac, include_raw_swd, ...
      include_swu)
   %WRITESYNTHETICPROMICE Write the minimal reader/builder NetCDF contract.
   if nargin < 3
      include_raw_swd = true;
   end
   if nargin < 4
      include_swu = true;
   end
   filename = fullfile(root, "hour", "KAN_M_hour.nc");

   % Four hourly rows exercise source-faithful raw and corrected radiation.
   nccreate(filename, 'time', 'Dimensions', {'time', 4});
   ncwrite(filename, 'time', (0:3)');
   ncwriteatt(filename, 'time', 'units', ...
      'hours since 2020-01-01 00:00:00');
   variables = ["cc", "usr", "dsr_cor", "usr_cor", "t_u", "p_u", ...
      "rh_u", "wspd_u", "dlr", "albedo"];
   values = {cfrac, [1; 2; 3; 4], [11; NaN; 33; -4], ...
      [1.1; NaN; 3.3; -0.4], [-10; -9; -8; -7], ...
      [900; 900; 900; 900], [80; 80; 80; 80], [5; 5; 5; 5], ...
      [200; 200; 200; 200], [0.8; 0.8; 0.8; 0.8]};
   % Omit both upward-shortwave source variants for derivation tests.
   if ~include_swu
      keep = ~ismember(variables, ["usr", "usr_cor"]);
      variables = variables(keep);
      values = values(keep);
   end
   if include_raw_swd
      variables(end + 1) = "dsr";
      values{end + 1} = [10; 20; 30; 40];
   end
   for k = 1:numel(variables)
      nccreate(filename, variables(k), 'Dimensions', {'time', 4});
      ncwrite(filename, variables(k), values{k});
   end

   % The reader's scalar location metadata comes from global attributes.
   ncwriteatt(filename, '/', 'site_id', 'KAN_M');
   ncwriteatt(filename, '/', 'latitude', 67.1);
   ncwriteatt(filename, '/', 'longitude', -48.5);
   ncwriteatt(filename, '/', 'altitude', 1840);
end

function filename = writeSyntheticThermistors( ...
      root, time, tice10m, temperatures, depths, subdir)
   %WRITESYNTHETICTHERMISTORS Write a minimal depth-tagged PROMICE string.
   if nargin < 6
      subdir = "hour";
   end
   filename = fullfile(root, subdir, "KAN_U_" + subdir + ".nc");
   n = numel(time);
   nccreate(filename, 'time', 'Dimensions', {'time', n});
   ncwrite(filename, 'time', time);
   units = "hours since 2020-01-01 00:00:00";
   if subdir == "day"
      units = "days since 2020-01-01 00:00:00";
   end
   ncwriteatt(filename, 'time', 'units', units);

   % t_i_10m is GEUS's derived target; the indexed string and optional depths
   % provide the native profile evidence used only for review classification.
   nccreate(filename, 't_i_10m', 'Dimensions', {'time', n});
   ncwrite(filename, 't_i_10m', tice10m);
   for k = 1:size(temperatures, 2)
      temperature_name = sprintf('t_i_%d', k);
      nccreate(filename, temperature_name, 'Dimensions', {'time', n});
      ncwrite(filename, temperature_name, temperatures(:, k));
      if ~isempty(depths)
         depth_name = sprintf('d_t_i_%d', k);
         nccreate(filename, depth_name, 'Dimensions', {'time', n});
         ncwrite(filename, depth_name, depths(:, k));
      end
   end

   % The production reader always attaches these scalar provenance attributes.
   ncwriteatt(filename, '/', 'site_id', 'KAN_U');
   ncwriteatt(filename, '/', 'latitude', 67.0);
   ncwriteatt(filename, '/', 'longitude', -47.0);
   ncwriteatt(filename, '/', 'altitude', 1840);
end
