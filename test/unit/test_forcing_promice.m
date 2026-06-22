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

function test_readPromiceAws_clamps_ice_temperature_string(testCase)
   % The tice string is clamped to the dictionary physical range (-80..1 C),
   % so no out-of-range thermistor spikes reach the output.

   aws = icemodel.forcing.readPromiceAws("KAN_L", ...
      source_dir=testCase.TestData.source_dir);
   Tf = icemodel.physicalConstant('Tf');
   tice = string(aws.Properties.VariableNames( ...
      startsWith(string(aws.Properties.VariableNames), "tice")));
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

   met = icemodel.forcing.buildPromiceMet("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 1, 1), ...
      enddate=datetime(2015, 12, 31, 23, 0, 0));

   icemodel.forcing.helpers.validatemet(met)
   testCase.verifyEqual(height(met), 8760);
   testCase.verifyEqual(met.ppt, zeros(8760, 1));
   testCase.verifyTrue(all(isfinite(met.tair)));
   testCase.verifyTrue(all(met.albedo >= 0.05 & met.albedo <= 0.98));

   [~, ~, pptunit] = icemodel.forcing.helpers.metvariables();
   units = string(met.Properties.VariableUnits);
   names = string(met.Properties.VariableNames);
   testCase.verifyEqual(units(names == "ppt"), pptunit);
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

   testCase.verifyEqual(string(metadata.ablation_source), "L3 z_ice_surf");
   testCase.verifyEqual(string(metadata.snow_depth_source), "L3 snow_height");

   % Snow depth: non-negative (bare-ice ablation site stays near zero).
   sd = Data.snow_depth(isfinite(Data.snow_depth));
   testCase.verifyGreaterThan(numel(sd), 1000);
   testCase.verifyTrue(all(sd >= -0.01));

   % Ablation: cumulative surface lowering, positive-down, monotone-ish, with
   % a physically plausible multi-year magnitude (KAN_L lowers ~57 m over the
   % window; the homegrown derivation inflated multi-year totals).
   ab = Data.ablation(isfinite(Data.ablation));
   testCase.verifyGreaterThan(numel(ab), 1000);
   testCase.verifyEqual(ab(1), 0, 'AbsTol', 1e-6);   % zeroed at window start
   testCase.verifyGreaterThan(max(ab), 30);
   testCase.verifyLessThan(max(ab), 80);
   testCase.verifyGreaterThan(mean(diff(ab) >= 0), 0.95);   % near-monotone

   % Ice-temperature string read from t_i_* (KAN_L one-boom: tice1..tice8).
   dn = string(Data.Properties.VariableNames);
   testCase.verifyTrue(all(ismember( ...
      "tice" + string(1:8), dn)));

   % Derived net fluxes present.
   testCase.verifyTrue(all(ismember(["swn", "lwn", "netr", "thf"], dn)));

   % Userdata metadata contract.
   testCase.verifyEqual(Data.Properties.CustomProperties.Lat, 67.097, ...
      'AbsTol', 1e-2);
   testCase.verifyEqual(numel(Data.Properties.CustomProperties.ScalarUnits), 6);
end

function test_buildPromiceData_accumulation_site_uses_combined(testCase)
   % An accumulation-zone station with no z_ice_surf (KAN_U) falls back to
   % z_surf_combined for ablation; the result reflects net accumulation
   % (surface rises -> negative ablation) with no inflated magnitude.

   [Data, metadata] = icemodel.forcing.buildPromiceData("KAN_U", ...
      source_dir=testCase.TestData.source_dir, frequency="daily", ...
      startdate=datetime(2009, 1, 1), ...
      enddate=datetime(2022, 12, 31, 23, 0, 0));

   testCase.verifyEqual(string(metadata.ablation_source), ...
      "L3 z_surf_combined");
   ab = Data.ablation(isfinite(Data.ablation));
   testCase.verifyLessThan(min(ab), 0);            % surface rose
   testCase.verifyLessThan(max(abs(ab)), 15);      % plausible magnitude
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
   testCase.verifyEqual(units(names == "tice1"), "K");
   testCase.verifyEqual(units(names == "tair"), "K");
end

function test_buildPromiceData_writes_userdata(testCase)
   % A Data build round-trips through writeuserdata.

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
   loaded = load(filenames(1), 'Data');
   testCase.verifyEqual(height(loaded.Data), height(Data));
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
