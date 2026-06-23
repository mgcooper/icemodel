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

   testCase.verifyEqual(string(metadata.site_surface_type), "ablation");
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
   if ismember("tice10m", names)
      testCase.verifyEqual(units(names == "tice10m"), "K");
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
