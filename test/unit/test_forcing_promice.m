function tests = test_forcing_promice
   %TEST_FORCING_PROMICE Verify the PROMICE forcing/evaluation builders.
   %
   % These tests read the PROMICE v3 NetCDF bundle from the raw-source
   % directory (the S03 external drive layout, or a local cache) and skip
   % cleanly when it is not available.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve the raw PROMICE source; skip the whole file when absent.

   candidates = [ ...
      "/Volumes/S03/DATA/greenland/geus/aws/v3", ...
      string(fullfile(icemodel.getpath('data'), 'forcing', 'promice'))];
   % Require readable PROMICE NetCDF files, not just the folder (a
   % spun-down mount passes isfolder but returns no files -> skip cleanly,
   % not error).
   hasdata = arrayfun(@(p) ...
      ~isempty(dir(fullfile(p, "hourly", "*_v03.nc"))), candidates);
   source_dir = candidates(hasdata);
   testCase.assumeTrue(~isempty(source_dir), ...
      'PROMICE v3 source data not available (S03 unmounted/spun down, no cache)');
   testCase.TestData.source_dir = source_dir(1);
end

function test_readPromiceAws_standard_channels(testCase)
   % The reader returns icemodel-named channels on a UTC hourly axis
   % with station metadata.

   [aws, metadata] = icemodel.forcing.readPromiceAws("kanm", ...
      source_dir=testCase.TestData.source_dir);

   needed = ["tair", "swd", "lwd", "albedo", "wspd", "rh", "psfc", ...
      "tsfc", "boom_height", "transducer_depth", "tice1"];
   testCase.verifyTrue(all(ismember(needed, ...
      string(aws.Properties.VariableNames))));
   testCase.verifyEqual(aws.Time.TimeZone, 'UTC');
   testCase.verifyEqual(minute(aws.Time(1)), 0);   % snapped to hour start
   testCase.verifyEqual(metadata.lat, 67.067, 'AbsTol', 1e-3);
   testCase.verifyGreaterThan(median(aws.rh, 'omitnan'), 50);   % percent
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

function test_buildPromiceMet_satisfies_met_contract(testCase)
   % A full-year KAN_M met build passes the met contract with a complete
   % hourly axis.

   met = icemodel.forcing.buildPromiceMet("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 1, 1), ...
      enddate=datetime(2015, 12, 31, 23, 0, 0));

   icemodel.forcing.helpers.validatemet(met)
   testCase.verifyEqual(height(met), 8760);
   testCase.verifyEqual(met.ppt, zeros(8760, 1));
   testCase.verifyTrue(all(isfinite(met.tair)));
   testCase.verifyTrue(all(met.albedo >= 0.05 & met.albedo <= 0.98));
end

function test_buildPromiceMet_matches_legacy_kanm_2015(testCase)
   % Tolerance comparison against the staged legacy artifact
   % demo/data/input/met/met_kanm_kanm_2015_1hr.mat. Bounds are
   % evidence-based (measured 2026-06-12, recorded in the owning
   % ExecPlan): most channels agree to rounding; swd carries upstream
   % source-version drift (GEUS reprocessed the tilt-corrected
   % radiation between the legacy .txt vintage and the v3 .nc), so its
   % gate is mean-deviation plus a generous max.

   legacy_file = fullfile(icemodel.internal.fullpath('demo'), 'data', ...
      'input', 'met', 'met_kanm_kanm_2015_1hr.mat');
   testCase.assumeTrue(isfile(legacy_file), 'legacy kanm artifact not staged');

   met = icemodel.forcing.buildPromiceMet("kanm", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 1, 1), ...
      enddate=datetime(2015, 12, 31, 23, 0, 0));
   legacy = load(legacy_file, 'met').('met');

   t_new = met.Time;
   t_new.TimeZone = '';
   [tf, loc] = ismember(t_new, legacy.Time);
   testCase.verifyEqual(sum(tf), 8760);

   maxdev = @(v) max(abs(met.(v)(tf) - legacy.(v)(loc(tf))), [], 'omitnan');
   meandev = @(v) mean(abs(met.(v)(tf) - legacy.(v)(loc(tf))), 'omitnan');

   testCase.verifyLessThan(maxdev("tair"), 1e-3);
   testCase.verifyLessThan(maxdev("psfc"), 2);
   testCase.verifyLessThan(maxdev("wspd"), 0.05);
   testCase.verifyLessThan(maxdev("rh"), 0.5);
   testCase.verifyLessThan(maxdev("lwd"), 2);
   testCase.verifyLessThan(maxdev("albedo"), 0.1);
   testCase.verifyLessThan(maxdev("tsfc"), 2);
   testCase.verifyLessThan(maxdev("swd"), 70);
   testCase.verifyLessThan(meandev("swd"), 1);
end

function test_buildPromiceData_evaluation_channels(testCase)
   % The Data builder produces despiked ablation, non-negative snow
   % depth, derived net fluxes, and the userdata CustomProperties.

   Data = icemodel.forcing.buildPromiceData("kanm", ...
      source_dir=testCase.TestData.source_dir, frequency="daily");

   % Ablation: trimmed before the curated 2012-05-04 record start.
   pre = Data.Time < datetime(2012, 5, 4, 'TimeZone', 'UTC');
   testCase.verifyTrue(all(~isfinite(Data.ablation(pre))));
   testCase.verifyGreaterThan(sum(isfinite(Data.ablation)), 1000);

   % Snow depth: non-negative where defined, with a plausible maximum.
   sd = Data.snow_depth(isfinite(Data.snow_depth));
   testCase.verifyTrue(all(sd >= 0));
   testCase.verifyLessThan(max(sd), 3);

   % Derived net fluxes present; observational gaps preserved (daily
   % means smooth partial-day gaps, but whole missing days stay NaN).
   testCase.verifyTrue(all(ismember(["swn", "lwn", "netr", "thf"], ...
      string(Data.Properties.VariableNames))));
   testCase.verifyGreaterThan(sum(~isfinite(Data.transducer_depth)), 0);

   % Userdata metadata contract.
   testCase.verifyEqual(Data.Properties.CustomProperties.Lat, 67.067, ...
      'AbsTol', 1e-3);
   testCase.verifyEqual(numel(Data.Properties.CustomProperties.ScalarUnits), 6);
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

function test_readPromiceAws_second_station_smoke(testCase)
   % The reader generalizes beyond the KAN sites (any v3 station file).

   [aws, metadata] = icemodel.forcing.readPromiceAws("nukl", ...
      source_dir=testCase.TestData.source_dir, ...
      startdate=datetime(2015, 6, 1), enddate=datetime(2015, 6, 30));

   testCase.verifyGreaterThan(height(aws), 24*28);
   testCase.verifyTrue(contains(lower(metadata.source_file), "nuk_l"));
   testCase.verifyTrue(all(isfinite(aws.tair)));
end
