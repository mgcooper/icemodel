function tests = test_forcing_mar
   %TEST_FORCING_MAR Verify the MAR forcing/Data builders.
   %
   % Reads MAR v3.11 yearly NetCDF from the raw-source directory (S03
   % external drive layout or local cache); skips cleanly when absent.
   %
   % Note on the legacy comparison: the legacy ak4 artifacts cannot be
   % reproduced cell-exactly. Their cell selection came from ncrowcol's
   % independent row/column nearest-match on a curvilinear grid (the
   % stored Lat/Lon metadata are mutually inconsistent), and no single
   % cell of the current MAR archive reproduces the stored series
   % (best whole-year match deviates 6.4 K in tair). The gates here are
   % therefore (a) exact self-consistency of the builder against the
   % raw NetCDF at the selected cell, and (b) statistical agreement
   % with the legacy artifact at the one-cell-offset scale. Recorded in
   % the owning ExecPlan (2026-06-12).
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve raw sources; build the shared single-point extraction once
   % (the per-test fixture cost dominates otherwise).

   candidates = [ ...
      "/Volumes/S03/DATA/greenland/mar3p11/RUH2", ...
      string(fullfile(icemodel.getpath('data'), 'forcing', 'mar'))];
   % Require actual MAR files to be readable, not just the folder to exist:
   % an external drive that is mounted but spun down passes isfolder yet
   % returns no files, which must skip cleanly rather than error.
   hasdata = arrayfun(@(p) isfolder(p) ...
      && ~isempty(dir(fullfile(p, 'MARv3.11*.nc'))), candidates);
   source_dir = candidates(hasdata);
   testCase.assumeTrue(~isempty(source_dir), ...
      'MAR source data not available (S03 unmounted/spun down, no cache)');
   testCase.TestData.source_dir = source_dir(1);

   modis_dir = "/Volumes/S03/DATA/greenland/geus/albedo/gris";
   if ~isfolder(modis_dir)
      modis_dir = "";
   end
   testCase.TestData.modis_dir = modis_dir;

   [met, metadata] = icemodel.forcing.buildMarMet([67.1556, -49.9226], ...
      2009, source_dir=testCase.TestData.source_dir, ...
      modis_dir=modis_dir);
   testCase.TestData.met = met;
   testCase.TestData.metadata = metadata;
end

function test_buildMarMet_satisfies_met_contract(testCase)
   % A full-year point build passes the met contract on a complete
   % hourly axis with ppt = snow + rain.

   met = testCase.TestData.met;
   icemodel.forcing.helpers.validatemet(met)
   testCase.verifyEqual(height(met), 8760);
   testCase.verifyEqual(met.ppt, met.snow + met.rain, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(isfinite(met.tair)));
   testCase.verifyGreaterThan(median(met.rh), 40);   % percent scale
end

function test_buildMarMet_self_consistent_with_raw_netcdf(testCase)
   % The extracted unclamped channels reproduce the raw NetCDF at the
   % selected cell exactly (unit conversion only).

   met = testCase.TestData.met;
   metadata = testCase.TestData.metadata;
   filename = metadata.source_files(1);

   for pair = {["tair", "TTH"], ["swd", "SWDH"], ["snow", "SFH"]}
      outname = pair{1}(1);
      marname = pair{1}(2);
      raw = icemodel.forcing.readMar3p11(filename, marname, ...
         start=metadata.grid_start, count=metadata.grid_count);
      testCase.verifyEqual(met.(outname), raw(:), 'AbsTol', 1e-9, ...
         sprintf('%s does not match the raw NetCDF', outname));
   end
end

function test_buildMarData_polygon_average(testCase)
   % A small polygon around the point averages multiple cells and
   % produces finite channels.

   metadata = testCase.TestData.metadata;
   proj = icemodel.forcing.helpers.psnProjection();
   [x0, y0] = projfwd(proj, metadata.lat, metadata.lon);
   half = 20e3;   % 20 km half-width box around the cell
   poly = polyshape(x0 + half*[-1 1 1 -1], y0 + half*[-1 -1 1 1]);

   % Equal-weight aggregation (no exactremap dependency).
   Data = icemodel.forcing.buildMarData(poly, 2009, ...
      source_dir=testCase.TestData.source_dir, remap="equal");

   testCase.verifyGreaterThan(height(Data), 0);
   testCase.verifyTrue(all(isfinite(Data.tair)));

   % Conservative area-weighted remap (the default), via exactremap on the
   % native MAR grid with the ice mask inpainted: same shape, finite, and
   % physically near the equal-weight mean (different weighting, same field).
   % Skipped only when the exactremap toolbox is absent.
   testCase.assumeTrue(~isempty(which('exactremap')), ...
      'exactremap toolbox not on path');
   Dc = icemodel.forcing.buildMarData(poly, 2009, ...
      source_dir=testCase.TestData.source_dir, remap="conservative");
   testCase.verifyEqual(height(Dc), height(Data));
   testCase.verifyTrue(all(isfinite(Dc.tair)));
   testCase.verifyLessThan(abs(mean(Dc.tair - Data.tair)), 5);
end

function test_buildMarMet_statistical_agreement_with_legacy(testCase)
   % Statistical agreement with the legacy ak4 artifact (one-cell
   % offset scale; see file header note).

   legacy_file = ['/Users/mattcooper/MATLAB/projects/runoff/data/' ...
      'icemodel/input/met/met_ak4_mar_2009_1hr.mat'];
   testCase.assumeTrue(isfile(legacy_file), 'legacy ak4 artifact not found');

   met = testCase.TestData.met;
   legacy = load(legacy_file, 'met').('met');
   t_new = met.Time;
   t_new.TimeZone = '';
   [tf, loc] = ismember(t_new, legacy.Time);
   testCase.verifyEqual(sum(tf), 8760);

   r = @(v) corr(met.(v)(tf), legacy.(v)(loc(tf)), 'rows', 'complete');
   bias = @(v) mean(met.(v)(tf) - legacy.(v)(loc(tf)), 'omitnan');

   testCase.verifyGreaterThan(r("tair"), 0.98);
   testCase.verifyGreaterThan(r("lwd"), 0.95);
   testCase.verifyGreaterThan(r("psfc"), 0.95);
   testCase.verifyLessThan(abs(bias("tair")), 2);
   testCase.verifyLessThan(abs(bias("swd")), 5);
end

function test_buildMarMet_modis_channel(testCase)
   % The MODIS albedo channel is added when modis_dir is supplied.

   testCase.assumeTrue(testCase.TestData.modis_dir ~= "", ...
      'GEUS MODIS albedo source not available');
   met = testCase.TestData.met;
   testCase.verifyTrue(ismember("modis", ...
      string(met.Properties.VariableNames)));
   testCase.verifyGreaterThan(sum(isfinite(met.modis)), 8000);
   testCase.verifyTrue(all(met.modis >= 0 & met.modis <= 1));
end

function test_data2met_orders_required_first(testCase)
   % data2met puts the required contract variables first.

   met = testCase.TestData.met;
   varnames = string(met.Properties.VariableNames);
   required = icemodel.forcing.helpers.metvariables();
   testCase.verifyEqual(varnames(1:numel(required)), required);
end
