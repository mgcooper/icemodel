function tests = test_forcing_mar
   %TEST_FORCING_MAR Verify the MAR forcing/Data builders.
   %
   % Reads the staged MAR v3.11 NetCDF subset under test/data/forcing; skips
   % cleanly when absent.
   %
   % The durable gates here are (a) the met contract on a full hourly axis
   % and (b) exact self-consistency of the builder against the raw NetCDF at
   % the selected cell. The new-vs-legacy ak4 statistical comparison (the
   % legacy artifacts cannot be reproduced cell-exactly) is a user-facing
   % script, not part of this formal suite: see
   % test/interactive/ablation_comparison/compare_forcing_vs_legacy.m. Recorded in the owning
   % ExecPlan (2026-06-12).
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve raw sources; build the shared single-point extraction once
   % (the per-test fixture cost dominates otherwise).

   % Bootstrap the repo test environment so exactremap is available when
   % this file is run directly.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.bootstrap_cleanup = cleanup;

   cfg = icemodel.config('getenv', true);
   forcing_root = string(fullfile(cfg.ICEMODEL_DATA_PATH, 'forcing'));
   source_dir = fullfile(forcing_root, 'mar');
   testCase.assumeTrue(isfolder(source_dir) ...
      && ~isempty(dir(fullfile(source_dir, 'MARv3.11*.nc'))), ...
      'MAR fixture data not available under test/data/forcing');
   testCase.TestData.source_dir = source_dir;
   testCase.TestData.year = 2012;

   modis_dir = fullfile(forcing_root, 'geus', 'albedo', 'gris');
   if ~isfolder(modis_dir)
      modis_dir = "";
   end
   testCase.TestData.modis_dir = modis_dir;

   [met, metadata, Data] = icemodel.forcing.buildMarMet( ...
      [67.1556, -49.9226], ...
      testCase.TestData.year, source_dir=testCase.TestData.source_dir, ...
      modis_dir=modis_dir, dt_out="");
   testCase.TestData.met = met;
   testCase.TestData.metadata = metadata;
   testCase.TestData.Data = Data;
end

function test_buildMarMet_defaults_to_15m(testCase)
   % The public met builder advertises the repository model-met cadence while
   % the shared native fixture above explicitly opts out for raw comparisons.
   met = icemodel.forcing.buildMarMet([67.1556, -49.9226], ...
      testCase.TestData.year, source_dir=testCase.TestData.source_dir, ...
      modis_dir=testCase.TestData.modis_dir);

   testCase.verifyEqual(seconds(median(diff(met.Time))), 900);
   testCase.verifyEqual(string(met.Properties.UserData.met_resample_policy), ...
      "interval_start_zero_order_hold");
   verifyResampleMissingLowerBound(testCase, met)
end

function verifyResampleMissingLowerBound(testCase, met)
   %VERIFYRESAMPLEMISSINGLOWERBOUND Enforce source-derived gap provenance.
   expected = met.Properties.UserData.met_resample_expected_missing_counts;
   names = intersect(string(fieldnames(expected)), ...
      string(met.Properties.VariableNames), 'stable');
   for name = reshape(names, 1, [])
      testCase.verifyGreaterThanOrEqual(nnz(~isfinite(met.(char(name)))), ...
         expected.(char(name)));
   end
end

function test_buildMarMet_satisfies_met_contract(testCase)
   % A point build passes the met contract on the source file's full
   % hourly axis with ppt = snow + rain.

   met = testCase.TestData.met;
   icemodel.forcing.helpers.validatemet(met)
   testCase.verifyEqual(seconds(median(diff(met.Time))), 3600);
   expected = ncinfo(testCase.TestData.metadata.source_files(1), 'TIME').Size(1) * 24;
   testCase.verifyEqual(height(met), expected);
   testCase.verifyEqual(met.ppt, met.snowf + met.rainf, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(isfinite(met.tair)));
   testCase.verifyGreaterThan(median(met.rh), 40);   % percent scale
   testCase.verifyTrue(isfinite( ...
      met.Properties.CustomProperties.Slope));

   % The wrapper returns source Data separately and one exact finalized met
   % metadata record, including the common contract and completion policy.
   metadata = testCase.TestData.metadata;
   Data = testCase.TestData.Data;
   testCase.verifyClass(Data, 'timetable');
   testCase.verifyEqual(metadata, met.Properties.UserData);
   testCase.verifyEqual(metadata.met_variables, ...
      string(met.Properties.VariableNames(:)));
   testCase.verifyTrue(metadata.fillwithmissing);
   testCase.verifyEqual(metadata.method, Data.Properties.UserData.method);
end

function test_buildMarMet_bounds_daily_diagnostics_and_stamps_shsn2(testCase)
   % Separately processed years hold bounded endpoints, and snow depth retains
   % the explicit SHSN2-not-SHSN3 provenance even for a one-year fixture.
   met = testCase.TestData.met;
   metadata = testCase.TestData.metadata;

   testCase.verifyGreaterThanOrEqual(min(met.cfrac), 0);
   testCase.verifyLessThanOrEqual(max(met.cfrac), 1);
   testCase.verifyEqual(string(metadata.mar_snowd_source), "SHSN2");
   testCase.verifyEqual(string(metadata.mar_snowd_semantics), ...
      "snow_pack_height_above_ice");
   testCase.verifyEqual(string(metadata.mar_snowd_shsn3_policy), ...
      "not_used_total_multilayer_snow_firn_thickness");
   testCase.verifyEqual(string(metadata.mar_snowd_qc_status), ...
      "insufficient_context");
end

function test_buildMarMet_self_consistent_with_raw_netcdf(testCase)
   % The extracted unclamped channels reproduce the raw NetCDF at the
   % selected cell exactly (unit conversion only).

   met = testCase.TestData.met;
   metadata = testCase.TestData.metadata;
   filename = metadata.source_files(1);

   % tair/swd are carried through unchanged; the precipitation channel (snowf)
   % is converted from the MAR mWE/h source rate to the canonical m s-1 rate
   % inside buildMarData (/3600), so it matches the raw NetCDF only after that
   % scaling.
   pptscale = struct('tair', 1, 'swd', 1, 'snowf', 1/3600);
   for pair = {["tair", "TTH"], ["swd", "SWDH"], ["snowf", "SFH"]}
      outname = pair{1}(1);
      marname = pair{1}(2);
      raw = icemodel.forcing.readMar3p11(filename, marname, ...
         start=metadata.grid_start, count=metadata.grid_count);
      testCase.verifyEqual(met.(outname), raw(:) * pptscale.(outname), ...
         'AbsTol', 1e-9, ...
         sprintf('%s does not match the raw NetCDF', outname));
   end
end

function test_buildMarMet_uses_native_daily_mass_diagnostics(testCase)
   % Public runoff/SMB preserve source hourly structure only where its daily
   % aggregate agrees; every replaced day uses native daily RU/SMB divided 24.
   met = testCase.TestData.met;
   metadata = testCase.TestData.metadata;
   filename = metadata.source_files(1);
   source_info = ncinfo(filename);
   testCase.assumeTrue(all(ismember(["RU", "SMB"], ...
      string({source_info.Variables.Name}))), ...
      'Reduced MAR fixture omits native daily RU/SMB products');

   [runoff_daily, runoff_units, Tdaily] = ...
      icemodel.forcing.readMar3p11(filename, "RU", ...
      start=metadata.grid_start, count=metadata.grid_count, ...
      sector=metadata.mar_qc_sector);
   [smb_daily, smb_units] = icemodel.forcing.readMar3p11( ...
      filename, "SMB", start=metadata.grid_start, ...
      count=metadata.grid_count, sector=metadata.mar_qc_sector);
   expected_runoff = icemodel.forcing.helpers.dailyToHourly( ...
      runoff_daily(:) / 24, Tdaily, met.Time, method="previous");
   expected_smb = icemodel.forcing.helpers.dailyToHourly( ...
      smb_daily(:) / 24, Tdaily, met.Time, method="previous");
   runoff_hourly = icemodel.forcing.readMar3p11(filename, "RUH", ...
      start=metadata.grid_start, count=metadata.grid_count);
   smb_hourly = icemodel.forcing.readMar3p11(filename, "SMBH", ...
      start=metadata.grid_start, count=metadata.grid_count);

   testCase.verifyEqual(runoff_units, 'mWE/day');
   testCase.verifyEqual(smb_units, 'mWE/day');
   testCase.verifyEqual(reshape(sum(reshape(met.runoff, 24, []), 1), ...
      [], 1), runoff_daily(:), 'AbsTol', 1e-5);
   testCase.verifyEqual(reshape(sum(reshape(met.smb, 24, []), 1), ...
      [], 1), smb_daily(:), 'AbsTol', 1e-5);
   for item = {"runoff", runoff_hourly(:), expected_runoff; ...
         "smb", smb_hourly(:), expected_smb}'
      channel = item{1};
      raw = item{2};
      daily_rate = item{3};
      status = metadata.("mar_qc_" + channel + "_day_status");
      for day = 1:numel(status)
         rows = (day - 1) * 24 + (1:24);
         if status(day) == 1
            testCase.verifyEqual(met.(channel)(rows), raw(rows), ...
               'AbsTol', 1e-12);
         elseif status(day) == 2
            testCase.verifyEqual(met.(channel)(rows), daily_rate(rows), ...
               'AbsTol', 1e-12);
         end
      end
   end
   testCase.verifyEqual(string(metadata.mar_qc_method), ...
      "daily_constrained_hourly");
   testCase.verifyEqual(string(metadata.mar_qc_status), "applied");
   testCase.verifyEqual(metadata.mar_qc_channels, ["runoff"; "smb"]);
   testCase.verifyEqual(metadata.mar_qc_sector, 1);
   testCase.verifyEqual(string(metadata.mar_qc_sector_name), ...
      "permanent_ice");
   testCase.verifyEqual(met.Properties.UserData, metadata);

   % The paired tundra sector is independently selectable and retains the
   % same daily time shape, proving the 2-category axis is not flattened.
   tundra = icemodel.forcing.readMar3p11(filename, "RU", ...
      start=metadata.grid_start, count=metadata.grid_count, sector=2);
   testCase.verifySize(tundra, size(runoff_daily));
   testCase.verifyNotEqual(tundra, runoff_daily);
end

function test_reduced_mar_source_marks_hourly_fallback(testCase)
   % Intentionally reduced sources without RU/SMB retain RUH/SMBH but state
   % explicitly that native-daily QC was not applicable.
   metadata = testCase.TestData.metadata;
   source_info = ncinfo(metadata.source_files(1));
   has_daily = all(ismember(["RU", "SMB"], ...
      string({source_info.Variables.Name})));
   testCase.assumeFalse(has_daily, ...
      'Source contains native daily RU/SMB and does not use the fallback');

   testCase.verifyEqual(string(metadata.mar_qc_method), ...
      "daily_constrained_hourly");
   testCase.verifyEqual(string(metadata.mar_qc_status), "not_applicable");
   testCase.verifyEmpty(metadata.mar_qc_channels);
   testCase.verifyEqual(string(metadata.mar_qc_fallback), ...
      "hourly_RUH_SMBH_retained_native_daily_unavailable");
   testCase.verifyEqual(string(metadata.mar_qc_basis), ...
      "MAR native daily RU/SMB unavailable; retained hourly RUH/SMBH");
   testCase.verifyEqual(string(metadata.mar_diagnostic_status), ...
      "not_available");
   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_melt_validation_status), "not_available");
   testCase.verifyEmpty(metadata.mar_diagnostic_channels);
   testCase.verifyEmpty(metadata.mar_diagnostic_native_variables);
end

function test_full_mar_source_exposes_optional_diagnostics(testCase)
   % The production schema maps SUH directly, holds daily SU/RZ as explicitly
   % combined terms, and uses daily ME only to validate public hourly melt.
   source_dir = string(getenv("ICEMODEL_MAR_DIR"));
   if source_dir == ""
      source_dir = "/Volumes/S03/DATA/greenland/mar3p11/RUH2";
   end
   filename = fullfile(source_dir, "MARv3.11-ERA5-15km-2019.nc");
   testCase.assumeTrue(isfile(filename), ...
      'Full MAR source is unavailable for optional-diagnostic integration');

   [Data, metadata] = icemodel.forcing.buildMarData( ...
      [67.1556, -49.9226], 2019, source_dir=source_dir, fillgaps=false);
   expected_names = ["subl", "subl_evap", "refreeze_deposition"];
   testCase.verifyTrue(all(ismember(expected_names, ...
      string(Data.Properties.VariableNames))));

   [subl, subl_units] = icemodel.forcing.readMar3p11(filename, "SUH", ...
      start=metadata.grid_start, count=metadata.grid_count);
   [subl_evap, daily_units, Tdaily] = ...
      icemodel.forcing.readMar3p11(filename, "SU", ...
      start=metadata.grid_start, count=metadata.grid_count, ...
      sector=metadata.mar_diagnostic_su_sector);
   refreeze_deposition = icemodel.forcing.readMar3p11(filename, "RZ", ...
      start=metadata.grid_start, count=metadata.grid_count, sector=1);
   melt_daily = icemodel.forcing.readMar3p11(filename, "ME", ...
      start=metadata.grid_start, count=metadata.grid_count, sector=1);
   expected_subl_evap = icemodel.forcing.helpers.dailyToHourly( ...
      subl_evap(:) / 24, Tdaily, Data.Time, method="previous");
   expected_refreeze = icemodel.forcing.helpers.dailyToHourly( ...
      refreeze_deposition(:) / 24, Tdaily, Data.Time, method="previous");

   testCase.verifyEqual(subl_units, 'mWE/h');
   testCase.verifyEqual(daily_units, 'mWE/day');
   testCase.verifyEqual(Data.subl, subl(:), 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.subl_evap, expected_subl_evap, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(Data.refreeze_deposition, expected_refreeze, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(reshape(sum(reshape(Data.melt, 24, []), 1), [], 1), ...
      melt_daily(:), 'AbsTol', 1e-5);
   testCase.verifyEqual(string(metadata.mar_diagnostic_status), "applied");
   testCase.verifyEqual(string( ...
      metadata.mar_diagnostic_melt_validation_status), "validated");
   testCase.verifyEqual(metadata.mar_diagnostic_native_variables, ...
      ["SUH"; "SU"; "RZ"; "ME"]);
   testCase.verifyEqual(metadata.mar_diagnostic_su_sector, 1);
   testCase.verifyEqual(string(metadata.mar_diagnostic_su_sector_name), ...
      "permanent_ice");
end

function test_buildMarData_preserves_native_melt_gap_with_fillgaps(testCase)
   % Generic met filling must not manufacture hourly MEH support or a validated
   % ME/MEH day when the native melt source has one missing posting.
   source_dir = string(tempname);
   mkdir(source_dir)
   testCase.addTeardown(@() rmdir(source_dir, 's'));
   writeTinyMarMeltGapSource(source_dir);

   % Compare the public default against the source-faithful no-fill control.
   [filled, metadata] = icemodel.forcing.buildMarData( ...
      [70, -40], 2012, source_dir=source_dir, fillgaps=true);
   [native, native_metadata] = icemodel.forcing.buildMarData( ...
      [70, -40], 2012, source_dir=source_dir, fillgaps=false);

   % Day two retains its missing 04:00 MEH posting. Native daily ME remains a
   % finite reference, but the residual and validation result are unverified.
   missing_row = 24 + 5;
   testCase.verifyTrue(isnan(filled.melt(missing_row)));
   testCase.verifyEqual(filled.melt, native.melt);
   testCase.verifyEqual(metadata.mar_diagnostic_melt_day_status, ...
      uint8([1; 3]));
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_melt_daily_reference_mwe, [0.024; 0.024], ...
      'AbsTol', 1e-12);
   testCase.verifyTrue(isnan( ...
      metadata.mar_diagnostic_melt_residual_mwe_day(2)));
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_melt_validated_day_count, 1);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_melt_mismatch_day_count, 0);
   testCase.verifyEqual( ...
      metadata.mar_diagnostic_melt_unverified_day_count, 1);

   % Every native mass diagnostic uses the same source-preserving default path;
   % the explicit no-fill control must agree channel-for-channel and in metadata.
   for name = ["runoff", "smb", "subl", "subl_evap", ...
         "refreeze_deposition"]
      testCase.verifyEqual(filled.(name), native.(name));
   end
   fields = string(fieldnames(metadata));
   diagnostic_fields = fields(startsWith(fields, "mar_diagnostic_"));
   for field = reshape(diagnostic_fields, 1, [])
      testCase.verifyEqual(metadata.(field), native_metadata.(field));
   end
end

function test_readMar3p11_retains_large_source_mass_fluxes(testCase)
   % Values above 999 mmWE/h are legitimate source pulses; only signed 1e30
   % fill magnitudes are missing. A tiny source isolates that reader contract.
   root = tempname;
   mkdir(root)
   testCase.addTeardown(@() rmdir(root, 's'));
   filename = fullfile(root, 'MAR-test-2012.nc');
   nccreate(filename, 'RUH', 'Dimensions', ...
      {'x', 1, 'y', 1, 'ATMXH', 24, 'TIME', 1});
   values = zeros(1, 1, 24, 1);
   values(1, 1, 1, 1) = 1114;
   values(1, 1, 2, 1) = -1702;
   values(1, 1, 3, 1) = -1e34;
   ncwrite(filename, 'RUH', values);
   ncwriteatt(filename, 'RUH', 'units', 'mmWE/h');
   nccreate(filename, 'RU', 'Dimensions', ...
      {'x', 1, 'y', 1, 'SECTOR', 2, 'TIME', 1});
   ncwrite(filename, 'RU', reshape([24, 240], 1, 1, 2, 1));
   ncwriteatt(filename, 'RU', 'units', 'mmWE/day');

   % Optional MAR products use the same native conversions but retain their
   % distinct hourly/daily support and surface-sector definitions.
   nccreate(filename, 'SUH', 'Dimensions', ...
      {'x', 1, 'y', 1, 'ATMXH', 24, 'TIME', 1});
   ncwrite(filename, 'SUH', reshape(linspace(-2, 2, 24), 1, 1, 24, 1));
   ncwriteatt(filename, 'SUH', 'units', 'mmWE/h');
   nccreate(filename, 'SU', 'Dimensions', ...
      {'x', 1, 'y', 1, 'SECTOR', 2, 'TIME', 1});
   ncwrite(filename, 'SU', reshape([24, 240], 1, 1, 2, 1));
   ncwriteatt(filename, 'SU', 'units', 'mmWE/day');
   for item = {"RZ", 48; "ME", 72}'
      nccreate(filename, item{1}, 'Dimensions', ...
         {'x', 1, 'y', 1, 'SECTOR1_1', 1, 'TIME', 1});
      ncwrite(filename, item{1}, reshape(item{2}, 1, 1, 1, 1));
      ncwriteatt(filename, item{1}, 'units', 'mmWE/day');
   end

   raw = icemodel.forcing.readMar3p11(filename, "RUH");
   blocks = icemodel.forcing.readMar3p11(filename, "RU", ...
      slabs={[1 1; 1 1], [1 1; 1 1]}, sector=[1 2]);

   testCase.verifyEqual(raw(1), 1.114, 'AbsTol', 1e-12);
   testCase.verifyEqual(raw(2), -1.702, 'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(raw(3)));
   testCase.verifyEqual(blocks{1}, 0.024, 'AbsTol', 1e-12);
   testCase.verifyEqual(blocks{2}, 0.240, 'AbsTol', 1e-12);

   [subl, subl_units, Thourly] = ...
      icemodel.forcing.readMar3p11(filename, "SUH");
   [subl_evap, daily_units, Tdaily] = ...
      icemodel.forcing.readMar3p11(filename, "SU", sector=1);
   refreeze = icemodel.forcing.readMar3p11(filename, "RZ", sector=1);
   melt = icemodel.forcing.readMar3p11(filename, "ME", sector=1);
   testCase.verifyEqual(subl_units, 'mWE/h');
   testCase.verifyEqual(daily_units, 'mWE/day');
   testCase.verifyEqual(subl([1, end]), [-0.002, 0.002], ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(subl_evap / 24, 0.001, 'AbsTol', 1e-12);
   testCase.verifyEqual(numel(Thourly), 24);
   testCase.verifyEqual(numel(Tdaily), 1);
   testCase.verifyEqual(refreeze, 0.048, 'AbsTol', 1e-12);
   testCase.verifyEqual(melt, 0.072, 'AbsTol', 1e-12);
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
   Data = icemodel.forcing.buildMarData(poly, testCase.TestData.year, ...
      source_dir=testCase.TestData.source_dir, remap="equal");

   testCase.verifyGreaterThan(height(Data), 0);
   testCase.verifyTrue(all(isfinite(Data.tair)));

   % Conservative area-weighted remap (the default), via exactremap on the
   % native MAR grid with the ice mask inpainted: same shape, finite, and
   % physically near the equal-weight mean (different weighting, same field).
   Dc = icemodel.forcing.buildMarData(poly, testCase.TestData.year, ...
      source_dir=testCase.TestData.source_dir, remap="conservative");
   testCase.verifyEqual(height(Dc), height(Data));
   testCase.verifyTrue(all(isfinite(Dc.tair)));
   testCase.verifyLessThan(abs(mean(Dc.tair - Data.tair)), 5);
end

function test_buildMarMet_modis_channel(testCase)
   % The MODIS albedo channel is added when modis_dir is supplied.

   testCase.assumeTrue(testCase.TestData.modis_dir ~= "", ...
      'GEUS MODIS albedo source not available');
   met = testCase.TestData.met;
   testCase.verifyTrue(ismember("modis", ...
      string(met.Properties.VariableNames)));
   testCase.verifyGreaterThan(sum(isfinite(met.modis)), 700);
   testCase.verifyTrue(all(met.modis >= 0 & met.modis <= 1));
end

function test_data2met_orders_required_first(testCase)
   % data2met puts the required contract variables first.

   met = testCase.TestData.met;
   varnames = string(met.Properties.VariableNames);
   required = icemodel.forcing.helpers.metvariables();
   testCase.verifyEqual(varnames(1:numel(required)), required);
end

function test_buildMarMet_precip_unit_harmonized(testCase)
   % The met precipitation channel carries the canonical m s-1 rate (the MAR
   % snow + rain channels emitted in m s-1 by buildMarData), labelled in
   % VariableUnits, so MAR met agrees with the ESM-SnowMIP m s-1 convention.
   % The ppt = snow + rain identity holds.

   met = testCase.TestData.met;
   [~, ~, pptunit] = icemodel.forcing.helpers.metvariables();
   units = string(met.Properties.VariableUnits);
   names = string(met.Properties.VariableNames);
   testCase.verifyEqual(units(names == "ppt"), pptunit);
   testCase.verifyEqual(met.ppt, met.snowf + met.rainf, 'AbsTol', 1e-15);

   % Ablation-zone annual water-equivalent precipitation is small as a m s-1
   % rate but a physically plausible total when integrated over the hour:
   % sum(ppt)*3600 is the annual depth in metres (a fraction of a metre).
   annual_depth = sum(met.ppt) * 3600;
   testCase.verifyGreaterThan(annual_depth, 0);
   testCase.verifyLessThan(annual_depth, 5);
end

function test_modis_polygon_is_area_weighted_roi_mean(testCase)
   % MODIS albedo for a POLYGON (catchment) request is the area-weighted ROI
   % mean over the GEUS 5 km grid - consistent with every other gridded
   % channel - not the single nearest MODIS cell. Asserted by reproducing the
   % conservative remap directly through readGeusModis and matching the
   % builder's polygon MODIS channel; and by confirming the POINT case is
   % still the unchanged nearest cell.

   testCase.assumeTrue(testCase.TestData.modis_dir ~= "", ...
      'GEUS MODIS albedo source not available');

   metadata = testCase.TestData.metadata;
   proj = icemodel.forcing.helpers.psnProjection();
   [x0, y0] = projfwd(proj, metadata.lat, metadata.lon);
   half = 15e3;   % 15 km half-width catchment box (spans several MODIS cells)
   poly = polyshape(x0 + half*[-1 1 1 -1], y0 + half*[-1 -1 1 1]);

   % Polygon Data build with the MODIS channel (conservative default).
   Data = icemodel.forcing.buildMarData(poly, testCase.TestData.year, ...
      source_dir=testCase.TestData.source_dir, ...
      modis_dir=testCase.TestData.modis_dir);
   testCase.verifyTrue(ismember("modis", ...
      string(Data.Properties.VariableNames)));

   % Reference: the area-weighted ROI mean straight from readGeusModis on the
   % same polygon, interpolated to the Data axis the same way the builder does.
   modisfile = dir(fullfile(testCase.TestData.modis_dir, ...
      sprintf('*_%d_*.nc', testCase.TestData.year)));
   [albedo_roi, Tdaily] = icemodel.forcing.readGeusModis( ...
      string(fullfile(modisfile.folder, modisfile.name)), poly, ...
      "nearest", remap="conservative");
   ref = icemodel.forcing.helpers.dailyToHourly(albedo_roi, Tdaily, Data.Time);
   testCase.verifyEqual(Data.modis, ref, 'AbsTol', 1e-9, ...
      'polygon MODIS channel is not the area-weighted ROI mean');

   % The ROI mean must differ from the single nearest cell (otherwise the
   % polygon path collapsed to the legacy nearest-cell behaviour).
   [albedo_point, ~] = icemodel.forcing.readGeusModis( ...
      string(fullfile(modisfile.folder, modisfile.name)), ...
      [metadata.lat, metadata.lon]);
   testCase.verifyGreaterThan( ...
      max(abs(albedo_roi - albedo_point), [], 'omitnan'), 0, ...
      'ROI mean is identical to the nearest cell (remap not applied)');
end

function test_modis_point_is_nearest_cell(testCase)
   % The POINT MODIS path is unchanged: nearest cell of the GEUS grid, the
   % single-cell hyperslab the legacy readGeusModis returned.

   testCase.assumeTrue(testCase.TestData.modis_dir ~= "", ...
      'GEUS MODIS albedo source not available');

   metadata = testCase.TestData.metadata;
   modisfile = dir(fullfile(testCase.TestData.modis_dir, ...
      sprintf('*_%d_*.nc', testCase.TestData.year)));
   filename = string(fullfile(modisfile.folder, modisfile.name));

   [albedo_new, ~] = icemodel.forcing.readGeusModis(filename, ...
      [metadata.lat, metadata.lon], "nearest");

   % Reproduce the legacy nearest-cell read directly from the NetCDF.
   LON = double(ncread(filename, 'lon'));
   LAT = double(ncread(filename, 'lat'));
   proj = icemodel.forcing.helpers.psnProjection();
   [X, Y] = projfwd(proj, LAT, LON);
   [xq, yq] = projfwd(proj, metadata.lat, metadata.lon);
   [~, idx] = min(hypot(X(:) - xq, Y(:) - yq));
   [i, j] = ind2sub(size(X), idx);
   info = ncinfo(filename, 'albedo');
   ndays = info.Size(end);
   albedo_legacy = squeeze(double(ncread(filename, 'albedo', ...
      [i j 1], [1 1 ndays])));

   testCase.verifyEqual(albedo_new, albedo_legacy, 'AbsTol', 1e-12, ...
      'point MODIS path is not the nearest cell');
end

function writeTinyMarMeltGapSource(source_dir)
   %WRITETINYMARMELTGAPSOURCE Create two complete MAR days with one MEH gap.
   filename = fullfile(source_dir, 'MARv3.11-test-2012.nc');

   % The four-cell geographic/native grids keep nearest-point selection fully
   % exercised without copying the large committed MAR fixture.
   nccreate(filename, 'X', 'Dimensions', {'X', 2});
   nccreate(filename, 'Y', 'Dimensions', {'Y', 2});
   ncwrite(filename, 'X', [0; 15]);
   ncwrite(filename, 'Y', [0; 15]);
   lon = [-40.1, -39.9; -40.1, -39.9];
   lat = [69.9, 69.9; 70.1, 70.1];
   for item = { ...
         'LON', lon; 'LAT', lat; 'SH', 1000 * ones(2); ...
         'SLO', 0.01 * ones(2); 'SRF', 4 * ones(2)}'
      nccreate(filename, item{1}, 'Dimensions', {'X', 2, 'Y', 2});
      ncwrite(filename, item{1}, item{2});
   end

   % Required hourly forcing plus optional SUH use native MAR units. MEH alone
   % carries one signed fill value on the second UTC day.
   hourly = { ...
      'TTH', -10, 'C'; 'QQH', 2, 'g/kg'; 'UUH', 3, 'm/s'; ...
      'VVH', 4, 'm/s'; 'SWDH', 100, 'W/m2'; 'LWDH', 200, 'W/m2'; ...
      'ALH', 0.6, '1'; 'SFH', 0.1, 'mmWE/h'; ...
      'RFH', 0.05, 'mmWE/h'; 'MEH', 1, 'mmWE/h'; ...
      'RUH', 1, 'mmWE/h'; 'SHFH', 10, 'W/m2'; ...
      'LHFH', 5, 'W/m2'; 'SMBH', -1, 'mmWE/h'; ...
      'SUH', 0.1, 'mmWE/h'};
   for k = 1:size(hourly, 1)
      values = hourly{k, 2} * ones(2, 2, 24, 2);
      if string(hourly{k, 1}) == "MEH"
         values(:, :, 5, 2) = 1e34;
      end
      nccreate(filename, hourly{k, 1}, 'Dimensions', ...
         {'X', 2, 'Y', 2, 'ATMXH', 24, 'TIME', 2});
      ncwrite(filename, hourly{k, 1}, values);
      ncwriteatt(filename, hourly{k, 1}, 'units', hourly{k, 3});
   end

   % Daily state channels provide the standard builder support for both days.
   daily = { ...
      'SHSN2', 1, 'm'; 'CC', 0.5, '1'; ...
      'ST', -10, 'C'; 'SP', 800, 'hPa'};
   for k = 1:size(daily, 1)
      nccreate(filename, daily{k, 1}, ...
         'Dimensions', {'X', 2, 'Y', 2, 'TIME', 2});
      ncwrite(filename, daily{k, 1}, daily{k, 2} * ones(2, 2, 2));
      ncwriteatt(filename, daily{k, 1}, 'units', daily{k, 3});
   end

   % Native daily mass products make all identities source-backed. Both sectors
   % agree so the selected surface category cannot affect this preservation test.
   mass = {'RU', 24; 'SMB', -24; 'SU', 2.4; 'RZ', 1.2; 'ME', 24};
   for k = 1:size(mass, 1)
      nccreate(filename, mass{k, 1}, ...
         'Dimensions', {'X', 2, 'Y', 2, 'SECTOR', 2, 'TIME', 2});
      ncwrite(filename, mass{k, 1}, mass{k, 2} * ones(2, 2, 2, 2));
      ncwriteatt(filename, mass{k, 1}, 'units', 'mmWE/day');
   end
end
