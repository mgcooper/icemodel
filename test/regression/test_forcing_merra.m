function tests = test_forcing_merra
   %TEST_FORCING_MERRA Verify the MERRA-2 Data builder.
   %
   % Reads the staged MERRA-2 daily fixture collection under test/data/forcing;
   % skips cleanly when absent.
   %
   % Note on swd: the builder uses the raw SWGDN downwelling channel rather
   % than the legacy SWGNT/(1-SNICEALB) derivation (which mixed the cell net
   % flux with the snow/ice tile albedo and inflated swd). The new-vs-legacy
   % ak4_merra statistical comparison is a user-facing script, not part of this
   % formal suite: see test/interactive/ablation_comparison/compare_forcing_vs_legacy.m. Documented
   % in the owning ExecPlan (2026-06-12).
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve sources and build the shared single-point extraction once.

   % Bootstrap the repo test environment so exactremap is available when
   % this file is run directly.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.bootstrap_cleanup = cleanup;

   cfg = icemodel.config('getenv', true);
   source_dir = string(fullfile( ...
      cfg.ICEMODEL_DATA_PATH, 'forcing', 'merra2'));
   testCase.assumeTrue(~isempty(dir(fullfile(source_dir, "slv", ...
      "*_Nx.*.nc4*"))), ...
      'MERRA-2 fixture data not available under test/data/forcing');
   testCase.TestData.source_dir = source_dir;
   testCase.TestData.year = 2012;
   testCase.TestData.n_source_days = numel(dir(fullfile(source_dir, "slv", ...
      "*_Nx.*.nc4*")));

   [Data, metadata] = icemodel.forcing.buildMerraData( ...
      [67.1556, -49.9226], testCase.TestData.year, source_dir=source_dir);
   testCase.TestData.Data = Data;
   testCase.TestData.metadata = metadata;
end

function test_buildMerraData_shape_and_channels(testCase)
   % A point build covers the staged hourly fixture with met, flux, and glacier
   % channels plus userdata CustomProperties.

   Data = testCase.TestData.Data;
   testCase.verifyEqual(testCase.TestData.n_source_days, 31);
   testCase.verifyEqual(height(Data), simulationHours(testCase.TestData.year));
   supported = (1:24 * testCase.TestData.n_source_days)';
   testCase.verifyTrue(all(ismember( ...
      ["tair", "psfc", "swd", "lwd", "rh", "wspd", "ppt", "snowf", ...
      "runoff", "albedo", "snowd", "swe", "shf", "lhf"], ...
      string(Data.Properties.VariableNames))));
   testCase.verifyTrue(all(isfinite(Data.tair(supported))));
   testCase.verifyTrue(all(isnan(Data.tair(supported(end) + 1:end))));
   testCase.verifyEqual(Data.Properties.CustomProperties.Lat, 67.0, ...
      'AbsTol', 0.3);
end

function test_buildMerraData_metadata_matches_and_persists(testCase)
   % The builder's public metadata and payload provenance must be one exact
   % record because stageRcmForcing persists Data and discards the second output.
   Data = testCase.TestData.Data;
   metadata = testCase.TestData.metadata;
   testCase.verifyEqual(Data.Properties.UserData, metadata);

   required = ["method", "grid_start", "checks", ...
      "merra_flux_sign_convention", "merra_source_time_coordinate", ...
      "merra_tavg3_source_grid_policy"];
   testCase.verifyTrue(all(isfield(metadata, required)));

   % Exercise the real writer. It may add cadence/location facts, but it must
   % preserve every builder field identically in both supported read locations.
   outdir = string(tempname);
   icemodel.helpers.ensureDirExists(outdir);
   cleanup = onCleanup(@() rmdir(outdir, 's'));
   filenames = icemodel.forcing.helpers.writeuserdata( ...
      Data, "merra_metadata", "merra2", outdir=outdir, naming="window");
   testCase.assertNumElements(filenames, 1);
   saved = load(filenames(1), 'Data', 'artifact_metadata');
   testCase.verifyEqual(saved.Data.Properties.UserData, ...
      saved.artifact_metadata);
   fields = string(fieldnames(metadata));
   for field = reshape(fields, 1, [])
      testCase.verifyTrue(isfield(saved.artifact_metadata, field));
      testCase.verifyEqual(saved.artifact_metadata.(field), ...
         metadata.(field));
   end
end

function test_buildMerraData_conservative_polygon_geographic(testCase)
   % Conservative polygon remap runs in MERRA's native geographic grid
   % (UseGeoCoords): a catchment build is full-length for the fixture, finite, and
   % physically near the equal-weight result (different weighting, same
   % field). MERRA's grid is coarse, so a small catchment overlaps few
   % cells; the result must still be valid.
   ak4 = '/Users/mattcooper/MATLAB/projects/runoff/data/ak4/ak4.mat';
   testCase.assumeTrue(isfile(ak4), 'ak4 catchment polygon not available');
   P = load(ak4).ak4.max.poly;

   % MERRA's grid is coarse, so a small catchment overlaps few cells and the
   % exactremap infill triangulation degenerates (benign collinearity warning,
   % result still valid); silence it for clean test output.
   src = testCase.TestData.source_dir;
   wstate = warning('off', 'all');
   cleanup = onCleanup(@() warning(wstate));
   Dc = icemodel.forcing.buildMerraData(P, testCase.TestData.year, source_dir=src, ...
      remap="conservative");
   De = icemodel.forcing.buildMerraData(P, testCase.TestData.year, source_dir=src, ...
      remap="equal");

   testCase.verifyEqual(height(Dc), height(De));
   supported = 1:24 * testCase.TestData.n_source_days;
   testCase.verifyTrue(all(isfinite(Dc.tair(supported))));
   testCase.verifyTrue(all(isnan(Dc.tair(supported(end) + 1:end))));
   testCase.verifyLessThan(abs(mean(Dc.tair(supported) ...
      - De.tair(supported))), 5);
end

function test_readMerra2_hyperslab_and_units(testCase)
   % The shared reader returns standard units (mass-flux kg m-2 s-1 -> mWE/h)
   % and a hyperslab read matching the full-grid read exactly.

   root = testCase.TestData.source_dir;
   slv = dir(fullfile(root, "slv", "*_Nx.*.nc4*"));
   f = string(fullfile(slv(1).folder, slv(1).name));

   [full, unit, Time] = icemodel.forcing.readMerra2(f, 'T2M');
   testCase.verifyEqual(unit, 'K');                 % already standard
   testCase.verifyEqual(size(full, 2), numel(Time));
   testCase.verifyEqual(Time(1), datetime(testCase.TestData.year, ...
      1, 1, 0, 30, 0, TimeZone="UTC"));
   testCase.verifyEqual(diff(Time), repmat(hours(1), numel(Time) - 1, 1));

   start = [5 7];
   count = [3 2];
   sub = icemodel.forcing.readMerra2(f, 'T2M', start=start, count=count);
   gsz = ncinfo(f, 'T2M').Size;
   [ii, jj] = ndgrid(start(1):start(1)+count(1)-1, ...
      start(2):start(2)+count(2)-1);
   lin = sub2ind(gsz(1:2), ii(:), jj(:));
   testCase.verifyEqual(sub, full(lin, :), 'AbsTol', 1e-9);

   % Mass-flux conversion: PRECTOTCORR kg m-2 s-1 -> mWE/h.
   flx = dir(fullfile(root, "flx", "*_Nx.*.nc4*"));
   ff = string(fullfile(flx(1).folder, flx(1).name));
   [~, punit] = icemodel.forcing.readMerra2(ff, 'PRECTOTCORR');
   testCase.verifyEqual(punit, 'mWE/h');
end

function test_readMerra2_preserves_instantaneous_and_tavg_coordinates(testCase)
   % The source-reader boundary preserves snapshot and interval-center stamps.
   root = string(tempname);
   mkdir(root)
   cleanup = onCleanup(@() rmdir(root, 's'));

   instantaneous = fullfile(root, "inst.nc");
   tavg1 = fullfile(root, "tavg1.nc");
   tavg3 = fullfile(root, "tavg3.nc");
   writeTimeFixture(instantaneous, ...
      "minutes since 2012-01-01 00:00:00", [0, 60, 120]);
   writeTimeFixture(tavg1, ...
      "minutes since 2012-01-01 00:30:00", [0, 60, 120]);
   writeTimeFixture(tavg3, ...
      "minutes since 2012-01-01 01:30:00", [0, 180, 360]);

   [~, ~, ti] = icemodel.forcing.readMerra2(instantaneous, "TEST");
   [~, ~, t1] = icemodel.forcing.readMerra2(tavg1, "TEST");
   [~, ~, t3] = icemodel.forcing.readMerra2(tavg3, "TEST");

   origin = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC");
   testCase.verifyEqual(ti, origin + hours((0:2)'));
   testCase.verifyEqual(t1, origin + minutes(30) + hours((0:2)'));
   testCase.verifyEqual(t3, origin + minutes(90) + hours((0:3:6)'));
end

function test_readMerra2Time_rejects_non_merra_units(testCase)
   % The coordinate-only reader fails explicitly when the units cannot establish
   % the native MERRA minute origin.
   filename = string(tempname) + ".nc";
   cleanup = onCleanup(@() delete(filename));
   writeTimeFixture(filename, ...
      "hours since 2012-01-01 00:00:00", [0, 1, 2]);

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readMerra2Time(filename), ...
      'icemodel:forcing:readMerra2Time:badUnits');
end

function test_buildMerraData_holds_tavg_support_and_does_not_extrapolate(testCase)
   % tavg3 diagnostics repeat over three hours; absent fixture days stay NaN.
   Data = testCase.TestData.Data;
   supported_rows = 24 * testCase.TestData.n_source_days;
   glc = reshape(Data.albedo(1:supported_rows), 3, []);

   testCase.verifyEqual(glc(2:3, :), repmat(glc(1, :), 2, 1), ...
      'AbsTol', 1e-12);
   testCase.verifyTrue(all(isnan(Data.tair(supported_rows + 1:end))));
   testCase.verifyTrue(all(isnan(Data.albedo(supported_rows + 1:end))));
   testCase.verifyEqual(Data.Time(1), ...
      datetime(testCase.TestData.year, 1, 1, TimeZone="UTC"));
   testCase.verifyEqual( ...
      string(Data.Properties.UserData.merra_source_time_coordinate), ...
      "native_at_reader");
   testCase.verifyEqual( ...
      string(Data.Properties.UserData.merra_time_relabel_policy), ...
      "time_averaged_center_to_interval_start");
   testCase.verifyEqual( ...
      string(Data.Properties.UserData.merra_time_upsample_policy), ...
      "zero_order_hold_over_declared_support");
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
      Data, Data.Properties.UserData));
   testCase.verifyEqual( ...
      Data.Properties.UserData.merra_tavg3_source_row_count, ...
      8 * testCase.TestData.n_source_days);
   testCase.verifyGreaterThan( ...
      Data.Properties.UserData.merra_tavg3_source_time_gap_count, 0);
end

function test_buildMerraData_calendar_from_files(testCase)
   % The calendar derives from the files present: requesting a year
   % outside the archive errors informatively instead of silently
   % assuming the legacy hardcoded 2008-2020 span.

   testCase.verifyError(@() icemodel.forcing.buildMerraData( ...
      [67.1556, -49.9226], 1999, ...
      source_dir=testCase.TestData.source_dir), ...
      'icemodel:forcing:buildMerraData:yearNotInArchive');
end

function test_buildMerraData_rejects_malformed_middle_native_time(testCase)
   % Every daily coordinate is decoded: a shifted middle file cannot inherit a
   % synthetic filename-derived grid proof merely because the endpoints are valid.
   shadow = string(tempname);
   mkdir(shadow)
   cleanup = onCleanup(@() rmdir(shadow, 's'));
   for collection = ["slv", "rad", "flx", "glc"]
      copyfile(fullfile(testCase.TestData.source_dir, collection), ...
         fullfile(shadow, collection));
   end
   glc = dir(fullfile(shadow, "glc", "*_Nx.*.nc4*"));
   file = fullfile(glc(ceil(numel(glc) / 2)).folder, ...
      glc(ceil(numel(glc) / 2)).name);
   native = ncread(file, 'time');
   native(1) = native(1) + 15;
   ncwrite(file, 'time', native);

   testCase.verifyError(@() icemodel.forcing.buildMerraData( ...
      [67.1556, -49.9226], testCase.TestData.year, source_dir=shadow), ...
      'icemodel:forcing:buildMerraData:badNativeTime');
end

function test_buildMerraData_mass_flux_units(testCase)
   % Diagnostic mass fluxes stay mWE/h rates and remain finite over the
   % staged source window; precipitation is the canonical m s-1 rate.

   Data = testCase.TestData.Data;
   supported = 1:24 * testCase.TestData.n_source_days;
   testCase.verifyTrue(all(isfinite(Data.runoff(supported))));
   testCase.verifyTrue(all(Data.runoff(supported) >= 0));
   testCase.verifyLessThan(sum(Data.runoff(supported)), 10);

   % Precipitation is now m s-1, so the source-window depth is sum(ppt) * 3600 [m].
   annual_depth = sum(Data.ppt(supported)) * 3600;
   testCase.verifyGreaterThanOrEqual(annual_depth, 0);
   testCase.verifyLessThan(annual_depth, 5);

   iu = strcmp(Data.Properties.VariableNames, 'ppt');
   testCase.verifyEqual(Data.Properties.VariableUnits{iu}, 'm s-1');
   ir = strcmp(Data.Properties.VariableNames, 'runoff');
   testCase.verifyEqual(Data.Properties.VariableUnits{ir}, 'mWE/h');
end

function test_buildMerraData_flux_sign_matches_native_source(testCase)
   % Native HFLUX/EFLUX are positive upward; the builder must invert each
   % channel exactly once and mark the canonical positive-toward-surface state.
   files = dir(fullfile(testCase.TestData.source_dir, "flx", ...
      "*_Nx.*.nc4*"));
   [~, order] = sort(string({files.name}));
   filename = string(fullfile(files(order(1)).folder, files(order(1)).name));
   start = testCase.TestData.metadata.grid_start;

   native_shf = reshape(icemodel.forcing.readMerra2( ...
      filename, "HFLUX", start=start, count=[1, 1]), [], 1);
   native_lhf = reshape(icemodel.forcing.readMerra2( ...
      filename, "EFLUX", start=start, count=[1, 1]), [], 1);
   Data = testCase.TestData.Data;

   testCase.verifyEqual(Data.shf(1:numel(native_shf)), -native_shf, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(Data.lhf(1:numel(native_lhf)), -native_lhf, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual( ...
      string(Data.Properties.UserData.merra_flux_sign_convention), ...
      "positive_toward_surface");
end

function test_data2met_from_merra_data(testCase)
   % MERRA Data converts to a contract-passing met timetable directly
   % (ppt is already a channel).

   met = icemodel.forcing.data2met(testCase.TestData.Data);
   icemodel.forcing.helpers.validatemet(met)
   required = icemodel.forcing.helpers.metvariables();
   varnames = string(met.Properties.VariableNames);
   testCase.verifyEqual(varnames(1:numel(required)), required);
end

function test_data2met_precip_unit_harmonized(testCase)
   % The MERRA precipitation channel is already the canonical m s-1 rate from
   % buildMerraData, so data2met carries it through unchanged and records
   % "m s-1" in VariableUnits.

   Data = testCase.TestData.Data;
   met = icemodel.forcing.data2met(Data);

   [~, ~, pptunit] = icemodel.forcing.helpers.metvariables();
   units = string(met.Properties.VariableUnits);
   names = string(met.Properties.VariableNames);
   testCase.verifyEqual(units(names == "ppt"), pptunit);

   % The met ppt equals the Data ppt (both m s-1; no boundary conversion).
   testCase.verifyEqual(met.ppt, Data.ppt, 'AbsTol', 1e-15);
end

function test_derivable_radiation_not_stored(testCase)
   % swu/netr are canonically derivable (icemodel.processmet) and must NOT
   % be carried in the Data output; the native net fluxes swn/lwn stay.

   names = string(testCase.TestData.Data.Properties.VariableNames);
   testCase.verifyFalse(any(ismember(["swu", "netr"], names)));
   testCase.verifyTrue(all(ismember(["swn", "lwn"], names)));
end

function test_buildMerraMet_satisfies_contract(testCase)
   % buildMerraMet returns a met-contract timetable (buildMerraData +
   % data2met), with required variables ordered first.

   [met, metadata, Data] = icemodel.forcing.buildMerraMet( ...
      [67.1556, -49.9226], ...
      testCase.TestData.year, ...
      source_dir=testCase.TestData.source_dir);
   icemodel.forcing.helpers.validatemet(met)
   required = icemodel.forcing.helpers.metvariables();
   testCase.verifyEqual(string(met.Properties.VariableNames(1:numel(required))), ...
      required);
   testCase.verifyEqual(seconds(median(diff(met.Time))), 900);
   % Every native interval contributes four held quarter-hours, including the
   % final interval through its exclusive end.
   testCase.verifyEqual(height(met), ...
      4 * simulationHours(testCase.TestData.year));
   testCase.verifyEqual(string(met.Properties.UserData.met_resample_policy), ...
      "interval_start_zero_order_hold");
   expected = met.Properties.UserData.met_resample_expected_missing_counts;
   % Default staging preserves native SNICEALB gaps instead of silently
   % applying a PROMICE-specific or generic albedo fill.
   testCase.verifyGreaterThan(expected.albedo, 0);
   testCase.verifyEqual(nnz(isnan(met.albedo)), expected.albedo);
   [forcing_ready, forcing_ready_reason] = ...
      icemodel.verification.setup.metForcingReady(met);
   testCase.verifyFalse(forcing_ready);
   testCase.verifyTrue(contains(string(forcing_ready_reason), "albedo"));
   names = intersect(string(fieldnames(expected)), ...
      string(met.Properties.VariableNames), 'stable');
   for name = reshape(names, 1, [])
      testCase.verifyGreaterThanOrEqual(nnz(~isfinite(met.(char(name)))), ...
         expected.(char(name)));
   end

   % The third output exposes the source product while metadata describes the
   % delivered met product and exactly matches its UserData payload.
   testCase.verifyClass(Data, 'timetable');
   testCase.verifyEqual(metadata, met.Properties.UserData);
   testCase.verifyEqual(metadata.met_variables, ...
      string(met.Properties.VariableNames(:)));
   testCase.verifyTrue(metadata.fillwithmissing);
   testCase.verifyEqual(metadata.method, Data.Properties.UserData.method);
end

function test_buildMerraMet_blank_dt_preserves_native(testCase)
   % An explicit blank override retains the native hourly MERRA model-met axis.
   met = icemodel.forcing.buildMerraMet([67.1556, -49.9226], ...
      testCase.TestData.year, source_dir=testCase.TestData.source_dir, ...
      dt_out="");

   icemodel.forcing.helpers.validatemet(met)
   testCase.verifyEqual(seconds(median(diff(met.Time))), 3600);
   testCase.verifyEqual(height(met), simulationHours(testCase.TestData.year));
end

function test_applyMerraTimeSupport_aligned_legacy_holds_tavg3(testCase)
   % A legacy hourly window beginning on a recoverable tavg3 source row can be
   % corrected exactly without reopening the native NetCDF collection.
   hourly_time = (datetime(2012, 1, 1, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable((1:6)', RowTimes=hourly_time, ...
      VariableNames="runoff");

   proof = provenMerraTimeMetadata(hourly_time);
   [corrected, metadata, diagnostics] = ...
      icemodel.forcing.helpers.applyMerraTimeSupport(Data, proof);

   testCase.verifyEqual(corrected.runoff, [1; 1; 1; 4; 4; 4]);
   testCase.verifyGreaterThan(diagnostics.replaced_count, 0);
   testCase.verifyTrue(diagnostics.metadata_changed);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(metadata));
end

function test_applyMerraTimeSupport_unproven_aligned_legacy_fails(testCase)
   % Clock alignment alone cannot prove that an old regularizer did not invent a
   % value at a missing native 3-hour stamp.
   hourly_time = (datetime(2012, 1, 1, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable((1:6)', RowTimes=hourly_time, ...
      VariableNames="runoff");

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.applyMerraTimeSupport(Data, struct()), ...
      'icemodel:forcing:applyMerraTimeSupport:unprovenSourceGrid');
end

function test_applyMerraTimeSupport_proven_missing_source_stays_nan(testCase)
   % An inventory-proven omitted 03:00 source row invalidates its whole support
   % block even when the legacy artifact carries finite interpolated values.
   hourly_time = (datetime(2012, 1, 1, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable((1:6)', RowTimes=hourly_time, ...
      VariableNames="runoff");
   metadata = provenMerraTimeMetadata(hourly_time);
   metadata.merra_tavg3_source_row_count = 1;
   metadata.merra_tavg3_source_time_gap_count = 1;
   metadata.merra_tavg3_missing_source_times = hourly_time(4);

   corrected = icemodel.forcing.helpers.applyMerraTimeSupport(Data, metadata);

   testCase.verifyEqual(corrected.runoff(1:3), ones(3, 1));
   testCase.verifyTrue(all(isnan(corrected.runoff(4:6))));
end

function test_applyMerraTimeSupport_clipped_legacy_fails_closed(testCase)
   % A legacy window starting inside a 3-hour block lacks its leading source row;
   % the helper must error rather than replace valid samples with NaN or guess.
   hourly_time = (datetime(2012, 1, 1, 1, 0, 0, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable((1:5)', RowTimes=hourly_time, ...
      VariableNames="runoff");
   before = Data;

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.applyMerraTimeSupport(Data, struct()), ...
      'icemodel:forcing:applyMerraTimeSupport:unrecoverableLeadingBlock');
   testCase.verifyEqual(Data, before);
end

function test_applyMerraTimeSupport_canonical_clipped_is_noop(testCase)
   % Canonical provenance proves the unavailable leading support block was
   % handled upstream, so a repeat repair preserves clipped values exactly.
   hourly_time = (datetime(2012, 1, 1, 1, 0, 0, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable([10; 10; 20; 20; 20], RowTimes=hourly_time, ...
      VariableNames="runoff");
   metadata = provenMerraTimeMetadata(hourly_time);
   Data.Properties.UserData = metadata;

   [current, returned, diagnostics] = ...
      icemodel.forcing.helpers.applyMerraTimeSupport(Data, metadata);

   testCase.verifyEqual(current, Data);
   testCase.verifyEqual(returned, metadata);
   testCase.verifyEqual(diagnostics.replaced_count, 0);
   testCase.verifyFalse(diagnostics.metadata_changed);
end

function test_applyMerraTimeSupport_marked_aligned_ramp_is_repaired(testCase)
   % Canonical markers alone cannot hide a numeric ramp when its aligned native
   % tavg3 source rows are available for an exact repair.
   hourly_time = (datetime(2012, 1, 1, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable((1:6)', RowTimes=hourly_time, ...
      VariableNames="runoff");
   metadata = provenMerraTimeMetadata(hourly_time);

   [corrected, returned, diagnostics] = ...
      icemodel.forcing.helpers.applyMerraTimeSupport(Data, metadata);

   testCase.verifyEqual(corrected.runoff, [1; 1; 1; 4; 4; 4]);
   testCase.verifyEqual(returned.merra_time_replaced_value_count, ...
      diagnostics.replaced_count);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(returned));
   testCase.verifyGreaterThan(diagnostics.replaced_count, 0);
   testCase.verifyTrue(diagnostics.metadata_changed);
end

function test_applyMerraTimeSupport_marked_clipped_ramp_fails(testCase)
   % A clipped ramp with stale canonical markers lacks the leading tavg3 source
   % row, so neither the markers nor the values prove a safe reconstruction.
   hourly_time = (datetime(2012, 1, 1, 1, 0, 0, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable((1:5)', RowTimes=hourly_time, ...
      VariableNames="runoff");

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.applyMerraTimeSupport( ...
      Data, canonicalMerraTimeMetadata()), ...
      'icemodel:forcing:applyMerraTimeSupport:unrecoverableLeadingBlock');
end

function test_hasConstantMerraTavg3Support_checks_15m_and_nans(testCase)
   % Twelve 15-minute rows share one three-hour block; paired NaNs are equal but
   % a mixed finite/missing sample or a ramp invalidates the support proof.
   met_time = (datetime(2012, 1, 1, TimeZone="UTC"):minutes(15): ...
      datetime(2012, 1, 1, 2, 45, 0, TimeZone="UTC"))';
   held = timetable(ones(12, 1), RowTimes=met_time, ...
      VariableNames="runoff");
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(held));

   ramped = held;
   ramped.runoff(end) = 2;
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(ramped));
   missing = held;
   missing.runoff(:) = NaN;
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(missing));
   missing.runoff(end) = 1;
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(missing));
end

function test_hasConstantMerraTavg3Support_validates_axis_and_types(testCase)
   % The linear adjacent-row proof requires sorted UTC row times and numeric
   % glacier channels; partial blocks and tables without those channels are safe.
   partial_time = [ ...
      datetime(2012, 1, 1, 0, 15, 0, TimeZone="UTC"); ...
      datetime(2012, 1, 1, 0, 30, 0, TimeZone="UTC"); ...
      datetime(2012, 1, 1, 3, 0, 0, TimeZone="UTC"); ...
      datetime(2012, 1, 1, 3, 15, 0, TimeZone="UTC")];
   partial = timetable([1; 1; 2; 2], RowTimes=partial_time, ...
      VariableNames="runoff");
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(partial));

   % Non-UTC and nonmonotonic axes fail instead of allowing separated rows from
   % one support block to escape the adjacent comparison.
   local_time = partial;
   local_time.Properties.RowTimes.TimeZone = "";
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(local_time));
   unsorted = partial([2, 1, 3, 4], :);
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(unsorted));

   % A malformed nonnumeric physical channel fails, while an unrelated table and
   % a non-timetable input have no MERRA tavg3 values requiring this proof.
   text_channel = timetable(["a"; "a"], ...
      RowTimes=partial_time(1:2), VariableNames="runoff");
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(text_channel));
   unrelated = timetable([1; 2], RowTimes=partial_time(1:2), ...
      VariableNames="tair");
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(unrelated));
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(table(1)));
end

function test_hasConstantMerraTavg3Support_scales_over_many_blocks(testCase)
   % A full 15-minute year exercises 2,928 independent support blocks. The
   % adjacent comparison remains linear; the former full-vector scan per block
   % made this production-sized shape quadratic.
   many_times = (datetime(2012, 1, 1, TimeZone="UTC"):minutes(15): ...
      datetime(2012, 12, 31, 23, 45, 0, TimeZone="UTC"))';
   block = floor((0:numel(many_times) - 1)' / 12);
   many = timetable(block, block + 10, RowTimes=many_times, ...
      VariableNames=["runoff", "swe"]);

   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(many));
   many.swe(end - 1) = many.swe(end - 1) + 1;
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasConstantMerraTavg3Support(many));
end

function test_hasProvenMerraTavg3SourceGrid_checks_exact_inventory(testCase)
   % The source-grid predicate accepts explicit gaps but rejects stale counts,
   % non-UTC lists, out-of-grid timestamps, and nonmonotonic application axes.
   hourly_time = (datetime(2012, 1, 1, TimeZone="UTC"):hours(1): ...
      datetime(2012, 1, 1, 5, 0, 0, TimeZone="UTC"))';
   Data = timetable(ones(6, 1), RowTimes=hourly_time, ...
      VariableNames="runoff");
   metadata = provenMerraTimeMetadata(hourly_time);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(Data, metadata));

   % One declared missing source stamp remains a complete, internally consistent
   % proof for the same expected grid.
   with_gap = metadata;
   with_gap.merra_tavg3_source_row_count = 1;
   with_gap.merra_tavg3_source_time_gap_count = 1;
   with_gap.merra_tavg3_missing_source_times = hourly_time(4);
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(Data, with_gap));

   bad_count = with_gap;
   bad_count.merra_tavg3_source_row_count = 2;
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(Data, bad_count));
   outside = with_gap;
   outside.merra_tavg3_missing_source_times = hourly_time(4) + hours(1);
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(Data, outside));
   local = with_gap;
   local.merra_tavg3_missing_source_times.TimeZone = "";
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(Data, local));
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
      Data([2, 1, 3:6], :), metadata));

   % Tables without tavg3 channels need no glacier-grid proof, while malformed
   % input types and policy markers fail once such a channel is present.
   unrelated = timetable(ones(6, 1), RowTimes=hourly_time, ...
      VariableNames="tair");
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
      unrelated, struct()));
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(table(1), metadata));
   bad_policy = metadata;
   bad_policy.merra_tavg3_source_grid_policy = 'unknown';
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(Data, bad_policy));
end

function test_hasCanonicalMerraTimeSupport_rejects_malformed(testCase)
   % Builders, repair tooling, and QA share one exact policy predicate.
   metadata = canonicalMerraTimeMetadata();
   testCase.verifyTrue( ...
      icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(metadata));

   missing = rmfield(metadata, 'merra_time_relabel_policy');
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(missing));
   nonscalar = metadata;
   nonscalar.merra_time_upsample_policy = ["a", "b"];
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(nonscalar));
   bad_support = metadata;
   bad_support.merra_collection_support_hours.glc = [3, 3];
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(bad_support));
   nonfinite = metadata;
   nonfinite.merra_collection_support_hours.glc = NaN;
   testCase.verifyFalse( ...
      icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(nonfinite));
end

function n = simulationHours(year)
   %SIMULATIONHOURS Number of hourly records in one requested model year.
   n = 24 * days(datetime(year + 1, 1, 1) - datetime(year, 1, 1));
end

function metadata = canonicalMerraTimeMetadata()
   %CANONICALMERRATIMEMETADATA Return the exact shared MERRA timing contract.
   metadata = struct( ...
      'merra_source_time_coordinate', 'native_at_reader', ...
      'merra_time_relabel_policy', ...
      'time_averaged_center_to_interval_start', ...
      'merra_time_upsample_policy', ...
      'zero_order_hold_over_declared_support', ...
      'merra_collection_support_hours', ...
       struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3));
end

function metadata = provenMerraTimeMetadata(times)
   %PROVENMERRATIMEMETADATA Add an exact gap-free native glc grid proof.
   metadata = canonicalMerraTimeMetadata();
   source_times = times(minute(times) == 0 & second(times) == 0 ...
      & mod(hour(times), 3) == 0);
   missing = NaT(0, 1, 'TimeZone', 'UTC');
   metadata.merra_tavg3_source_grid_policy = ...
      'native_glc_timestamp_inventory';
   metadata.merra_tavg3_expected_source_row_count = numel(source_times);
   metadata.merra_tavg3_source_row_count = numel(source_times);
   metadata.merra_tavg3_source_time_gap_count = 0;
   metadata.merra_tavg3_missing_source_times = missing;
end

function writeTimeFixture(filename, time_units, time_values)
   %WRITETIMEFIXTURE Create one-cell NetCDF with an explicit native time axis.
   ntime = numel(time_values);
   nccreate(filename, 'TEST', ...
      'Dimensions', {'lon', 1, 'lat', 1, 'time', ntime});
   ncwrite(filename, 'TEST', reshape(1:ntime, 1, 1, []));
   ncwriteatt(filename, 'TEST', 'units', 'K');
   nccreate(filename, 'time', 'Dimensions', {'time', ntime});
   ncwrite(filename, 'time', time_values);
   ncwriteatt(filename, 'time', 'units', time_units);
end
