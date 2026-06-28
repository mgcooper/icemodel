function tests = test_forcing_merra
   %TEST_FORCING_MERRA Verify the MERRA-2 Data builder.
   %
   % Reads the staged MERRA-2 daily fixture collection under data/test/forcing;
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

   source_dir = string(fullfile(icemodel.internal.fullpath('data'), 'test', ...
      'forcing', 'merra2'));
   testCase.assumeTrue(~isempty(dir(fullfile(source_dir, "slv", ...
      "*_Nx.*.nc4*"))), ...
      'MERRA-2 fixture data not available under data/test/forcing');
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
   testCase.verifyTrue(all(ismember( ...
      ["tair", "psfc", "swd", "lwd", "rh", "wspd", "ppt", "snowf", ...
      "runoff", "albedo", "snowd", "swe", "shf", "lhf"], ...
      string(Data.Properties.VariableNames))));
   testCase.verifyTrue(all(isfinite(Data.tair)));
   testCase.verifyEqual(Data.Properties.CustomProperties.Lat, 67.0, ...
      'AbsTol', 0.3);
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
   testCase.verifyTrue(all(isfinite(Dc.tair)));
   testCase.verifyLessThan(abs(mean(Dc.tair - De.tair)), 5);
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

function test_buildMerraData_calendar_from_files(testCase)
   % The calendar derives from the files present: requesting a year
   % outside the archive errors informatively instead of silently
   % assuming the legacy hardcoded 2008-2020 span.

   testCase.verifyError(@() icemodel.forcing.buildMerraData( ...
      [67.1556, -49.9226], 1999, ...
      source_dir=testCase.TestData.source_dir), ...
      'icemodel:forcing:buildMerraData:yearNotInArchive');
end

function test_buildMerraData_mass_flux_units(testCase)
   % Diagnostic mass fluxes stay mWE/h rates and remain finite over the
   % staged source window; precipitation is the canonical m s-1 rate.

   Data = testCase.TestData.Data;
   testCase.verifyTrue(all(isfinite(Data.runoff)));
   testCase.verifyTrue(all(Data.runoff >= 0));
   testCase.verifyLessThan(sum(Data.runoff), 10);

   % Precipitation is now m s-1, so the source-window depth is sum(ppt) * 3600 [m].
   annual_depth = sum(Data.ppt) * 3600;
   testCase.verifyGreaterThanOrEqual(annual_depth, 0);
   testCase.verifyLessThan(annual_depth, 5);

   iu = strcmp(Data.Properties.VariableNames, 'ppt');
   testCase.verifyEqual(Data.Properties.VariableUnits{iu}, 'm s-1');
   ir = strcmp(Data.Properties.VariableNames, 'runoff');
   testCase.verifyEqual(Data.Properties.VariableUnits{ir}, 'mWE/h');
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

   met = icemodel.forcing.buildMerraMet([67.1556, -49.9226], ...
      testCase.TestData.year, ...
      source_dir=testCase.TestData.source_dir);
   icemodel.forcing.helpers.validatemet(met)
   required = icemodel.forcing.helpers.metvariables();
   testCase.verifyEqual(string(met.Properties.VariableNames(1:numel(required))), ...
      required);
   testCase.verifyEqual(height(met), simulationHours(testCase.TestData.year));
end

function n = simulationHours(year)
   %SIMULATIONHOURS Number of hourly records in one requested model year.
   n = 24 * days(datetime(year + 1, 1, 1) - datetime(year, 1, 1));
end
