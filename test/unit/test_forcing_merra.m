function tests = test_forcing_merra
   %TEST_FORCING_MERRA Verify the MERRA-2 Data builder.
   %
   % Reads the MERRA-2 daily collection files from the raw-source
   % directory (S03 layout or local cache); skips cleanly when absent.
   %
   % Note on swd: the legacy ak4_merra artifacts derived downwelling
   % shortwave as SWGNT / (1 - SNICEALB), mixing the cell net flux with
   % the snow/ice tile albedo, which inflates swd (legacy 2013 annual
   % mean 203 W m-2). The builder uses the raw SWGDN downwelling
   % channel instead; the comparison gate is correlation-only for swd
   % and the bias is documented in the owning ExecPlan (2026-06-12).
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve sources and build the shared single-point extraction once.

   candidates = [ ...
      "/Volumes/S03/DATA/merra2/1hrly/ncfiles", ...
      string(fullfile(icemodel.getpath('data'), 'forcing', 'merra2'))];
   % Require readable MERRA files in the slv collection, not just the
   % folder (a spun-down mount passes isfolder but returns no files ->
   % skip cleanly, not error).
   hasdata = arrayfun(@(p) ~isempty(dir(fullfile(p, "slv", "*_Nx.*.nc4*"))), ...
      candidates);
   source_dir = candidates(hasdata);
   testCase.assumeTrue(~isempty(source_dir), ...
      'MERRA-2 source data not available (S03 unmounted/spun down, no cache)');
   testCase.TestData.source_dir = source_dir(1);

   [Data, metadata] = icemodel.forcing.buildMerraData( ...
      [67.1556, -49.9226], 2013, source_dir=source_dir(1));
   testCase.TestData.Data = Data;
   testCase.TestData.metadata = metadata;
end

function test_buildMerraData_shape_and_channels(testCase)
   % A one-year point build covers the full hourly year with met,
   % flux, and glacier channels plus userdata CustomProperties.

   Data = testCase.TestData.Data;
   testCase.verifyEqual(height(Data), 8760);
   testCase.verifyTrue(all(ismember( ...
      ["tair", "psfc", "swd", "lwd", "rh", "wspd", "ppt", "snowf", ...
      "runoff", "albedo", "snowd", "swe", "shf", "lhf"], ...
      string(Data.Properties.VariableNames))));
   testCase.verifyTrue(all(isfinite(Data.tair)));
   testCase.verifyGreaterThan(median(Data.rh), 50);   % percent scale
   testCase.verifyEqual(Data.Properties.CustomProperties.Lat, 67.0, ...
      'AbsTol', 0.3);
end

function test_buildMerraData_conservative_polygon_geographic(testCase)
   % Conservative polygon remap runs in MERRA's native geographic grid
   % (UseGeoCoords): a real catchment build is full-length, finite, and
   % physically near the equal-weight result (different weighting, same
   % field). MERRA's grid is coarse, so a small catchment overlaps few
   % cells; the result must still be valid.
   testCase.assumeTrue(~isempty(which('exactremap')), ...
      'exactremap toolbox not on path');
   ak4 = '/Users/mattcooper/MATLAB/projects/runoff/data/ak4/ak4.mat';
   testCase.assumeTrue(isfile(ak4), 'ak4 catchment polygon not available');
   P = load(ak4).ak4.max.poly;

   % MERRA's grid is coarse, so a small catchment overlaps few cells and the
   % exactremap infill triangulation degenerates (benign collinearity warning,
   % result still valid); silence it for clean test output.
   src = testCase.TestData.source_dir;
   wstate = warning('off', 'all');
   cleanup = onCleanup(@() warning(wstate));
   Dc = icemodel.forcing.buildMerraData(P, 2013, source_dir=src, ...
      remap="conservative");
   De = icemodel.forcing.buildMerraData(P, 2013, source_dir=src, ...
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

function test_buildMerraData_statistical_agreement_with_legacy(testCase)
   % Statistical agreement with the legacy ak4_merra artifact (nearest
   % cell vs the legacy catchment-interpolated point; rh additionally
   % carries the vapor-kernel change; swd the SWGDN decision - see the
   % file header note).

   legacy_file = ['/Users/mattcooper/MATLAB/projects/runoff/data/' ...
      'icemodel/input/userdata/ak4_merra_2013.mat'];
   testCase.assumeTrue(isfile(legacy_file), 'legacy ak4 artifact not found');

   Data = testCase.TestData.Data;
   legacy = load(legacy_file, 'Data').('Data');
   t_new = Data.Time;
   t_new.TimeZone = '';
   [tf, loc] = ismember(t_new, legacy.Time);
   testCase.verifyEqual(sum(tf), 8760);

   r = @(v) corr(Data.(v)(tf), legacy.(v)(loc(tf)), 'rows', 'complete');
   testCase.verifyGreaterThan(r("tair"), 0.999);
   testCase.verifyGreaterThan(r("psfc"), 0.999);
   testCase.verifyGreaterThan(r("runoff"), 0.99);
   testCase.verifyGreaterThan(r("ppt"), 0.98);
   testCase.verifyGreaterThan(r("albedo"), 0.99);
   testCase.verifyGreaterThan(r("swd"), 0.95);
   testCase.verifyGreaterThan(r("rh"), 0.75);

   bias = mean(Data.tair(tf) - legacy.tair(loc(tf)), 'omitnan');
   testCase.verifyLessThan(abs(bias), 1);
end

function test_buildMerraData_mass_flux_units(testCase)
   % Mass fluxes are mWE/h rates with plausible ablation-zone totals.

   Data = testCase.TestData.Data;
   testCase.verifyGreaterThan(sum(Data.runoff), 0.05);
   testCase.verifyLessThan(sum(Data.runoff), 10);
   testCase.verifyGreaterThan(sum(Data.ppt), 0.05);
   testCase.verifyLessThan(sum(Data.ppt), 5);
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

   met = icemodel.forcing.buildMerraMet([67.1556, -49.9226], 2013, ...
      source_dir=testCase.TestData.source_dir);
   icemodel.forcing.helpers.validatemet(met)
   required = icemodel.forcing.helpers.metvariables();
   testCase.verifyEqual(string(met.Properties.VariableNames(1:numel(required))), ...
      required);
   testCase.verifyEqual(height(met), 8760);
end
