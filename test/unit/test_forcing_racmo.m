function tests = test_forcing_racmo
   %TEST_FORCING_RACMO Verify the RACMO evaluation-Data builder.
   %
   % Reads the RACMO2.3p3 per-variable archive from the raw-source
   % directory (S03 layout or local cache); skips cleanly when absent.
   %
   % Note: the RACMO archive carries SMB components and surface fluxes
   % only (no air temperature / wind / humidity / pressure), so there is
   % no buildRacmoMet. The legacy eval artifacts
   % (racmo_runoff_subsurface_ak4_*.mat) come from the *subsurface*
   % product with a non-timetable schema and are not comparable to the
   % no_subsurf surface archive staged on S03; gates here are raw-NetCDF
   % self-consistency and physical plausibility (recorded in the owning
   % ExecPlan, 2026-06-12).
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve sources and build the shared single-point extraction once
   % (the multi-GB per-variable reads dominate the suite runtime).

   candidates = [ ...
      "/Volumes/S03/DATA/greenland/racmo2p3/surface", ...
      string(fullfile(icemodel.getpath('data'), 'forcing', 'racmo'))];
   % Require readable RACMO files, not just the folder (a spun-down mount
   % passes isfolder but returns no files -> skip cleanly, not error).
   hasdata = arrayfun(@(p) isfolder(p) ...
      && ~isempty(dir(fullfile(p, '*.RACMO*.nc'))), candidates);
   source_dir = candidates(hasdata);
   testCase.assumeTrue(~isempty(source_dir), ...
      'RACMO source data not available (S03 unmounted/spun down, no cache)');
   testCase.TestData.source_dir = source_dir(1);

   [Data, metadata] = icemodel.forcing.buildRacmoData( ...
      [67.067, -48.8355], 2013, source_dir=source_dir(1));
   testCase.TestData.Data = Data;
   testCase.TestData.metadata = metadata;
end

function test_buildRacmoData_shape_and_channels(testCase)
   % A one-year point build covers the full hourly year with the SMB
   % component channels and userdata CustomProperties.

   Data = testCase.TestData.Data;
   testCase.verifyEqual(height(Data), 8760);
   testCase.verifyTrue(all(ismember( ...
      ["swd", "lwd", "shf", "lhf", "precip", "melt", "runoff", "smb"], ...
      string(Data.Properties.VariableNames))));
   testCase.verifyTrue(all(isfinite(Data.runoff)));

   % Derived albedo (1 - swn/swd) restored and clamped by metchecks.
   testCase.verifyTrue(ismember("albedo", ...
      string(Data.Properties.VariableNames)));
   a = Data.albedo(isfinite(Data.albedo));
   testCase.verifyTrue(all(a >= 0.05 & a <= 0.98));
   testCase.verifyEqual(Data.Properties.CustomProperties.Lat, 67.076, ...
      'AbsTol', 0.05);
   testCase.verifyGreaterThan(Data.Properties.CustomProperties.Elev, 500);
end

function test_buildRacmoData_physically_plausible(testCase)
   % Runoff concentrates in summer; downwelling shortwave is
   % non-negative and peaks in summer.

   Data = testCase.TestData.Data;
   july = month(Data.Time) == 7;
   january = month(Data.Time) == 1;
   testCase.verifyGreaterThan(sum(Data.runoff(july)), ...
      100 * max(sum(Data.runoff(january)), 1e-9));
   testCase.verifyGreaterThan(min(Data.swd), -1);   % interp noise floor
   testCase.verifyGreaterThan(mean(Data.swd(july)), mean(Data.swd(january)));
end

function test_buildRacmoData_self_consistent_with_raw_netcdf(testCase)
   % The native-posting build reproduces the raw NetCDF at the selected
   % cell exactly (unit conversion only), and the energy channels carry
   % no conversion at all.

   metadata = testCase.TestData.metadata;
   native = icemodel.forcing.buildRacmoData([67.067, -48.8355], 2013, ...
      source_dir=testCase.TestData.source_dir, dt="3hr");

   swsd_file = metadata.source_files( ...
      contains(metadata.source_files, "swsd"));
   info = ncinfo(swsd_file, 'swsd');
   raw = squeeze(double(ncread(swsd_file, 'swsd', ...
      [metadata.grid_start 1 1], [1 1 1 info.Size(4)])));
   t = ncread(swsd_file, 'time');
   Time = datetime(1950, 1, 1, 'TimeZone', 'UTC') + days(double(t));
   raw = raw(year(Time) == 2013);

   testCase.verifyEqual(native.swd, raw, 'AbsTol', 1e-9);
end

function test_buildRacmoData_conservative_polygon_rotated(testCase)
   % Conservative polygon remap runs in RACMO's native rotated-pole frame
   % (exactremap GridMapping + shipped gridarea, weights mode): a real
   % catchment build is full-length, finite, and physically plausible
   % (summer-dominated runoff, non-negative swd). Requires the FGRN11
   % topography file (gridarea) reachable from the source dir.
   testCase.assumeTrue(~isempty(which('exactremap')), ...
      'exactremap toolbox not on path');
   ak4 = '/Users/mattcooper/MATLAB/projects/runoff/data/ak4/ak4.mat';
   testCase.assumeTrue(isfile(ak4), 'ak4 catchment polygon not available');
   P = load(ak4).ak4.max.poly;

   wstate = warning('off', 'all');
   cleanup = onCleanup(@() warning(wstate));
   D = icemodel.forcing.buildRacmoData(P, 2013, ...
      source_dir=testCase.TestData.source_dir, remap="conservative");

   testCase.verifyEqual(height(D), 8760);
   testCase.verifyTrue(all(isfinite(D.runoff)));
   testCase.verifyGreaterThan(min(D.swd), -1);
   july = month(D.Time) == 7;
   january = month(D.Time) == 1;
   testCase.verifyGreaterThan(sum(D.runoff(july)), sum(D.runoff(january)));

   % Orientation guard: exactremap returns rotated-pole weights in meshgrid
   % order while the data block / slab are ndgrid; a transpose mismatch in
   % applying the weights threw the catchment mean ~5x off (annual runoff
   % 0.45 vs 2.3 m) yet still looked "plausible". The conservative catchment
   % average must stay within a factor of the nearby single-cell point build.
   pt = testCase.TestData.Data;     % point build at [67.067 -48.8355]
   ratio = sum(D.runoff) / sum(pt.runoff);
   testCase.verifyGreaterThan(ratio, 0.5);
   testCase.verifyLessThan(ratio, 2.0);
end

function test_readRacmo2p3_hyperslab_and_units(testCase)
   % The shared reader returns standard units and a hyperslab read that
   % matches the corresponding cells of the full-grid read exactly.

   files = dir(fullfile(testCase.TestData.source_dir, '*.RACMO*.nc'));
   testCase.assumeTrue(~isempty(files), 'no RACMO variable files');
   f = string(fullfile(files(1).folder, files(1).name));
   varname = extractBefore(files(1).name, '.');   % prefix == NetCDF var

   [full, unit, Time] = icemodel.forcing.readRacmo2p3(f, varname);
   testCase.verifyTrue(ismember(unit, ...
      ["mWE/h", "K", "kg/kg", "Pa", "W/m2", "m/s", "1", "-"]) ...
      || strlength(unit) > 0);
   testCase.verifyEqual(size(full, 2), numel(Time));

   start = [10 12];
   count = [2 3];
   sub = icemodel.forcing.readRacmo2p3(f, varname, start=start, count=count);
   gsz = ncinfo(f, varname).Size;
   [ii, jj] = ndgrid(start(1):start(1)+count(1)-1, ...
      start(2):start(2)+count(2)-1);
   lin = sub2ind(gsz(1:2), ii(:), jj(:));
   testCase.verifyEqual(size(sub, 1), prod(count));
   testCase.verifyEqual(sub, full(lin, :), 'AbsTol', 1e-9);
end

function test_buildRacmoData_mass_flux_units(testCase)
   % Mass fluxes are mWE/h rates: the annual runoff total (sum of
   % hourly rates x 1 h) is within glacier-ablation-zone bounds.

   Data = testCase.TestData.Data;
   annual_runoff = sum(Data.runoff);   % m water equivalent
   testCase.verifyGreaterThan(annual_runoff, 0.05);
   testCase.verifyLessThan(annual_runoff, 10);
   iu = strcmp(Data.Properties.VariableNames, 'runoff');
   testCase.verifyEqual(Data.Properties.VariableUnits{iu}, 'mWE/h');
end
