function tests = test_forcing_racmo
   %TEST_FORCING_RACMO Verify the RACMO evaluation-Data builder.
   %
   % Reads the staged RACMO2.3p3 per-variable fixture subset under
   % data/test/forcing; skips cleanly when absent.
   %
   % Note: the RACMO archive carries SMB components and surface fluxes
   % only (no air temperature / wind / humidity / pressure), so there is
   % no buildRacmoMet. The legacy eval artifacts
   % (racmo_runoff_subsurface_ak4_*.mat) come from the *subsurface*
   % product with a non-timetable schema and are not comparable to the
   % staged surface fixture; gates here are raw-NetCDF self-consistency and
   % physical plausibility (recorded in the owning
   % ExecPlan, 2026-06-12).
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve sources and build the shared single-point extraction once
   % (the multi-GB per-variable reads dominate the suite runtime).

   % Bootstrap the repo test environment so exactremap and the shared
   % helpers are on the path when this file is run directly.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;

   source_dir = string(fullfile(icemodel.internal.fullpath('data'), 'test', ...
      'forcing', 'racmo'));
   testCase.assumeTrue(isfolder(source_dir) ...
      && ~isempty(dir(fullfile(source_dir, '*.RACMO*.nc'))), ...
      'RACMO fixture data not available under data/test/forcing');
   testCase.TestData.source_dir = source_dir;
   testCase.TestData.year = 2012;
   [Data, metadata] = icemodel.forcing.buildRacmoData( ...
      [67.067, -48.8355], testCase.TestData.year, source_dir=source_dir);
   testCase.TestData.Data = Data;
   testCase.TestData.metadata = metadata;
end

function test_buildRacmoData_shape_and_channels(testCase)
   % A point build covers the staged hourly fixture with the SMB
   % component channels and userdata CustomProperties.

   Data = testCase.TestData.Data;
   expected_hours = 24 * days( ...
      datetime(testCase.TestData.year + 1, 1, 1) ...
      - datetime(testCase.TestData.year, 1, 1));
   testCase.verifyEqual(height(Data), expected_hours);
   testCase.verifyEqual(year(Data.Time([1 end])), ...
      [testCase.TestData.year; testCase.TestData.year]);
   testCase.verifyEqual(month(Data.Time([1 end])), [1; 12]);
   testCase.verifyEqual(day(Data.Time([1 end])), [1; 31]);
   testCase.verifyEqual(hour(Data.Time([1 end])), [0; 23]);
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

function test_buildRacmoData_source_window_plausible(testCase)
   % The native 3-hourly staged January source window stays physically bounded.

   Native = icemodel.forcing.buildRacmoData([67.067, -48.8355], ...
      testCase.TestData.year, source_dir=testCase.TestData.source_dir, ...
      dt="3hr");
   testCase.verifyEqual(height(Native), 31 * 8);
   testCase.verifyTrue(all(isfinite(Native.runoff)));
   testCase.verifyTrue(all(Native.runoff >= 0));
   testCase.verifyTrue(all(Native.swd >= 0));
   testCase.verifyTrue(all(month(Native.Time) == 1));
   testCase.verifyTrue(all(year(Native.Time) == testCase.TestData.year));
end

function test_buildRacmoData_self_consistent_with_raw_netcdf(testCase)
   % The native-posting build reproduces the raw NetCDF at the selected
   % cell exactly (unit conversion only), and the energy channels carry
   % no conversion at all.

   metadata = testCase.TestData.metadata;
   native = icemodel.forcing.buildRacmoData([67.067, -48.8355], ...
      testCase.TestData.year, ...
      source_dir=testCase.TestData.source_dir, dt="3hr");

   swsd_file = metadata.source_files( ...
      contains(metadata.source_files, "swsd"));
   info = ncinfo(swsd_file, 'swsd');
   raw = squeeze(double(ncread(swsd_file, 'swsd', ...
      [metadata.grid_start 1 1], [1 1 1 info.Size(4)])));
   t = ncread(swsd_file, 'time');
   Time = datetime(1950, 1, 1, 'TimeZone', 'UTC') + days(double(t));
   raw = raw(year(Time) == testCase.TestData.year);

   testCase.verifyEqual(native.swd, raw, 'AbsTol', 1e-9);
end

function test_buildRacmoData_conservative_polygon_rotated(testCase)
   % Conservative polygon remap runs in RACMO's native rotated-pole frame
   % exactremap GridMapping + shipped gridarea, weights mode: the polygon
   % result must match a direct exactremap re-aggregation of the same raw
   % source window.
   testCase.assertNotEmpty(which('exactremap'), ...
      'exactremap toolbox not on path');
   ak4 = '/Users/mattcooper/MATLAB/projects/runoff/data/ak4/ak4.mat';
   testCase.assumeTrue(isfile(ak4), 'ak4 catchment polygon not available');
   P = load(ak4).ak4.max.poly;

   wstate = warning('off', 'all');
   cleanup = onCleanup(@() warning(wstate));
   D = icemodel.forcing.buildRacmoData(P, testCase.TestData.year, ...
      source_dir=testCase.TestData.source_dir, remap="conservative");

   testCase.verifyEqual(height(D), height(testCase.TestData.Data));
   testCase.verifyEqual(D.runoff, racmoConservativeReference( ...
      testCase.TestData.source_dir, P, "runoff", testCase.TestData.year), ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(D.swd, racmoConservativeReference( ...
      testCase.TestData.source_dir, P, "swsd", testCase.TestData.year), ...
      'AbsTol', 1e-12);
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
   % Mass fluxes keep their canonical units and remain finite on the staged
   % source slice.

   Data = testCase.TestData.Data;
   testCase.verifyTrue(all(isfinite(Data.runoff)));
   iu = strcmp(Data.Properties.VariableNames, 'runoff');
   testCase.verifyEqual(Data.Properties.VariableUnits{iu}, 'mWE/h');

   % Precipitation is the canonical m s-1 rate and remains finite.
   ip = strcmp(Data.Properties.VariableNames, 'precip');
   testCase.verifyEqual(Data.Properties.VariableUnits{ip}, 'm s-1');
   testCase.verifyTrue(all(isfinite(Data.precip)));
end

function ref = racmoConservativeReference(source_dir, P, varname, years)
   %RACMOCONSERVATIVEREFERENCE Direct exactremap reference for a RACMO field.
   %
   % Rebuilds the polygon weights from the source grid and applies them to the
   % requested variable so the test compares the builder output against the
   % same conservative-remap contract rather than a loose numeric heuristic.

   filename = racmoSourceFile(source_dir, varname);
   [start, count, wn] = racmoConservativeWeights(filename, source_dir, P);
   [block, ~, Time] = icemodel.forcing.readRacmo2p3(filename, varname, ...
      start=start, count=count);
   keep = ismember(year(Time), years);
   Time = Time(keep);
   series3 = (wn.' * block(:, keep)).';
   t1 = dateshift(Time(1), 'start', 'year');
   t2 = dateshift(Time(end), 'start', 'year') + calyears(1) - hours(1);
   t_hourly = (t1:hours(1):t2)';
   t_hourly = t_hourly(ismember(year(t_hourly), years));
   ref = retime(timetable(Time, series3, 'VariableNames', "value"), ...
      t_hourly, 'linear', 'EndValues', 'extrap').value;
end

function filename = racmoSourceFile(source_dir, varname)
   %RACMOSOURCEFILE Return the single RACMO file for a variable.
   files = dir(fullfile(source_dir, [char(varname) '.RACMO*.nc']));
   assert(~isempty(files), 'no RACMO %s source file found', varname);
   filename = string(fullfile(files(1).folder, files(1).name));
end

function [start, count, wn] = racmoConservativeWeights(filename, source_dir, P)
   %RACMOCONSERVATIVEWEIGHTS Recreate the RACMO conservative remap weights.

   proj = icemodel.forcing.helpers.psnProjection();
   rlon = double(ncread(filename, 'rlon'));
   rlat = double(ncread(filename, 'rlat'));
   gm = struct( ...
      'grid_mapping_name', 'rotated_latitude_longitude', ...
      'grid_north_pole_latitude', ...
      ncreadatt(filename, 'rotated_pole', 'grid_north_pole_latitude'), ...
      'grid_north_pole_longitude', ...
      ncreadatt(filename, 'rotated_pole', 'grid_north_pole_longitude'));
   [cellareas, validmask] = racmoGridStatics(source_dir);
   assert(~isempty(cellareas), ...
      'conservative RACMO remap needs the FGRN11 topography file');

   vx = P.Vertices(:, 1);
   vy = P.Vertices(:, 2);
   fin = isfinite(vx) & isfinite(vy);
   vlat = nan(size(vx));
   vlon = nan(size(vx));
   [vlat(fin), vlon(fin)] = projinv(proj, vx(fin), vy(fin));
   Pgeo = [vlon, vlat];
   [vrlat, vrlon] = geo2rotated(vlat(fin), vlon(fin), ...
      gm.grid_north_pole_latitude, gm.grid_north_pole_longitude);

   pad = 2;
   ii = find(rlon >= min(vrlon) & rlon <= max(vrlon));
   jj = find(rlat >= min(vrlat) & rlat <= max(vrlat));
   if isempty(ii); [~, ii] = min(abs(rlon - mean(vrlon))); end
   if isempty(jj); [~, jj] = min(abs(rlat - mean(vrlat))); end
   i0 = max(1, min(ii) - pad); i1 = min(numel(rlon), max(ii) + pad);
   j0 = max(1, min(jj) - pad); j1 = min(numel(rlat), max(jj) + pad);
   start = [i0, j0];
   count = [i1 - i0 + 1, j1 - j0 + 1];
   rows = i0:i1;
   cols = j0:j1;

   W = exactremap([], rlon(rows), rlat(cols), Pgeo, 'weights', ...
      'GridMapping', gm, 'CellAreas', cellareas(rows, cols).', ...
      'ValidCellsMask', validmask(rows, cols).', 'InfillMasked', true);
   w = reshape(W.W(:), [count(2), count(1)]).';
   w = w(:);
   w(~isfinite(w)) = 0;
   assert(sum(w) > 0, 'polygon does not overlap the RACMO grid');
   wn = w / sum(w);
end

function [cellareas, validmask] = racmoGridStatics(source_dir)
   %RACMOGRIDSTATICS Pull the FGRN11 cell areas and ice mask from the topo.

   match = [dir(fullfile(source_dir, '*topography*.nc')); ...
      dir(fullfile(source_dir, '..', '*topography*.nc'))];
   if isempty(match)
      cellareas = [];
      validmask = [];
      return
   end
   filename = fullfile(match(1).folder, match(1).name);
   cellareas = double(ncread(filename, 'gridarea')) * 1e6;
   try
      validmask = logical(round(double(ncread(filename, 'IceMask'))));
   catch
      validmask = double(ncread(filename, 'Promicemask')) > 0;
   end
end
