function tests = test_forcing_helpers
   %TEST_FORCING_HELPERS Verify the icemodel.forcing.helpers contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Per-test temporary output directory for the write helpers.

   testCase.TestData.outdir = string(tempname);
   mkdir(testCase.TestData.outdir)
end

function teardown(testCase)
   % Remove the temporary output directory.

   if isfolder(testCase.TestData.outdir)
      rmdir(testCase.TestData.outdir, 's')
   end
end

%% metfilename

function test_metfilename_window_form(testCase)
   % Window form encodes YYYYMMDD start/end stamps and the dt suffix.

   name = icemodel.forcing.helpers.metfilename("kanm", "mar", ...
      datetime(2015, 10, 1), datetime(2016, 9, 30), 3600);
   testCase.verifyEqual(name, "met_kanm_mar_20151001_20160930_1hr.mat");
end

function test_metfilename_yearly_form(testCase)
   % Per-year legacy form encodes the year and supports both dt inputs.

   name = icemodel.forcing.helpers.metfilename("kanm", "mar", 2016, [], 900);
   testCase.verifyEqual(name, "met_kanm_mar_2016_15m.mat");

   name = icemodel.forcing.helpers.metfilename("kanm", "mar", 2016, [], "1hr");
   testCase.verifyEqual(name, "met_kanm_mar_2016_1hr.mat");
end

function test_metfilename_rejects_bad_dt(testCase)
   % Unsupported timesteps must error, mirroring createMetFileNames.

   testCase.verifyError(@() icemodel.forcing.helpers.metfilename( ...
      "kanm", "mar", 2016, [], 1800), ...
      'icemodel:forcing:metfilename:unsupportedTimestep');
end

function test_metfilename_roundtrip_with_createMetFileNames(testCase)
   % The writer-side names must match what the reader-side
   % icemodel.createMetFileNames builds from equivalent opts.

   opts = struct('sitename', 'kanm', 'forcings', 'mar', ...
      'simyears', 2016, 'dt', 3600);
   expected = icemodel.createMetFileNames(opts);
   actual = icemodel.forcing.helpers.metfilename("kanm", "mar", 2016, [], 3600);
   testCase.verifyEqual(char(actual), expected{1});

   opts.startdate = datetime(2015, 10, 1);
   opts.enddate = datetime(2016, 9, 30);
   expected = icemodel.createMetFileNames(opts);
   actual = icemodel.forcing.helpers.metfilename("kanm", "mar", ...
      opts.startdate, opts.enddate, 3600);
   testCase.verifyEqual(char(actual), expected{1});
end

%% windFromComponents

function test_windFromComponents_cardinal_directions(testCase)
   % Meteorological convention: direction the wind blows FROM, clockwise
   % from north, in (0, 360].

   u = [0; -1; 0; 1];   % toward: north, west, south, east
   v = [-1; 0; 1; 0];
   [wspd, wdir] = icemodel.forcing.helpers.windFromComponents(u, v);

   testCase.verifyEqual(wspd, ones(4, 1), 'AbsTol', 1e-12);
   testCase.verifyEqual(wdir, [360; 90; 180; 270], 'AbsTol', 1e-12);
end

%% metchecks

function test_metchecks_counts_and_fills_gaps(testCase)
   % NaN counts reflect the pre-fill state; gaps fill linearly with
   % nearest-value end fill.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 4, 0, 0))';
   tair = [270; 271; NaN; 273; NaN];
   met = timetable(Time, tair);

   [met, checks] = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(checks.numnan.tair, 2);
   testCase.verifyEqual(met.tair, [270; 271; 272; 273; 273], 'AbsTol', 1e-12);
end

function test_metchecks_clamps_ranges(testCase)
   % Recognized variables clamp to the documented physical ranges.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   albedo = [1.2; 0.5; 0.01];
   rh = [105; 50; 2];
   wspd = [0; 5; 5];
   tsfc = [280; 270; 271];   % kelvin (min > 100)
   met = timetable(Time, albedo, rh, wspd, tsfc);

   met = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(met.albedo, [0.98; 0.5; 0.05], 'AbsTol', 1e-12);
   testCase.verifyEqual(met.rh, [99.99; 50; 5], 'AbsTol', 1e-12);
   testCase.verifyEqual(met.wspd, [0.1; 5; 5], 'AbsTol', 1e-12);
   testCase.verifyEqual(met.tsfc, [273.16; 270; 271], 'AbsTol', 1e-12);
end

function test_metchecks_fills_wdir_circularly(testCase)
   % A gap between 350 and 10 degrees must fill near north (360/0), not
   % at the linear midpoint 180. This is the intentional fix relative to
   % the legacy runoff metchecks.

   Time = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 1, 2, 0, 0))';
   wdir = [350; NaN; 10];
   met = timetable(Time, wdir);

   met = icemodel.forcing.helpers.metchecks(met);

   testCase.verifyEqual(met.wdir(2), 360, 'AbsTol', 1e-9);
   testCase.verifyEqual(met.wdir([1 3]), [350; 10], 'AbsTol', 1e-9);
end

%% dailyToHourly

function test_dailyToHourly_linear_ramp(testCase)
   % A daily ramp interpolates to the exact hourly line, with linear
   % extrapolation past the final daily sample.

   t_daily = [datetime(2016, 1, 1); datetime(2016, 1, 2)];
   t_hourly = (datetime(2016, 1, 1):hours(1):datetime(2016, 1, 2, 23, 0, 0))';
   daily = [0; 24];

   hourly = icemodel.forcing.helpers.dailyToHourly(daily, t_daily, t_hourly);

   testCase.verifyEqual(hourly, (0:47)', 'AbsTol', 1e-9);
end

%% interpRcm

function test_interpRcm_recovers_planar_field(testCase)
   % Linear interpolation of a planar field is exact inside the hull,
   % and the timestep loop (interpolant-value reuse) must scale values.

   [X, Y] = meshgrid(0:4, 0:4);
   x = X(:);
   y = Y(:);
   v1 = 2*x + 3*y;
   v = [v1, 2*v1];   % cells x time
   xq = [1.5; 2.5];
   yq = [2.5; 0.5];

   V = icemodel.forcing.helpers.interpRcm(x, y, v, xq, yq, method="linear");

   expected = [2*xq + 3*yq, 2*(2*xq + 3*yq)]';   % time x cells
   testCase.verifyEqual(V, expected, 'AbsTol', 1e-9);
end

%% gridLocation point methods

function test_gridLocation_nearest_picks_nearest_cell(testCase)
   % "nearest" returns the single nearest cell; collapse yields its value.

   [X, Y] = ndgrid(0:6, 0:6);
   loc = [2.3, 3.7];   % nearest cell is (2, 4)
   [start, count, collapse, ~, loctype] = ...
      icemodel.forcing.helpers.gridLocation(X, Y, loc, "nearest");

   testCase.verifyEqual(loctype, "point");
   testCase.verifyEqual(count, [1 1]);
   testCase.verifyEqual([X(start(1), start(2)), Y(start(1), start(2))], [2 4]);

   block = 2 * X(start(1), start(2)) + 3 * Y(start(1), start(2));   % 1 cell
   testCase.verifyEqual(collapse(block), 16, 'AbsTol', 1e-12);
end

function test_gridLocation_natural_reproduces_planar_field(testCase)
   % "natural" collapses a neighbourhood by natural-neighbour interpolation,
   % which is exact for a linear field: f = 2x + 3y at (2.3, 3.7) = 15.7.

   [X, Y] = ndgrid(0:6, 0:6);
   loc = [2.3, 3.7];
   [start, count, collapse] = ...
      icemodel.forcing.helpers.gridLocation(X, Y, loc, "natural");

   rows = start(1):start(1) + count(1) - 1;
   cols = start(2):start(2) + count(2) - 1;
   xs = X(rows, cols);
   ys = Y(rows, cols);
   block = 2 * xs(:) + 3 * ys(:);   % planar field over the neighbourhood cells

   testCase.verifyEqual(collapse(block), 2 * 2.3 + 3 * 3.7, 'AbsTol', 1e-9);
end

%% conservative polygon remap (exactremap-guarded)

function test_remapPolygon_conservative_areaavg(testCase)
   % Conservative area-weighted remap of a linear field f = 10x + y over
   % the square [1,3]x[1,3]: the exact area-weighted mean of a linear
   % field equals its value at the centroid (2,2) = 22 (scaled x2 -> 44).
   % Inputs are ndgrid (the icemodel builder convention); remapPolygon
   % transposes to exactremap's meshgrid contract and reshapes the
   % cells-by-time block to the 3-D form exactremap weights correctly.
   testCase.assumeTrue(~isempty(which('exactremap')), ...
      'exactremap toolbox not on path');

   [X, Y] = ndgrid(0:5, 0:5);
   f = 10 * X + Y;
   block = [f(:), 2 * f(:)];   % two timesteps scaled [1, 2]
   P = polyshape([1 3 3 1], [1 1 3 3]);

   s = icemodel.forcing.helpers.remapPolygon(X, Y, block, P);
   testCase.verifyEqual(s, [22; 44], 'AbsTol', 1e-6);
end

function test_gridLocation_conservative_vs_equal(testCase)
   % The polygon collapse: "equal" weights only centre-in-polygon cells
   % (here just one -> 22), "conservative" area-weights the overlap (the
   % analytic linear-field mean over [1,3]^2, also 22 at the centroid).
   testCase.assumeTrue(~isempty(which('exactremap')), ...
      'exactremap toolbox not on path');

   [X, Y] = ndgrid(0:6, 0:6);
   P = polyshape([1 3 3 1], [1 1 3 3]);

   [s, c, collapse] = icemodel.forcing.helpers.gridLocation( ...
      X, Y, P, "nearest", remap="conservative");
   rows = s(1):s(1) + c(1) - 1;
   cols = s(2):s(2) + c(2) - 1;
   xs = X(rows, cols);
   ys = Y(rows, cols);
   block = 10 * xs(:) + ys(:);
   v = collapse(block);
   testCase.verifyEqual(v, 22, 'AbsTol', 1e-6);
end

%% validatemet

function test_validatemet_rejects_missing_variable(testCase)
   % Required met variables must be present.

   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met = removevars(met, 'psfc');

   testCase.verifyError(@() icemodel.forcing.helpers.validatemet(met), ...
      'icemodel:forcing:validatemet:missingVariables');
end

function test_validatemet_rejects_irregular_time(testCase)
   % The time axis must have a uniform timestep.

   met = makeSyntheticMet(datetime(2016, 1, 1), 24);
   met.Time(3) = met.Time(3) + minutes(1);

   testCase.verifyError(@() icemodel.forcing.helpers.validatemet(met), ...
      'icemodel:forcing:validatemet:irregularTimeAxis');
end

%% writemet

function test_writemet_window_writes_and_loads(testCase)
   % Window naming produces one file holding the full timetable as met.

   met = makeSyntheticMet(datetime(2016, 1, 1), 48);

   filenames = icemodel.forcing.helpers.writemet(met, "tst", "src", ...
      outdir=testCase.TestData.outdir);

   testCase.verifyEqual(numel(filenames), 1);
   [~, name, ext] = fileparts(filenames);
   testCase.verifyEqual(name + ext, "met_tst_src_20160101_20160102_1hr.mat");

   loaded = load(filenames, 'met');
   testCase.verifyEqual(loaded.met, met);
end

function test_writemet_yearly_splits_years(testCase)
   % Yearly naming produces one per-year file with that year's rows.

   met = makeSyntheticMet(datetime(2015, 12, 31), 48);   % spans 2015/2016

   filenames = icemodel.forcing.helpers.writemet(met, "tst", "src", ...
      outdir=testCase.TestData.outdir, naming="yearly");

   testCase.verifyEqual(numel(filenames), 2);
   testCase.verifyTrue(endsWith(filenames(1), "met_tst_src_2015_1hr.mat"));
   testCase.verifyTrue(endsWith(filenames(2), "met_tst_src_2016_1hr.mat"));

   loaded = load(filenames(2), 'met');
   testCase.verifyEqual(unique(year(loaded.met.Time)), 2016);
   testCase.verifyEqual(height(loaded.met), 24);
end

%% writeuserdata

function test_writeuserdata_per_year_with_metadata(testCase)
   % Data files split per year, keep the Data variable name, and carry
   % the location CustomProperties.

   Data = makeSyntheticData(datetime(2015, 12, 31), 48);

   filenames = icemodel.forcing.helpers.writeuserdata(Data, "tst", "src", ...
      outdir=testCase.TestData.outdir);

   testCase.verifyEqual(numel(filenames), 2);
   testCase.verifyTrue(endsWith(filenames(1), "tst_src_2015.mat"));
   testCase.verifyTrue(endsWith(filenames(2), "tst_src_2016.mat"));

   loaded = load(filenames(1), 'Data');
   testCase.verifyEqual(unique(year(loaded.Data.Time)), 2015);
   testCase.verifyEqual(loaded.Data.Properties.CustomProperties.Lat, 67.0);
end

function test_writeuserdata_rejects_missing_metadata(testCase)
   % The CustomProperties location-metadata contract is enforced.

   Data = makeSyntheticData(datetime(2016, 1, 1), 24);
   Data = removeCustomProperties(Data);

   testCase.verifyError(@() icemodel.forcing.helpers.writeuserdata( ...
      Data, "tst", "src", outdir=testCase.TestData.outdir), ...
      'icemodel:forcing:writeuserdata:missingMetadata');
end

%% psnProjection

function test_psnProjection_matches_legacy_projsipsn(testCase)
   % psnProjection (EPSG:3413) must reproduce the legacy projsipsn.mat
   % forward projection. Skips when matfunclib (which carries
   % projsipsn.mat) is not on the path.

   legacy_file = which('projsipsn.mat');
   testCase.assumeTrue(~isempty(legacy_file), ...
      'projsipsn.mat not on path (matfunclib not loaded)');

   legacy = load(legacy_file).('projsipsn');
   lat = [61.1; 67.0; 72.5];
   lon = [-48.3; -49.5; -40.0];

   [x_old, y_old] = projfwd(legacy, lat, lon);
   [x_new, y_new] = projfwd(icemodel.forcing.helpers.psnProjection(), lat, lon);

   testCase.verifyEqual(x_new, x_old, 'AbsTol', 0.01);
   testCase.verifyEqual(y_new, y_old, 'AbsTol', 0.01);
end

%% Local fixture helpers

function met = makeSyntheticMet(t0, nsteps)
   %MAKESYNTHETICMET Hourly timetable satisfying the met contract.
   Time = (t0:hours(1):(t0 + hours(nsteps - 1)))';
   n = numel(Time);
   tair = 265 + 5*sin(2*pi*(1:n)'/24);
   swd = max(0, 400*sin(2*pi*(1:n)'/24));
   lwd = 250*ones(n, 1);
   albedo = 0.8*ones(n, 1);
   wspd = 4 + sin(2*pi*(1:n)'/24);
   rh = 70*ones(n, 1);
   psfc = 1e5*ones(n, 1);
   ppt = zeros(n, 1);
   met = timetable(Time, tair, swd, lwd, albedo, wspd, rh, psfc, ppt);
end

function Data = makeSyntheticData(t0, nsteps)
   %MAKESYNTHETICDATA Data timetable with location CustomProperties.
   Data = makeSyntheticMet(t0, nsteps);
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = -2.1e5;
   Data.Properties.CustomProperties.Y = -2.5e6;
   Data.Properties.CustomProperties.Lat = 67.0;
   Data.Properties.CustomProperties.Lon = -49.5;
   Data.Properties.CustomProperties.Elev = 1270;
   Data.Properties.CustomProperties.Slope = 0.01;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];
end

function Data = removeCustomProperties(Data)
   %REMOVECUSTOMPROPERTIES Rebuild DATA without CustomProperties.
   Data = timetable(Data.Time, Data.tair, 'VariableNames', {'tair'});
end
