function [Data, metadata] = buildMarData(location, years, kwargs)
   %BUILDMARDATA Build a Data timetable from MAR v3.11 yearly NetCDF files.
   %
   %  [Data, metadata] = icemodel.forcing.buildMarData(location, years)
   %  [Data, metadata] = ... buildMarData(_, source_dir=..., ...
   %     modis_dir=..., fillgaps=true)
   %
   % Extracts the icemodel forcing/evaluation channels from MAR v3.11
   % for any years available in the source directory, at a point or
   % averaged over a polygon:
   %
   %  - location = [lat lon] (1x2, degrees): the nearest MAR cell is
   %    extracted directly (the legacy point workflow that produced the
   %    ak4/behar met and userdata artifacts).
   %  - location = Nx2 [lat lon] (N>1): a LIST of points. Returns a 1xN cell
   %    of Data timetables (metadata is a 1xN struct array). The grid is read
   %    ONCE and every yearly file is opened ONCE, slicing each point's
   %    hyperslab from that single open. N=1 is the single-point path.
   %  - location = polyshape (vertices in EPSG:3413 meters): the MAR cells
   %    are averaged over the polygon. remap="conservative" (default) uses
   %    exact overlap-area weighting via the exactremap toolbox; remap="equal"
   %    is a plain mean of in-polygon cell centers. The remap runs in the
   %    NATIVE MAR projection (a regular 15 km grid); the polygon is mapped
   %    from EPSG:3413 into native coordinates via the shipped LON/LAT.
   %
   % Channels (canonical units; daily MAR channels interpolated hourly):
   %  hourly: tair [K], shum [kg/kg] (dropped after rh derivation), swd,
   %  lwd, shf, lhf [W m-2], albedo [-], snowf, rainf [m s-1, the canonical
   %  water-equivalent precipitation rates], melt, runoff, smb [mWE/h];
   %  daily: snowd [m], cfrac [-], tsfc [K], psfc [Pa]; derived: wspd
   %  [m s-1], wdir [deg] (from UUH/VVH), rh [%] (icemodel.vapor kernel);
   %  optional: modis [-] (GEUS MODIS daily albedo).
   %
   % Inputs
   %  location - [lat lon] point or polyshape (see above)
   %  years    - calendar years to extract (one MAR file per year)
   %
   % Name-value
   %  source_dir : directory holding MAR yearly files (*-<YYYY>.nc).
   %      Defaults to the gitignored cache data/forcing/mar. Reference
   %      layout: /Volumes/S03/DATA/greenland/mar3p11/RUH2.
   %  modis_dir  : directory with GEUS Greenland_Reflectivity_<YYYY>_
   %      5km_C6.nc files; when given, a daily MODIS albedo channel is
   %      added at the requested location (point: nearest/natural; polygon:
   %      conservative/equal catchment mean, same as the gridded channels).
   %      Reference layout: /Volumes/S03/DATA/greenland/geus/albedo/gris.
   %  fillgaps   : gap-fill through metchecks (default true, the legacy
   %      RCM-Data behavior; MAR output is gap-free in practice)
   %
   % Outputs
   %  Data     - hourly timetable with userdata CustomProperties (X, Y,
   %             Lat, Lon, Elev, Slope, ScalarUnits)
   %  metadata - provenance: files read, cell/polygon info, policies
   %
   % Observation heights (important for the turbulent-flux scheme): MAR's
   % hourly diagnostics are at the standard meteorological heights -
   % temperature/humidity (TTH/QQH -> tair/rh) at 2 m, wind (UUH/VVH ->
   % wspd) at 10 m. They are intentionally on different levels. The model
   % must be told these heights: icemodel.setopts sets opts.z_tair = 2 and
   % opts.z_wind = 10 for forcings = "mar" (z_relh = z_tair). If you build
   % MAR forcing for a custom run, keep opts.z_wind = 10 / opts.z_tair = 2.
   %
   % Legacy: reimplements runoff/functions/saveMarData.m (the original
   % retained, unchanged, as the legacy reference workflow). Derivable
   % radiation terms (swu, lwu, swn, lwn, netr) that the legacy
   % computeDerivedValues stored are NOT carried here - icemodel.processmet
   % recomputes them on load from swd/albedo/tsfc/lwd.
   %
   % See also: icemodel.forcing.readMar3p11, icemodel.forcing.data2met,
   %  icemodel.forcing.buildMarMet, icemodel.forcing.helpers.writeuserdata

   arguments
      location
      years (1, :) double {mustBeInteger}
      kwargs.source_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.fillgaps (1, 1) logical = true
      kwargs.method (1, 1) string {mustBeMember(kwargs.method, ...
         ["nearest", "natural"])} = "nearest"
      kwargs.remap (1, 1) string {mustBeMember(kwargs.remap, ...
         ["equal", "conservative"])} = "conservative"
   end

   source_dir = kwargs.source_dir;
   if source_dir == ""
      source_dir = string(getenv("ICEMODEL_MAR_DIR"));
      if source_dir == ""
         source_dir = "/Volumes/S03/DATA/greenland/mar3p11/RUH2";
      end
   end
   if ~isfolder(source_dir)
      error('icemodel:forcing:buildMarData:sourceNotFound', ...
         ['MAR source directory not found: %s. Pass source_dir or stage ' ...
         'the MAR yearly files (reference layout: ' ...
         '/Volumes/S03/DATA/greenland/mar3p11/RUH2).'], source_dir)
   end

   % MAR channel table: output name <- MAR variable.
   hourly_vars = {
      'tair', 'TTH'  ; 'shum', 'QQH'  ; 'uwind', 'UUH' ; 'vwind', 'VVH'
      'swd',  'SWDH' ; 'lwd',  'LWDH' ; 'albedo', 'ALH'
      'snowf', 'SFH' ; 'rainf', 'RFH' ; 'melt', 'MEH'  ; 'runoff', 'RUH'
      'shf',  'SHFH' ; 'lhf',  'LHFH' ; 'smb',  'SMBH'
      };
   daily_vars = {
      'snowd', 'SHSN2'; 'cfrac', 'CC'; 'tsfc', 'ST'; 'psfc', 'SP'
      };

   % Accept one location (1x2 point or polyshape, returns a single Data
   % timetable) OR a list of N points (Nx2 [lat lon], returns a 1xN cell
   % array of Data timetables). N=1 is just the single-point path, so there
   % is ONE code path: resolve every point's grid hyperslab up front (grid
   % metadata + interpolants read ONCE), then loop years opening each yearly
   % file ONCE and reading every point's hyperslab from that single open.
   [locations, batch] = locationList(location);
   npts = numel(locations);

   files = strings(numel(years), 1);
   for n = 1:numel(years)
      files(n) = locateMarFile(source_dir, years(n));
   end

   % Per-point grid hyperslab + collapse rule + site metadata, resolved once
   % against a single grid read (marGridInfo / scatteredInterpolant built
   % ONCE for the whole point list rather than per point per year).
   grid = resolveGrid(files(1), kwargs.method, kwargs.remap);
   [slabs, collapses, sites] = deal(cell(1, npts), cell(1, npts), cell(1, npts));
   for p = 1:npts
      [slabs{p}, collapses{p}, sites{p}] = resolvePoint(grid, locations{p});
   end

   % Read each yearly file ONCE per variable, slicing every point's hyperslab
   % from that single open, then assemble per point.
   parts = cell(numel(years), npts);
   for n = 1:numel(years)
      blocks = extractOneYear(files(n), hourly_vars, daily_vars, slabs);
      for p = 1:npts
         parts{n, p} = assemblePart(blocks, hourly_vars, daily_vars, ...
            collapses{p}, p);
      end
   end

   data_out = cell(1, npts);
   meta_out = cell(1, npts);
   for p = 1:npts
      [data_out{p}, meta_out{p}] = finalizeMarData(vertcat(parts{:, p}), ...
         files, slabs{p}, sites{p}, years, locations{p}, kwargs);
   end

   % A single location returns a single Data timetable + metadata struct; a
   % point list returns a 1xN cell of timetables + a 1xN metadata struct array.
   metadata = [meta_out{:}];
   if batch
      Data = data_out;
   else
      Data = data_out{1};
   end
end

%% Local functions
function filename = locateMarFile(source_dir, yyyy)
   %LOCATEMARFILE Resolve the MAR yearly file for one calendar year.
   match = dir(fullfile(source_dir, sprintf('*-%d.nc', yyyy)));
   if numel(match) ~= 1
      error('icemodel:forcing:buildMarData:fileNotFound', ...
         'expected one MAR file matching *-%d.nc in %s, found %d', ...
         yyyy, source_dir, numel(match))
   end
   filename = string(fullfile(match.folder, match.name));
end

function [locations, batch] = locationList(location)
   %LOCATIONLIST Normalize the location input to a 1xN cell of locations.
   % A polyshape or a single [lat lon] row is one location (batch=false,
   % single Data timetable). An Nx2 [lat lon] (N>1) is a point list
   % (batch=true, 1xN cell of Data timetables); a single row stays scalar.
   if isnumeric(location) && size(location, 2) == 2 && size(location, 1) > 1
      locations = num2cell(location, 2)';
      batch = true;
   else
      locations = {location};
      batch = false;
   end
end

function grid = resolveGrid(filename, method, remap)
   %RESOLVEGRID Read the MAR grid + native-coordinate interpolants ONCE.
   %
   % Spatial selection/remap is done in the NATIVE MAR projection, where the
   % grid is exactly regular (the EPSG:3413 reprojection is curvilinear and
   % would be rejected as irregular by the conservative remap). The query
   % (point or polygon, given as [lat lon] / EPSG:3413) is mapped into native
   % coordinates with the paired LON/LAT and Xnat/Ynat values carried by each
   % MAR grid cell. The grid and native-coordinate maps are built once and
   % reused for every point in a batch.
   grid = icemodel.forcing.marGridInfo(filename);
   grid.method = method;
   grid.remap = remap;
   grid.proj = icemodel.forcing.helpers.psnProjection();
   grid.toNativeX = scatteredInterpolant(grid.LON(:), grid.LAT(:), ...
      grid.Xnat(:), 'natural', 'nearest');
   grid.toNativeY = scatteredInterpolant(grid.LON(:), grid.LAT(:), ...
      grid.Ynat(:), 'natural', 'nearest');
end

function [slab, collapse, site] = resolvePoint(grid, location)
   %RESOLVEPOINT Map one point or polygon onto a MAR grid hyperslab.
   %
   % Returns the bounding hyperslab as a [start; count] 2x2 (for the batch
   % reader), the collapse function handle that reduces a hyperslab block to
   % the target series (nearest cell, natural-neighbour point, equal-weight
   % polygon mean, or conservative area-weighted polygon remap), and the site
   % summary (nearest cell / in-polygon mean metadata). For the conservative
   % polygon remap the MAR ice mask (SRF == 4) is passed as the valid-cells
   % mask so off-ice cells are inpainted from on-ice neighbours.
   if isnumeric(location)
      assert(isequal(size(location), [1 2]), ...
         'point location must be [lat lon]')
      lat = location(1);
      lon = location(2);
      location = [grid.toNativeX(lon, lat), grid.toNativeY(lon, lat)];
   elseif isa(location, 'polyshape')
      [vlat, vlon] = projinv(grid.proj, location.Vertices(:, 1), ...
         location.Vertices(:, 2));
      location = polyshape(grid.toNativeX(vlon, vlat), ...
         grid.toNativeY(vlon, vlat));
   end
   [start, count, collapse, inslab, site.type] = ...
      icemodel.forcing.helpers.gridLocation(grid.Xnat, grid.Ynat, ...
      location, grid.method, remap=grid.remap, validmask=(grid.srf == 4));
   slab = [start; count];

   slabmean = @(field) icemodel.forcing.helpers.slabMean( ...
      field, start, count, inslab);
   site.lat = slabmean(grid.LAT);
   site.lon = slabmean(grid.LON);
   site.elev = slabmean(grid.elev);
   site.slope = slabmean(grid.slope);
   [site.x, site.y] = projfwd(grid.proj, site.lat, site.lon);
   site.srf_warning = any(slabmean(grid.srf) ~= 4);
end

function blocks = extractOneYear(filename, hourly_vars, daily_vars, slabs)
   %EXTRACTONEYEAR Read one MAR year, every point's hyperslab per variable.
   % Opens each yearly file ONCE per variable (the batch reader slices every
   % point's hyperslab from a single open) and returns the raw cells-by-time
   % blocks per variable per point plus the time axes, deferring the per-point
   % collapse to assemblePart. A single point is just a one-element slab list.
   blocks = struct('hourly', {cell(size(hourly_vars, 1), 1)}, ...
      'daily', {cell(size(daily_vars, 1), 1)}, 'Time', [], 'Tdaily', []);

   for n = 1:size(hourly_vars, 1)
      [data, ~, Time] = icemodel.forcing.readMar3p11(filename, ...
         hourly_vars{n, 2}, slabs=slabs);
      blocks.hourly{n} = data;
      if n == 1
         blocks.Time = Time;
      end
   end

   for n = 1:size(daily_vars, 1)
      [data, ~, Tdaily] = icemodel.forcing.readMar3p11(filename, ...
         daily_vars{n, 2}, slabs=slabs);
      blocks.daily{n} = data;
      if n == 1
         blocks.Tdaily = Tdaily;
      end
   end
end

function part = assemblePart(blocks, hourly_vars, daily_vars, collapse, p)
   %ASSEMBLEPART Collapse one point's blocks into one MAR year timetable.
   % COLLAPSE (from gridLocation) reduces each variable's hyperslab block
   % (cells x time) to the target series (nearest cell, natural-neighbour
   % point, or polygon mean). p indexes the point in each block's cell list.
   part = timetable(blocks.Time);
   for n = 1:size(hourly_vars, 1)
      part.(hourly_vars{n, 1}) = collapse(blocks.hourly{n}{p});
   end
   for n = 1:size(daily_vars, 1)
      part.(daily_vars{n, 1}) = icemodel.forcing.helpers.dailyToHourly( ...
         collapse(blocks.daily{n}{p}), blocks.Tdaily, part.Time);
   end
end

function [Data, metadata] = finalizeMarData(Data, files, slab, site, ...
      years, location, kwargs)
   %FINALIZEMARDATA Post-process one point's assembled MAR Data + metadata.
   % Identical to the legacy single-point tail: optional MODIS channel,
   % derived wind/RH, precip rate, metchecks, units, userdata CustomProperties,
   % and the provenance struct. start/count come from the point's slab.
   start = slab(1, :);
   count = slab(2, :);

   % Optional GEUS MODIS albedo at the requested location: nearest/natural for
   % a point, conservative (or equal) catchment mean for a polygon (the same
   % selection as the gridded channels). Added on the assembled axis via the
   % shared helper so MAR/MERRA/RACMO resolve MODIS identically.
   if kwargs.modis_dir ~= ""
      Data.modis = icemodel.forcing.helpers.modisAlbedoChannel( ...
         kwargs.modis_dir, years, location, kwargs.method, kwargs.remap, ...
         Data.Time);
   end

   % Derived channels: wind from components, RH from the canonical vapor
   % kernel; the component and humidity inputs then drop out.
   [Data.wspd, Data.wdir] = icemodel.forcing.helpers.windFromComponents( ...
      Data.uwind, Data.vwind);
   Data.rh = icemodel.vapor.relative_humidity_from_specific_humidity( ...
      Data.shum, Data.psfc, Data.tair);
   Data = removevars(Data, {'shum', 'uwind', 'vwind'});

   % Precipitation to the canonical water-equivalent rate m s-1. MAR posts
   % snowfall/rainfall as mWE/h (an hourly water-equivalent depth rate), so
   % dividing by 3600 s/h yields m s-1. Channels use the canonical snowf/rainf
   % names (icemodel.forcing.helpers.metvariables optional split; data2met
   % derives ppt = rainf + snowf). The diagnostic mass fluxes (melt/runoff/smb)
   % keep their natural mWE/h rate.
   for ch = ["snowf", "rainf"]
      Data.(ch) = Data.(ch) / 3600;
   end

   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);

   % Per-variable units from the shared canonical map. Precipitation is m s-1
   % (converted above); the diagnostic mass fluxes are mWE/h rates (cumulative
   % sums need the timestep in hours).
   Data.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(Data.Properties.VariableNames));

   % Userdata CustomProperties.
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = site.x;
   Data.Properties.CustomProperties.Y = site.y;
   Data.Properties.CustomProperties.Lat = site.lat;
   Data.Properties.CustomProperties.Lon = site.lon;
   Data.Properties.CustomProperties.Elev = site.elev;
   Data.Properties.CustomProperties.Slope = site.slope;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];

   metadata = struct( ...
      'source_files', files, ...
      'location_type', site.type, ...
      'method', kwargs.method, ...
      'remap', kwargs.remap, ...
      'grid_start', start, ...
      'grid_count', count, ...
      'n_cells', prod(count), ...
      'lat', site.lat, 'lon', site.lon, ...
      'elev', site.elev, ...
      'humidity_kernel', ...
      "icemodel.vapor.relative_humidity_from_specific_humidity", ...
      'checks', checks);
   if site.srf_warning
      warning('icemodel:forcing:buildMarData:surfaceNotIce', ...
         'MAR surface type at the requested location is not ice sheet')
   end
end
