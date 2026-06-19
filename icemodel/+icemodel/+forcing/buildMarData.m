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
   %  - location = polyshape (vertices in EPSG:3413 meters): the MAR cells
   %    are averaged over the polygon. remap="conservative" (default) uses
   %    exact overlap-area weighting via the exactremap toolbox; remap="equal"
   %    is a plain mean of in-polygon cell centers. The remap runs in the
   %    NATIVE MAR projection (a regular 15 km grid); the polygon is mapped
   %    from EPSG:3413 into native coordinates via the shipped LON/LAT.
   %
   % Channels (standard units; daily MAR channels interpolated hourly):
   %  hourly: tair [K], shum [kg/kg] (dropped after rh derivation), swd,
   %  lwd, shf, lhf [W m-2], albedo [-], snow, rain, melt, runoff, smb
   %  [mWE/h]; daily: snowd [m], cfrac [-], tsfc [K], psfc [Pa]; derived:
   %  wspd [m s-1], wdir [deg] (from UUH/VVH), rh [%] (icemodel.vapor
   %  kernel); optional: modis [-] (GEUS MODIS daily albedo).
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
   %      added (point location only). Reference layout:
   %      /Volumes/S03/DATA/greenland/geus/albedo/gris.
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
      source_dir = string(fullfile(icemodel.getpath('data'), ...
         'forcing', 'mar'));
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
      'snow', 'SFH'  ; 'rain', 'RFH'  ; 'melt', 'MEH'  ; 'runoff', 'RUH'
      'shf',  'SHFH' ; 'lhf',  'LHFH' ; 'smb',  'SMBH'
      };
   daily_vars = {
      'snowd', 'SHSN2'; 'cfrac', 'CC'; 'tsfc', 'ST'; 'psfc', 'SP'
      };

   parts = cell(numel(years), 1);
   files = strings(numel(years), 1);
   for n = 1:numel(years)
      files(n) = locateMarFile(source_dir, years(n));
      if n == 1
         [start, count, collapse, site] = ...
            resolveLocation(files(1), location, kwargs.method, kwargs.remap);
      end
      parts{n} = extractOneYear(files(n), hourly_vars, daily_vars, ...
         start, count, collapse);
      if kwargs.modis_dir ~= ""
         parts{n}.modis = modisChannel(kwargs.modis_dir, years(n), ...
            site.lat, site.lon, parts{n}.Time);
      end
   end
   Data = vertcat(parts{:});

   % Derived channels: wind from components, RH from the canonical vapor
   % kernel; the component and humidity inputs then drop out.
   [Data.wspd, Data.wdir] = icemodel.forcing.helpers.windFromComponents( ...
      Data.uwind, Data.vwind);
   Data.rh = icemodel.vapor.relative_humidity_from_specific_humidity( ...
      Data.shum, Data.psfc, Data.tair);
   Data = removevars(Data, {'shum', 'uwind', 'vwind'});

   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);

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

function [start, count, collapse, site] = resolveLocation( ...
      filename, location, method, remap)
   %RESOLVELOCATION Map a point or polygon onto a MAR grid hyperslab.
   %
   % Returns the bounding hyperslab (start/count over the grid dims), the
   % collapse function handle that reduces a hyperslab block to the target
   % series (nearest cell, natural-neighbour point, equal-weight polygon
   % mean, or conservative area-weighted polygon remap), and the site
   % summary (nearest cell / in-polygon mean metadata). For the
   % conservative polygon remap the MAR ice mask (SRF == 4) is passed as
   % the valid-cells mask so off-ice cells are inpainted from on-ice
   % neighbours.
   %
   % Spatial selection/remap is done in the NATIVE MAR projection, where the
   % grid is exactly regular (the EPSG:3413 reprojection is curvilinear and
   % would be rejected as irregular by the conservative remap). The query
   % (point or polygon, given as [lat lon] / EPSG:3413) is mapped into native
   % coordinates with the shipped LON/LAT <-> native correspondence.
   grid = icemodel.forcing.marGridInfo(filename);
   proj = icemodel.forcing.helpers.psnProjection();
   toNativeX = scatteredInterpolant(grid.LON(:), grid.LAT(:), ...
      grid.Xnat(:), 'natural', 'nearest');
   toNativeY = scatteredInterpolant(grid.LON(:), grid.LAT(:), ...
      grid.Ynat(:), 'natural', 'nearest');

   if isnumeric(location)
      assert(isequal(size(location), [1 2]), ...
         'point location must be [lat lon]')
      lat = location(1);
      lon = location(2);
      location = [toNativeX(lon, lat), toNativeY(lon, lat)];
   elseif isa(location, 'polyshape')
      [vlat, vlon] = projinv(proj, location.Vertices(:, 1), ...
         location.Vertices(:, 2));
      location = polyshape(toNativeX(vlon, vlat), toNativeY(vlon, vlat));
   end
   [start, count, collapse, inslab, site.type] = ...
      icemodel.forcing.helpers.gridLocation(grid.Xnat, grid.Ynat, ...
      location, method, remap=remap, validmask=(grid.srf == 4));

   slabmean = @(field) icemodel.forcing.helpers.slabMean( ...
      field, start, count, inslab);
   site.lat = slabmean(grid.LAT);
   site.lon = slabmean(grid.LON);
   site.elev = slabmean(grid.elev);
   site.slope = slabmean(grid.slope);
   [site.x, site.y] = projfwd(proj, site.lat, site.lon);
   site.srf_warning = any(slabmean(grid.srf) ~= 4);
end

function part = extractOneYear(filename, hourly_vars, daily_vars, ...
      start, count, collapse)
   %EXTRACTONEYEAR Read and assemble one MAR year at the target location.
   % COLLAPSE (from gridLocation) reduces each variable's hyperslab block
   % (cells x time) to the target series (nearest cell, natural-neighbour
   % point, or polygon mean).

   for n = 1:size(hourly_vars, 1)
      [data, ~, Time] = icemodel.forcing.readMar3p11(filename, ...
         hourly_vars{n, 2}, start=start, count=count);
      if n == 1
         part = timetable(Time);
      end
      part.(hourly_vars{n, 1}) = collapse(data);
   end

   for n = 1:size(daily_vars, 1)
      [data, ~, Tdaily] = icemodel.forcing.readMar3p11(filename, ...
         daily_vars{n, 2}, start=start, count=count);
      part.(daily_vars{n, 1}) = icemodel.forcing.helpers.dailyToHourly( ...
         collapse(data), Tdaily, part.Time);
   end
end

function modis = modisChannel(modis_dir, yyyy, lat, lon, Time)
   %MODISCHANNEL Daily GEUS MODIS albedo interpolated to the time axis.
   match = dir(fullfile(modis_dir, sprintf('*_%d_*.nc', yyyy)));
   if numel(match) ~= 1
      error('icemodel:forcing:buildMarData:modisNotFound', ...
         'expected one MODIS file for %d in %s, found %d', ...
         yyyy, modis_dir, numel(match))
   end
   [albedo, Tdaily] = icemodel.forcing.readGeusModis( ...
      string(fullfile(match.folder, match.name)), lat, lon);
   modis = icemodel.forcing.helpers.dailyToHourly(albedo, Tdaily, Time);
end
