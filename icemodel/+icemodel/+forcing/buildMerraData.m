function [Data, metadata] = buildMerraData(location, years, kwargs)
   %BUILDMERRADATA Build a Data timetable from MERRA-2 daily NetCDF files.
   %
   %  [Data, metadata] = icemodel.forcing.buildMerraData(location, years)
   %  [Data, metadata] = ... buildMerraData(_, source_dir=..., ...
   %     modis_dir=..., fillgaps=true)
   %
   % Extracts the icemodel forcing/evaluation channels from the MERRA-2
   % single-level collections at a point or averaged over a polygon. The
   % archive is organized as one NetCDF per day per collection
   % (MERRA2_*.tavg1_2d_<col>_Nx.YYYYMMDD.nc4*), with collections:
   %
   %    slv (hourly): T2M -> tair [K], QV2M -> shum [kg/kg] (dropped
   %        after rh derivation), U2M/V2M -> wspd/wdir, PS -> psfc [Pa]
   %    rad (hourly): SWGDN -> swd, SWGNT -> swn, LWGAB -> lwd,
   %        LWGNT -> lwn [W m-2]
   %    flx (hourly): HFLUX -> shf, EFLUX -> lhf [W m-2],
   %        PRECTOTCORR -> ppt, PRECSNO -> snowf, EVAP -> evap [mWE/h]
   %    glc (3-hourly): RUNOFF -> runoff [mWE/h], SNICEALB -> albedo,
   %        SNOWDP_GL -> snowd [m], SNOMAS_GL -> swe [kg m-2]
   %
   % Derived: wspd/wdir from U2M/V2M; rh [%] via the canonical icemodel.vapor
   % kernel. The derivable radiation terms (swu, lwu, netr) are NOT stored -
   % icemodel.processmet recomputes them on load from swd/albedo/tsfc/lwd.
   % Mass fluxes convert kg m-2 s-1 -> meters water equivalent per hour.
   % MERRA-2 time-average stamps sit at the bin center (00:30 hourly, 01:30
   % 3-hourly); everything is interpolated onto the on-the-hour axis, the
   % legacy convention.
   %
   % Legacy: reimplements runoff/functions/saveMerraData.m (the original
   % retained, unchanged, as the legacy reference workflow). Note the legacy
   % code reconstructed swd as SWGNT/(1-SNICEALB); this builder reads the
   % native SWGDN downwelling flux directly.
   %
   % The available period is whatever days exist in the source
   % directory - the calendar derives from the files themselves (this
   % replaces the legacy hardcoded 2008-2020 calendar).
   %
   % Inputs
   %  location - [lat lon] point or polyshape in EPSG:3413 meters
   %  years    - calendar years to extract
   %
   % Name-value
   %  source_dir : directory holding the collection subdirectories
   %      (flx/, glc/, rad/, slv/). Defaults to the gitignored cache
   %      data/forcing/merra2. Reference layout:
   %      /Volumes/S03/DATA/merra2/1hrly/ncfiles.
   %  modis_dir : optional GEUS MODIS albedo directory (see buildMarData)
   %  fillgaps  : gap-fill through metchecks (default true)
   %
   % Outputs
   %  Data     - hourly timetable with userdata CustomProperties
   %  metadata - provenance: collections, file counts, cell info
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.buildRacmoData, icemodel.forcing.data2met

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
         'forcing', 'merra2'));
   end
   if ~isfolder(source_dir)
      error('icemodel:forcing:buildMerraData:sourceNotFound', ...
         ['MERRA-2 source directory not found: %s. Pass source_dir or ' ...
         'stage the collection subdirectories (reference layout: ' ...
         '/Volumes/S03/DATA/merra2/1hrly/ncfiles).'], source_dir)
   end

   % Channel map: {collection, MERRA variable, output name}.
   channels = {
      'slv', 'T2M',         'tair'
      'slv', 'QV2M',        'shum'
      'slv', 'U2M',         'uwind'
      'slv', 'V2M',         'vwind'
      'slv', 'PS',          'psfc'
      'rad', 'SWGDN',       'swd'
      'rad', 'SWGNT',       'swn'
      'rad', 'LWGAB',       'lwd'
      'rad', 'LWGNT',       'lwn'
      'flx', 'HFLUX',       'shf'
      'flx', 'EFLUX',       'lhf'
      'flx', 'PRECTOTCORR', 'ppt'
      'flx', 'PRECSNO',     'snowf'
      'flx', 'EVAP',        'evap'
      'glc', 'RUNOFF',      'runoff'
      'glc', 'SNICEALB',    'albedo'
      'glc', 'SNOWDP_GL',   'snowd'
      'glc', 'SNOMAS_GL',   'swe'
      };
   collections = unique(channels(:, 1), 'stable');

   % File inventory per collection, keyed by the YYYYMMDD name token.
   % The calendar derives from the files present (no hardcoded period).
   inventory = struct();
   for c = collections'
      files = dir(fullfile(source_dir, c{1}, '*_Nx.*.nc4*'));
      assert(~isempty(files), ...
         'no MERRA-2 %s files found under %s', c{1}, source_dir)
      tokens = regexp({files.name}, '_Nx\.(\d{8})\.', 'tokens', 'once');
      dates = datetime(string(cellfun(@(t) t{1}, tokens, ...
         'UniformOutput', false)), 'InputFormat', 'yyyyMMdd', ...
         'TimeZone', 'UTC');
      keep = ismember(year(dates), years);
      if ~any(keep)
         error('icemodel:forcing:buildMerraData:yearNotInArchive', ...
            'no MERRA-2 %s files for years %s (archive: %d-%d)', ...
            c{1}, mat2str(years), year(min(dates)), year(max(dates)))
      end
      dates = dates(:);
      [inventory.(c{1}).dates, order] = sort(dates(keep));
      kept = files(keep);
      inventory.(c{1}).files = string(fullfile( ...
         {kept(order).folder}, {kept(order).name}))';
   end

   % Grid and target cells from the first slv file (all collections
   % share the clipped regular lat/lon grid). MERRA-2 variables read as
   % [nlon nlat ntime], so the coordinate grids use the same [lon lat]
   % orientation to keep gridLocation indices aligned with ncread.
   first = inventory.slv.files(1);
   lat = double(ncread(first, 'lat'));
   lon = double(ncread(first, 'lon'));
   [LON, LAT] = ndgrid(lon, lat);
   proj = icemodel.forcing.helpers.psnProjection();
   [X, Y] = projfwd(proj, LAT, LON);

   % Ice mask for conservative remap. MERRA-2 glacier-tile (glc) variables
   % are valid only over land-ice cells (non-glacier cells carry the
   % _FillValue), so SNOMAS_GL validity gives a static land-ice mask; off-
   % ice cells are inpainted from on-ice neighbours. (Degrades safely to
   % all-valid if a build of MERRA ever stores 0 rather than fill there.)
   validmask = [];
   if kwargs.remap == "conservative"
      validmask = merraIceMask(inventory.glc.files(1), size(X));
   end

   % Conservative polygon remap runs in MERRA's NATIVE geographic grid
   % (regular lon/lat) with exactremap UseGeoCoords=true, which computes
   % true ellipsoidal overlap areas - the correct conservative weighting for
   % a lat/lon grid (reprojecting to EPSG:3413 first would make the grid
   % irregular). Point/nearest and equal-weight stay in the projected grid.
   if isa(location, 'polyshape') && kwargs.remap == "conservative"
      [vlat, vlon] = projinv(proj, location.Vertices(:, 1), ...
         location.Vertices(:, 2));
      [start, count, collapse, inslab, loctype] = ...
         icemodel.forcing.helpers.gridLocation(LON, LAT, ...
         polyshape(vlon, vlat), kwargs.method, remap="conservative", ...
         validmask=validmask, usegeocoords=true);
   else
      if isnumeric(location)
         assert(isequal(size(location), [1 2]), ...
            'point location must be [lat lon]')
         [xq, yq] = projfwd(proj, location(1), location(2));
         location = [xq, yq];
      end
      [start, count, collapse, inslab, loctype] = ...
         icemodel.forcing.helpers.gridLocation(X, Y, location, ...
         kwargs.method, remap=kwargs.remap, validmask=validmask);
   end

   % Read each channel: per-day hyperslabs concatenated, stamped at the
   % MERRA bin centers, collapsed to the target (nearest cell, natural-
   % neighbour point, or polygon mean), then interpolated onto the hourly
   % axis.
   Time = hourlyAxis(years);
   Data = timetable(Time);
   for n = 1:size(channels, 1)
      [col, ncname, outname] = channels{n, :};
      [block, stamps] = readChannelSeries(inventory.(col), ncname, ...
         start, count);
      Data.(outname) = interp1(stamps, collapse(block), Time, ...
         'linear', 'extrap');
   end

   % Mass fluxes (PRECTOTCORR/PRECSNO/EVAP/RUNOFF) arrive already converted
   % from kg m-2 s-1 to mWE/h by icemodel.forcing.readMerra2.

   % Derived channels: wind speed/direction and relative humidity. The
   % derivable radiation terms (swu = swd - swn, netr = swn + lwn, lwu) are
   % intentionally NOT stored - icemodel.processmet recomputes them on load
   % from swd/albedo/tsfc/lwd. Only the native MERRA inputs (incl. the net
   % fluxes swn/lwn and the SNICEALB albedo) are carried.
   [Data.wspd, Data.wdir] = icemodel.forcing.helpers.windFromComponents( ...
      Data.uwind, Data.vwind);
   Data.rh = icemodel.vapor.relative_humidity_from_specific_humidity( ...
      Data.shum, Data.psfc, Data.tair);
   Data = removevars(Data, {'shum', 'uwind', 'vwind'});

   % Optional GEUS MODIS albedo channel (point locations).
   site_lat = icemodel.forcing.helpers.slabMean(LAT, start, count, inslab);
   site_lon = icemodel.forcing.helpers.slabMean(LON, start, count, inslab);
   if kwargs.modis_dir ~= ""
      Data.modis = modisChannel(kwargs.modis_dir, years, ...
         site_lat, site_lon, Data.Time);
   end

   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);

   % Userdata CustomProperties (MERRA carries no terrain height in
   % these collections; Elev is NaN).
   [site_x, site_y] = projfwd(proj, site_lat, site_lon);
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = site_x;
   Data.Properties.CustomProperties.Y = site_y;
   Data.Properties.CustomProperties.Lat = site_lat;
   Data.Properties.CustomProperties.Lon = site_lon;
   Data.Properties.CustomProperties.Elev = NaN;
   Data.Properties.CustomProperties.Slope = NaN;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];

   metadata = struct( ...
      'source_dir', source_dir, ...
      'collections', {collections'}, ...
      'n_files', numel(inventory.slv.files), ...
      'location_type', loctype, ...
      'method', kwargs.method, ...
      'remap', kwargs.remap, ...
      'grid_start', start, ...
      'grid_count', count, ...
      'n_cells', prod(count), ...
      'lat', site_lat, 'lon', site_lon, ...
      'humidity_kernel', ...
      "icemodel.vapor.relative_humidity_from_specific_humidity", ...
      'mass_flux_units', "mWE/h (rate)", ...
      'checks', checks);
end

%% Local functions
function t_hourly = hourlyAxis(years)
   %HOURLYAXIS Full on-the-hour axis covering the requested years.
   parts = cell(numel(years), 1);
   for n = 1:numel(years)
      t0 = datetime(years(n), 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
      parts{n} = (t0:hours(1):(t0 + calyears(1) - hours(1)))';
   end
   t_hourly = vertcat(parts{:});
end

function [block, stamps] = readChannelSeries(coll, ncname, start, count)
   %READCHANNELSERIES Concatenate one channel's hyperslab over the daily files.
   %
   % Returns the raw cells-by-time block (cells flattened column-major over
   % the hyperslab, matching gridLocation's collapse) plus the bin-center
   % timestamps: 24 hourly samples at :30 for tavg1 collections, 8
   % three-hourly samples at 1:30 for tavg3 (glc). The caller applies the
   % collapse (nearest / natural / polygon mean).
   % Per-file hyperslab read + standard-unit conversion + fill-masking is
   % delegated to the shared reader icemodel.forcing.readMerra2 (so mass
   % fluxes arrive already in mWE/h); this loop only concatenates the daily
   % files and stamps the bin centers.
   n_files = numel(coll.files);
   info = ncinfo(coll.files(1), ncname);
   n_per_day = info.Size(3);
   ncells = prod(count(1:2));

   block = nan(ncells, n_per_day * n_files);
   for k = 1:n_files
      block(:, (k-1)*n_per_day + 1:k*n_per_day) = ...
         icemodel.forcing.readMerra2(coll.files(k), ncname, ...
         start=start, count=count);
   end

   step = 24 / n_per_day;
   offsets = hours(step/2:step:24);
   stamps = reshape((coll.dates + offsets)', [], 1);
end

function mask = merraIceMask(glcfile, gridsize)
   %MERRAICEMASK Static land-ice mask from MERRA-2 glacier-tile validity.
   %
   % MERRA-2 tavg3_2d_glc variables (RUNOFF, SNOMAS_GL, ...) are defined
   % on the glacier tile; non-glacier cells carry the _FillValue. A cell
   % is land-ice where SNOMAS_GL is finite and below the fill magnitude.
   v = double(ncread(glcfile, 'SNOMAS_GL', [1 1 1], [Inf Inf 1]));
   mask = reshape(isfinite(v) & v < 1e14, gridsize);
end

function modis = modisChannel(modis_dir, years, lat, lon, Time)
   %MODISCHANNEL GEUS MODIS daily albedo interpolated to the time axis.
   modis = nan(numel(Time), 1);
   for yyyy = years
      match = dir(fullfile(modis_dir, sprintf('*_%d_*.nc', yyyy)));
      if numel(match) ~= 1
         error('icemodel:forcing:buildMerraData:modisNotFound', ...
            'expected one MODIS file for %d in %s, found %d', ...
            yyyy, modis_dir, numel(match))
      end
      [albedo, Tdaily] = icemodel.forcing.readGeusModis( ...
         string(fullfile(match.folder, match.name)), lat, lon);
      inyear = year(Time) == yyyy;
      modis(inyear) = icemodel.forcing.helpers.dailyToHourly( ...
         albedo, Tdaily, Time(inyear));
   end
end
