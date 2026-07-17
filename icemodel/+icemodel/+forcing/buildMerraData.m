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
   %    flx (hourly): HFLUX -> shf, EFLUX -> lhf [W m-2, sign flipped to
   %        the icemodel positive-toward-surface convention],
   %        PRECTOTCORR -> ppt, PRECSNO -> snowf [m s-1, the canonical
   %        water-equivalent precipitation rate], EVAP -> evap [mWE/h]
   %    glc (3-hourly): RUNOFF -> runoff [mWE/h], SNICEALB -> albedo,
   %        SNOWDP_GL -> snowd [m], SNOMAS_GL -> swe [kg m-2]
   %
   % Derived: wspd/wdir from U2M/V2M; rh [%] via the canonical icemodel.vapor
   % kernel. The derivable radiation terms (swu, lwu, netr) are NOT stored -
   % icemodel.processmet recomputes them on load from swd/albedo/tsfc/lwd.
   % Mass fluxes convert kg m-2 s-1 -> meters water equivalent per hour.
   % MERRA-2 posts tavg samples at the bin center (00:30 hourly, 01:30
   % 3-hourly). This application builder explicitly relabels only those
   % time-averaged collections to the INTERVAL START and holds each mean over
   % its declared support. The source reader preserves native coordinates.
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
   %  location - [lat lon] point, polyshape (EPSG:3413 m), or an Nx2 [lat lon]
   %             list of points. A point list returns a 1xN cell of Data
   %             timetables (metadata a 1xN struct array); the inventory, grid,
   %             and ice mask are read ONCE and every daily file opened ONCE,
   %             slicing each point's hyperslab from that single open. N=1 is
   %             the single-point path.
   %  years    - calendar years to extract
   %
   % Name-value
   %  source_dir : directory holding the collection subdirectories
   %      (flx/, glc/, rad/, slv/). Defaults to the gitignored cache
   %      data/forcing/merra2. Reference layout:
   %      /Volumes/S03/DATA/merra2/1hrly/ncfiles.
   %  modis_dir : optional GEUS MODIS albedo directory (see buildMarData)
   %  fillgaps  : opt in to legacy metchecks gap filling (default false)
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
      kwargs.fillgaps (1, 1) logical = false
      kwargs.method (1, 1) string {mustBeMember(kwargs.method, ...
         ["nearest", "natural"])} = "nearest"
      kwargs.remap (1, 1) string {mustBeMember(kwargs.remap, ...
         ["equal", "conservative"])} = "conservative"
   end

   source_dir = kwargs.source_dir;
   if source_dir == ""
      source_dir = string(getenv("ICEMODEL_MERRA_DIR"));
      if source_dir == ""
         source_dir = "/Volumes/S03/DATA/merra2/1hrly/ncfiles";
      end
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
   support_hours = struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3);

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

   % Read native collection coordinates once per collection, validate their
   % official center stamps, and convert only these declared tavg products to
   % application-level interval starts. Every channel in a collection shares
   % the resulting axis.
   interval_starts = struct();
   for c = collections'
      collection = c{1};
      interval_starts.(collection) = averagedIntervalStarts( ...
         inventory.(collection), support_hours.(collection));
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

   % Accept one location (1x2 point or polyshape, returns a single Data
   % timetable) OR a list of N points (Nx2 [lat lon], returns a 1xN cell of
   % Data timetables). N=1 is the single-point path. The file inventory, grid,
   % and ice mask above are read ONCE for the whole list; each point's grid
   % hyperslab + collapse rule is resolved here, then every channel's daily
   % files are opened ONCE and each point's hyperslab sliced from that open.
   [locations, batch] = ...
      icemodel.forcing.helpers.normalizeLocations(location);
   npts = numel(locations);

   grid = struct('X', X, 'Y', Y, 'LON', LON, 'LAT', LAT, 'proj', proj, ...
      'validmask', validmask, 'method', kwargs.method, 'remap', kwargs.remap);
   [slabs, collapses, sites] = deal(cell(1, npts), cell(1, npts), cell(1, npts));
   for p = 1:npts
      [slabs{p}, collapses{p}, sites{p}] = resolvePoint(grid, locations{p});
   end

   % Read each channel: daily hyperslabs concatenated per point (one open per
   % day per channel for ALL points), then collapsed and held over the source
   % interval support on the hourly application axis.
   Time = hourlyAxis(years);
   source_grid = merraTavg3SourceGrid(Time, interval_starts.glc);
   series = cell(size(channels, 1), npts);
   for n = 1:size(channels, 1)
      [col, ncname, ~] = channels{n, :};
      blocks = readChannelSeries(inventory.(col), ncname, slabs, ...
         interval_starts.(col));
      for p = 1:npts
         values = reshape(collapses{p}(blocks{p}), ...
            numel(interval_starts.(col)), []);
         series{n, p} = holdIntervalAverages(interval_starts.(col), ...
            values, Time, support_hours.(col));
      end
   end

   data_out = cell(1, npts);
   meta_out = cell(1, npts);
   for p = 1:npts
      Data = timetable(Time);
      for n = 1:size(channels, 1)
         Data.(channels{n, 3}) = series{n, p};
      end
      [data_out{p}, meta_out{p}] = finalizeMerraData(Data, sites{p}, ...
         slabs{p}, years, locations{p}, source_dir, collections, ...
         numel(inventory.slv.files), proj, support_hours, source_grid, kwargs);
   end

   metadata = [meta_out{:}];
   if batch
      Data = data_out;
   else
      Data = data_out{1};
   end
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

function blocks = readChannelSeries(coll, ncname, slabs, stamps)
   %READCHANNELSERIES Concatenate one channel's hyperslabs over the daily files.
   %
   % Returns a 1xN cell of raw cells-by-time blocks (one per point; cells
   % flattened column-major over each point's hyperslab, matching
   % gridLocation's collapse). The caller applies each point's collapse and
   % support-aware hold using the validated interval-start axis.
   % Per-file hyperslab read + standard-unit conversion + fill-masking is
   % delegated to the shared reader icemodel.forcing.readMerra2 (so mass
   % fluxes arrive already in mWE/h); this loop opens each daily file ONCE
   % (reading EVERY point's hyperslab from that open). A single point is just a
   % one-element slab list.
   n_files = numel(coll.files);
   npts = numel(slabs);
   info = ncinfo(coll.files(1), ncname);
   n_per_day = info.Size(3);
   assert(numel(stamps) == n_per_day * n_files, ...
      'MERRA interval-start axis does not match the channel record count')

   blocks = cell(1, npts);
   for p = 1:npts
      blocks{p} = nan(prod(slabs{p}(2, :)), n_per_day * n_files);
   end
   for k = 1:n_files
      day = icemodel.forcing.readMerra2(coll.files(k), ncname, slabs=slabs);
      cols = (k-1)*n_per_day + 1:k*n_per_day;
      for p = 1:npts
         blocks{p}(:, cols) = day{p};
      end
   end

end

function starts = averagedIntervalStarts(coll, support_hours)
   %AVERAGEDINTERVALSTARTS Validate native centers and relabel tavg support.
   % The official tavg1/tavg3 coordinates are centered at support/2 and repeat
   % at the support cadence. Decode every daily coordinate so a malformed middle
   % file or filename/native-day mismatch cannot acquire a synthetic proof.
   expected_offsets = hours((support_hours / 2):support_hours: ...
      (24 - support_hours / 2))';
   starts = NaT(numel(coll.files) * numel(expected_offsets), 1, ...
      'TimeZone', 'UTC');
   for k = 1:numel(coll.files)
      native = icemodel.forcing.helpers.readMerra2Time(coll.files(k));
      expected = coll.dates(k) + expected_offsets;
      if ~isequal(native, expected)
         error('icemodel:forcing:buildMerraData:badNativeTime', ...
            ['MERRA tavg%d file does not match its native interval-center ' ...
            'coordinate: %s'], support_hours, coll.files(k))
      end
      rows = (k - 1) * numel(expected_offsets) + (1:numel(expected_offsets));
      starts(rows) = native - hours(support_hours / 2);
   end
end

function output = holdIntervalAverages(starts, values, target, support_hours)
   %HOLDINTERVALAVERAGES Hold each mean over [start,start+support).
   assert(size(values, 1) == numel(starts), ...
      'MERRA value rows do not match the source time axis')
   assert(all(diff(starts) >= hours(support_hours)), ...
      'MERRA source intervals overlap or are out of order')

   % Missing source files create no assignment, so their hourly rows remain
   % NaN. ismember also handles nonconsecutive requested years without treating
   % the omitted years as a continuous target axis.
   output = nan(numel(target), size(values, 2));
   for offset = 0:support_hours - 1
      [present, target_rows] = ismember(starts + hours(offset), target);
      output(target_rows(present), :) = values(present, :);
   end
end

function proof = merraTavg3SourceGrid(target, starts)
   %MERRATAVG3SOURCEGRID Record the exact native glc stamps used by this build.
   expected = target(mod(hour(target), 3) == 0);
   if numel(unique(starts)) ~= numel(starts)
      error('icemodel:forcing:buildMerraData:duplicateGlcTime', ...
         'MERRA glc inventory contains duplicate native timestamps')
   end
   present = ismember(expected, starts);
   missing = expected(~present);
   proof = struct( ...
      'merra_tavg3_source_grid_policy', ...
      'native_glc_timestamp_inventory', ...
      'merra_tavg3_expected_source_row_count', numel(expected), ...
      'merra_tavg3_source_row_count', nnz(present), ...
      'merra_tavg3_source_time_gap_count', nnz(~present), ...
      'merra_tavg3_missing_source_times', missing);
end

function [slab, collapse, site] = resolvePoint(grid, location)
   %RESOLVEPOINT Map one point or polygon onto a MERRA grid hyperslab.
   %
   % Returns the hyperslab as a [start; count] 2x2 (for the batch reader),
   % the collapse handle, and the site lat/lon (slab mean). Conservative
   % polygon remap runs in MERRA's NATIVE geographic grid (regular lon/lat)
   % with exactremap UseGeoCoords=true, which computes true ellipsoidal
   % overlap areas - the correct conservative weighting for a lat/lon grid
   % (reprojecting to EPSG:3413 first would make the grid irregular).
   % Point/nearest and equal-weight stay in the projected grid.
   if isa(location, 'polyshape') && grid.remap == "conservative"
      [vlat, vlon] = projinv(grid.proj, location.Vertices(:, 1), ...
         location.Vertices(:, 2));
      [start, count, collapse, inslab, site.type] = ...
         icemodel.forcing.helpers.gridLocation(grid.LON, grid.LAT, ...
         polyshape(vlon, vlat), grid.method, remap="conservative", ...
         validmask=grid.validmask, usegeocoords=true);
   else
      if isnumeric(location)
         assert(isequal(size(location), [1 2]), ...
            'point location must be [lat lon]')
         [xq, yq] = projfwd(grid.proj, location(1), location(2));
         location = [xq, yq];
      end
      [start, count, collapse, inslab, site.type] = ...
         icemodel.forcing.helpers.gridLocation(grid.X, grid.Y, location, ...
         grid.method, remap=grid.remap, validmask=grid.validmask);
   end
   slab = [start; count];
   site.lat = icemodel.forcing.helpers.slabMean(grid.LAT, start, count, inslab);
   site.lon = icemodel.forcing.helpers.slabMean(grid.LON, start, count, inslab);
end

function [Data, metadata] = finalizeMerraData(Data, site, slab, years, ...
      location, source_dir, collections, n_files, proj, support_hours, ...
      source_grid, kwargs)
   %FINALIZEMERRADATA Post-process one point's assembled MERRA Data + metadata.
   % Identical to the legacy single-point tail: derived wind/RH, optional
   % MODIS channel, precip rate, metchecks, units, userdata CustomProperties,
   % and the provenance struct. start/count come from the point's slab.
   start = slab(1, :);
   count = slab(2, :);

   % Mass fluxes (PRECTOTCORR/PRECSNO/EVAP/RUNOFF) arrive already converted
   % from kg m-2 s-1 to mWE/h by icemodel.forcing.readMerra2.

   % MERRA-2 HFLUX/EFLUX are positive upward from the surface; icemodel uses
   % positive toward the surface for the staged forcing/evaluation contract.
   for ch = ["shf", "lhf"]
      Data.(ch) = -Data.(ch);
   end

   % Accumulate source-sign and time-support proof in the single metadata
   % record that is finalized and attached after all payload processing.
   artifact_metadata = Data.Properties.UserData;
   if isempty(artifact_metadata) || ~isstruct(artifact_metadata)
      artifact_metadata = struct();
   end
   artifact_metadata.merra_flux_sign_convention = ...
      'positive_toward_surface';
   proof_fields = fieldnames(source_grid);
   for k = 1:numel(proof_fields)
      field = proof_fields{k};
      artifact_metadata.(field) = source_grid.(field);
   end
   artifact_metadata.merra_source_time_coordinate = 'native_at_reader';
   artifact_metadata.merra_time_relabel_policy = ...
      'time_averaged_center_to_interval_start';
   artifact_metadata.merra_time_upsample_policy = ...
      'zero_order_hold_over_declared_support';
   artifact_metadata.merra_collection_support_hours = support_hours;
   [Data, artifact_metadata] = ...
      icemodel.forcing.helpers.applyMerraTimeSupport( ...
      Data, artifact_metadata);

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

   % Optional GEUS MODIS albedo channel at the requested location:
   % nearest/natural for a point, conservative (or equal) catchment mean
   % (area-weighted ROI mean) for a polygon.
   modis_metadata = struct();
   if kwargs.modis_dir ~= ""
      [modis, modis_metadata] = ...
         icemodel.forcing.helpers.modisAlbedoChannel( ...
         kwargs.modis_dir, years, location, kwargs.method, kwargs.remap, ...
         Data.Time);
      % Keep provenance current at the same source-selection boundary;
      % stageRcmForcing derives met from this Data without a second GEUS read.
      modis_fields = fieldnames(modis_metadata);
      for k = 1:numel(modis_fields)
         artifact_metadata.(modis_fields{k}) = ...
            modis_metadata.(modis_fields{k});
      end
   end

   % Precipitation to the canonical water-equivalent rate m s-1. MERRA mass
   % fluxes arrive as mWE/h, so dividing by 3600 s/h yields m s-1 for the
   % precipitation channels (ppt, snowf). The diagnostic mass fluxes
   % (evap/runoff) keep their natural mWE/h rate; swe is a store (kg m-2).
   for ch = ["ppt", "snowf"]
      Data.(ch) = Data.(ch) / 3600;
   end

   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);
   % Attach optional MODIS only after generic gap filling so missing source years
   % remain missing even for direct builder calls that request filled met inputs.
   if isfield(modis_metadata, 'modis_coverage_years') ...
         && ~isempty(modis_metadata.modis_coverage_years)
      Data.modis = modis;
   end

   % Per-variable units from the shared canonical map. Precipitation is m s-1
   % (converted above); the diagnostic mass fluxes are mWE/h rates and swe is
   % a kg m-2 store.
   Data.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(Data.Properties.VariableNames));

   % Attach the shared location schema. These MERRA collections carry no
   % terrain height or surface slope, so both source fields remain NaN.
   [site_x, site_y] = projfwd(proj, site.lat, site.lon);
   location_metadata = struct( ...
      'lat_wgs84', site.lat, 'lon_wgs84', site.lon, ...
      'x_epsg3413', site_x, 'y_epsg3413', site_y, 'elev_m', NaN);
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, location_metadata);

   % Finalize the exact record returned to the caller and persisted on Data.
   % Starting from the accumulated source proof avoids parallel public and
   % payload metadata records that can silently drift apart.
   metadata = artifact_metadata;
   metadata.source_dir = source_dir;
   metadata.collections = collections';
   metadata.n_files = n_files;
   metadata.location_type = site.type;
   metadata.method = kwargs.method;
   metadata.remap = kwargs.remap;
   metadata.grid_start = start;
   metadata.grid_count = count;
   metadata.n_cells = prod(count);
   metadata.lat = site.lat;
   metadata.lon = site.lon;
   metadata.humidity_kernel = ...
      "icemodel.vapor.relative_humidity_from_specific_humidity";
   metadata.mass_flux_units = ...
      "precip m s-1; diagnostic fluxes mWE/h (rate)";
   metadata.merra_collection_support_hours = support_hours;
   metadata.checks = checks;
   metadata = icemodel.forcing.helpers.columnizeMetadata(metadata);
   Data.Properties.UserData = metadata;
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
