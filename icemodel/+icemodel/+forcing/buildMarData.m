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
   %  water-equivalent precipitation rates], melt [mWE/h], runoff and smb
   %  [mWE/h, hourly RUH/SMBH where their UTC-day sums match native daily
   %  RU/SMB, otherwise native daily RU/SMB divided by 24], subl [mWE/h,
   %  native hourly SUH when available];
   %  daily: snowd [m], cfrac [-], tsfc [K], psfc [Pa]; derived: wspd
   %  [m s-1], wdir [deg] (from UUH/VVH), rh [%] (icemodel.vapor kernel);
   %  optional diagnostics: subl_evap [mWE/h, native daily SU/24],
   %  refreeze_deposition [mWE/h, native daily RZ/24], and daily ME used only
   %  to validate hourly MEH; modis [-] (GEUS MODIS daily albedo).
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
   %  metadata - provenance: files read, cell/polygon info, MAR native-daily
   %             diagnostic policy, and build checks
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
      'snowf', 'SFH' ; 'rainf', 'RFH' ; 'melt', 'MEH'
      'runoff', 'RUH'; 'shf',  'SHFH' ; 'lhf',  'LHFH'; 'smb', 'SMBH'
      };
   daily_vars = {
      'snowd', 'SHSN2', "none", false; 'cfrac', 'CC', "none", false
      'tsfc', 'ST', "none", false; 'psfc', 'SP', "none", false
      };

   % Accept one location (1x2 point or polyshape, returns a single Data
   % timetable) OR a list of N points (Nx2 [lat lon], returns a 1xN cell
   % array of Data timetables). N=1 is just the single-point path, so there
   % is ONE code path: resolve every point's grid hyperslab up front (grid
   % metadata + interpolants read ONCE), then loop years opening each yearly
   % file ONCE and reading every point's hyperslab from that single open.
   [locations, batch] = ...
      icemodel.forcing.helpers.normalizeLocations(location);
   npts = numel(locations);

   files = strings(numel(years), 1);
   for n = 1:numel(years)
      files(n) = locateMarFile(source_dir, years(n));
   end

   % Read each yearly header once, then reuse the compact schema matrix for all
   % required and optional availability decisions. Full MAR files are large;
   % repeated ncinfo calls would add avoidable staging overhead even though no
   % field data are read here.
   schema_names = ["RU", "SMB", "SUH", "SU", "RZ", "ME"];
   schema_available = false(numel(files), numel(schema_names));
   for n = 1:numel(files)
      info = ncinfo(files(n));
      schema_available(n, :) = ismember( ...
         schema_names, string({info.Variables.Name}));
   end
   has_all = @(variables) all(schema_available(:, ...
      ismember(schema_names, variables)), 'all');
   if has_all(["RU", "SMB"])
      daily_vars = [daily_vars; ...
         {'runoff_daily', 'RU', "site", true; ...
         'smb_daily', 'SMB', "site", true}];
   end

   % Optional MAR mass-balance diagnostics are emitted only when every selected
   % source year declares the native field. This keeps mixed/reduced archives
   % explicit and avoids filling a missing source year with invented support.
   if has_all("SUH")
      hourly_vars = [hourly_vars; {'subl', 'SUH'}];
   end
   if has_all("SU")
      daily_vars = [daily_vars; {'subl_evap', 'SU', "site", true}];
   end
   if has_all("RZ")
      daily_vars = [daily_vars; ...
         {'refreeze_deposition', 'RZ', "ice", true}];
   end
   if has_all("ME")
      daily_vars = [daily_vars; ...
         {'melt_daily_reference', 'ME', "ice", true}];
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
      sectors = cellfun(@(site) site.sector, sites);
      blocks = extractOneYear( ...
         files(n), hourly_vars, daily_vars, slabs, sectors);
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
      [~, nearest] = min(hypot(grid.Xnat(:) - location(1), ...
         grid.Ynat(:) - location(2)));
      site.sector = 1 + double(grid.srf(nearest) ~= 4);
   elseif isa(location, 'polyshape')
      [vlat, vlon] = projinv(grid.proj, location.Vertices(:, 1), ...
         location.Vertices(:, 2));
      location = polyshape(grid.toNativeX(vlon, vlat), ...
         grid.toNativeY(vlon, vlat));
      % Polygon remapping already uses the ice-sheet valid mask, so its mass
      % diagnostics must come from MAR's permanent-ice surface sector.
      site.sector = 1;
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

function blocks = extractOneYear( ...
      filename, hourly_vars, daily_vars, slabs, sectors)
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
      sector_mode = string(daily_vars{n, 3});
      if sector_mode ~= "none"
         % SU follows the target surface sector; singleton RZ/ME products are
         % explicitly read from their permanent-ice sector.
         selected_sectors = sectors;
         if sector_mode == "ice"
            selected_sectors = ones(size(sectors));
         end
         [data, units, Tdaily] = icemodel.forcing.readMar3p11(filename, ...
            daily_vars{n, 2}, slabs=slabs, sector=selected_sectors);
      else
         [data, units, Tdaily] = icemodel.forcing.readMar3p11(filename, ...
            daily_vars{n, 2}, slabs=slabs);
      end
      if daily_vars{n, 4}
         assert(string(units) == "mWE/day", ...
            'MAR daily mass diagnostics must use mWE/day after conversion')
      end
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
      daily = collapse(blocks.daily{n}{p});
      if daily_vars{n, 4}
         % Daily mass rates use exact UTC-day support. RU/SMB remain beside
         % their hourly products for selective QC; SU/RZ become explicitly
         % combined diagnostics, while ME is a private MEH validation channel.
         part.(daily_vars{n, 1}) = ...
            icemodel.forcing.helpers.dailyToHourly( ...
            daily / 24, blocks.Tdaily, part.Time, method="previous");
      else
         if string(daily_vars{n, 1}) == "cfrac"
            % Cloud cover is a bounded diagnostic. Endpoint holding plus an
            % explicit range postcondition prevents year-end extrapolation.
            part.(daily_vars{n, 1}) = ...
               icemodel.forcing.helpers.dailyToHourly( ...
               daily, blocks.Tdaily, part.Time, bounds=[0 1]);
         else
            part.(daily_vars{n, 1}) = ...
               icemodel.forcing.helpers.dailyToHourly( ...
               daily, blocks.Tdaily, part.Time);
         end
      end
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
   modis_metadata = struct();
   if kwargs.modis_dir ~= ""
      [modis, modis_metadata] = ...
         icemodel.forcing.helpers.modisAlbedoChannel( ...
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
   % derives ppt = rainf + snowf). Melt keeps its native mWE/h rate.
   for ch = ["snowf", "rainf"]
      Data.(ch) = Data.(ch) / 3600;
   end

   % Native daily delayed RU and daily SMB constrain the hourly diagnostics.
   % Retain them outside metchecks so missing source days cannot be synthesized
   % by generic gap filling; the selective helper below preserves raw hourly
   % structure only where its complete UTC-day aggregate is source-consistent.
   names = string(Data.Properties.VariableNames);
   if all(ismember(["runoff_daily", "smb_daily"], names))
      replacements = struct('runoff', Data.runoff_daily, ...
         'smb', Data.smb_daily);
      Data = removevars(Data, {'runoff_daily', 'smb_daily'});
   else
      % Some intentionally reduced test/legacy sources contain only hourly
      % RUH/SMBH. Preserve those source values but mark native-daily QC as not
      % applicable; production MAR archives carry RU/SMB and take the branch
      % above. This is a source-schema compatibility hook, not a silent mask.
      replacements = struct();
   end
   if ismember("melt_daily_reference", names)
      melt_daily_rate = Data.melt_daily_reference;
      Data = removevars(Data, 'melt_daily_reference');
   else
      % Reduced/legacy sources without daily ME retain hourly MEH and receive
      % explicit not-available validation provenance below.
      melt_daily_rate = zeros(0, 1);
   end
   mass_names = intersect(["runoff", "smb"], ...
      string(Data.Properties.VariableNames), 'stable');
   diagnostic_names = intersect( ...
      ["melt", "subl", "subl_evap", "refreeze_deposition"], ...
      string(Data.Properties.VariableNames), 'stable');
   preserved_names = [mass_names, diagnostic_names];
   native_mass_diagnostics = Data(:, preserved_names);
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);

   % Generic met gap filling must not turn missing native mass samples into
   % apparent source support. Restore RUH/SMBH, native hourly MEH-derived melt,
   % and every optional diagnostic before source-aware QC/provenance classifies
   % their temporal support.
   for channel = reshape(preserved_names, 1, [])
      Data.(channel) = native_mass_diagnostics.(channel);
   end
   % MODIS is optional source data: attach it after generic metchecks so a
   % direct builder call with gap filling cannot invent values in missing years.
   if isfield(modis_metadata, 'modis_coverage_years') ...
         && ~isempty(modis_metadata.modis_coverage_years)
      Data.modis = modis;
   end

   % Apply source-specific QC after generic met checking so neither the SHSN2
   % discontinuity mask nor a missing daily mass reference is gap-filled.
   [Data, snow_metadata] = ...
      icemodel.forcing.helpers.applyMarSnowDepthQualityControl(Data);
   [Data, qc_metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      Data, replacements, sector=site.sector);
   snow_fields = fieldnames(snow_metadata);
   for k = 1:numel(snow_fields)
      qc_metadata.(snow_fields{k}) = snow_metadata.(snow_fields{k});
   end
   qc_metadata = icemodel.forcing.helpers.marDiagnosticMetadata( ...
      Data, melt_daily_rate, qc_metadata, sector=site.sector);

   % Per-variable units from the shared canonical map. Precipitation is m s-1
   % (converted above); the diagnostic mass fluxes are mWE/h rates (cumulative
   % sums need the timestep in hours).
   Data.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(Data.Properties.VariableNames));

   % Attach the shared location schema while retaining MAR's native terrain
   % elevation, projected coordinates, and real surface slope.
   location_metadata = struct( ...
      'lat_wgs84', site.lat, 'lon_wgs84', site.lon, ...
      'x_epsg3413', site.x, 'y_epsg3413', site.y, ...
      'elev_m', site.elev, 'slope', site.slope);
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, location_metadata);

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
   % The channel helper already resolved exact source years while reading the
   % physical data, so copying its canonical fields adds no source access.
   modis_fields = fieldnames(modis_metadata);
   for k = 1:numel(modis_fields)
      metadata.(modis_fields{k}) = modis_metadata.(modis_fields{k});
   end
   qc_fields = fieldnames(qc_metadata);
   for k = 1:numel(qc_fields)
      metadata.(qc_fields{k}) = qc_metadata.(qc_fields{k});
   end
   metadata = icemodel.forcing.helpers.columnizeMetadata(metadata);
   Data.Properties.UserData = metadata;
   if site.srf_warning
      warning('icemodel:forcing:buildMarData:surfaceNotIce', ...
         'MAR surface type at the requested location is not ice sheet')
   end
end
