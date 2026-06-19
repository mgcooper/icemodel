function [Data, metadata] = buildRacmoData(location, years, kwargs)
   %BUILDRACMODATA Build a Data timetable from RACMO2.3 NetCDF files.
   %
   %  [Data, metadata] = icemodel.forcing.buildRacmoData(location, years)
   %  [Data, metadata] = ... buildRacmoData(_, source_dir=..., dt="3hr")
   %
   % Extracts the RACMO2.3p3 surface energy/mass-balance channels at a
   % point or averaged over a polygon. The RACMO archive is organized
   % as one multi-year 3-hourly NetCDF per variable
   % (<var>.RACMO23p3_*_FGRN11_*.3H*.nc on the FGRN11 rotated-pole
   % grid); the available variables are SMB components and surface
   % fluxes, NOT full meteorological forcing (no air temperature, wind,
   % humidity, or pressure), so RACMO Data files serve evaluation and
   % met-swap of flux channels rather than met-file creation.
   %
   % Channels (file prefix -> output, standard units):
   %    swsd -> swd, lwsd -> lwd, swsn -> swn, lwsn -> lwn   [W m-2]
   %    senf -> shf, latf -> lhf                             [W m-2]
   %    precip -> precip, snowmelt -> melt, runoff -> runoff,
   %    smb -> smb, refreeze -> refreeze, subl -> subl,
   %    sndiv -> sndiv, meltin -> meltin                     [mWE/h]
   % Derived:  albedo [-] = 1 - swn/swd (RACMO ships no albedo variable).
   % Optional: modis [-] (GEUS MODIS daily albedo, when modis_dir is given).
   %
   % Derivable radiation terms (swu, lwu, netr) are NOT stored: icemodel.
   % processmet recomputes them on load from swd/albedo/tsfc/lwd. Only the
   % native RACMO inputs (incl. the net fluxes swn/lwn) and the non-derivable
   % albedo are carried.
   %
   % Mass fluxes convert from kg m-2 s-1 to meters water equivalent per
   % hour (x 3600 / 1000); they are rates, so cumulative sums must
   % multiply by the timestep in hours (1 for dt="1hr", 3 for "3hr").
   %
   % Legacy: reimplements runoff/functions/saveRacmoData.m (the original
   % retained, unchanged, as the legacy reference workflow).
   %
   % Inputs
   %  location - [lat lon] point or polyshape in EPSG:3413 meters
   %  years    - calendar years to keep (subset of the archive span)
   %
   % Name-value
   %  source_dir : directory with the per-variable RACMO files. Defaults
   %      to the gitignored cache data/forcing/racmo. Reference layout:
   %      /Volumes/S03/DATA/greenland/racmo2p3/surface.
   %  modis_dir : directory with GEUS Greenland_Reflectivity_<YYYY>_5km_C6.nc
   %      files; when given, adds a daily MODIS albedo channel at the site.
   %      Reference layout: /Volumes/S03/DATA/greenland/geus/albedo/gris.
   %  dt : "1hr" (default; linear interpolation to hourly, the legacy
   %      behavior) or "3hr" (native posting)
   %  method : point sampling "nearest" (default) | "natural"
   %  remap : polygon aggregation "conservative" (default) | "equal"
   %      ("conservative" uses exactremap with FGRN11 gridarea as the true
   %      cell areas and the ice mask as the valid-cells mask)
   %
   % Outputs
   %  Data     - timetable with userdata CustomProperties (X, Y, Lat,
   %             Lon, Elev, Slope, ScalarUnits)
   %  metadata - provenance: files read, grid hyperslab, cell count
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.helpers.gridLocation,
   %  icemodel.forcing.helpers.writeuserdata

   arguments
      location
      years (1, :) double {mustBeInteger}
      kwargs.source_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.method (1, 1) string {mustBeMember(kwargs.method, ...
         ["nearest", "natural"])} = "nearest"
      kwargs.remap (1, 1) string {mustBeMember(kwargs.remap, ...
         ["equal", "conservative"])} = "conservative"
      kwargs.dt (1, 1) string {mustBeMember(kwargs.dt, ...
         ["1hr", "3hr"])} = "1hr"
   end

   source_dir = kwargs.source_dir;
   if source_dir == ""
      source_dir = string(fullfile(icemodel.getpath('data'), ...
         'forcing', 'racmo'));
   end
   if ~isfolder(source_dir)
      error('icemodel:forcing:buildRacmoData:sourceNotFound', ...
         ['RACMO source directory not found: %s. Pass source_dir or ' ...
         'stage the per-variable files (reference layout: ' ...
         '/Volumes/S03/DATA/greenland/racmo2p3/surface).'], source_dir)
   end

   % Channel table: file prefix -> output name.
   channels = {
      'swsd', 'swd'     ; 'lwsd', 'lwd'      ; 'swsn', 'swn'
      'lwsn', 'lwn'     ; 'senf', 'shf'      ; 'latf', 'lhf'
      'precip', 'precip'; 'snowmelt', 'melt' ; 'runoff', 'runoff'
      'smb', 'smb'      ; 'refreeze', 'refreeze'; 'subl', 'subl'
      'sndiv', 'sndiv'  ; 'meltin', 'meltin'
      };

   % Resolve the available files and the grid from the first one.
   [files, found] = locateRacmoFiles(source_dir, channels(:, 1));
   first = files(find(found, 1));
   assert(~isempty(first), 'no RACMO variable files found in %s', source_dir)

   LAT = double(ncread(first, 'lat'));
   LON = double(ncread(first, 'lon'));
   proj = icemodel.forcing.helpers.psnProjection();
   [X, Y] = projfwd(proj, LAT, LON);

   % Conservative polygon remap runs in RACMO's NATIVE rotated-pole frame
   % (the FGRN11 rlon/rlat grid is regular there; reprojecting to EPSG:3413
   % is curvilinear). exactremap's rotated-pole support handles the rotation
   % from the CF grid mapping and weights cells by the shipped true cell
   % areas (gridarea); off-ice cells (IceMask) are inpainted. Point/nearest
   % and equal-weight stay in the projected grid.
   if isa(location, 'polyshape') && kwargs.remap == "conservative"
      [start, count, collapse, inslab, loctype] = ...
         resolveRacmoConservative(first, source_dir, location, proj);
   else
      if isnumeric(location)
         assert(isequal(size(location), [1 2]), ...
            'point location must be [lat lon]')
         [xq, yq] = projfwd(proj, location(1), location(2));
         location = [xq, yq];
      end
      [start, count, collapse, inslab, loctype] = ...
         icemodel.forcing.helpers.gridLocation(X, Y, location, kwargs.method, ...
         remap=kwargs.remap);
   end

   % Time axis (shared by all variables): days since 1950-01-01, native
   % 3-hourly posting.
   t = ncread(first, 'time');
   t_units = ncreadatt(first, 'time', 'units');
   assert(startsWith(t_units, 'days since 1950-01-01'), ...
      'unexpected RACMO time units: %s', t_units)
   Time = datetime(1950, 1, 1, 'TimeZone', 'UTC') + days(double(t));
   keep = ismember(year(Time), years);
   assert(any(keep), ...
      'requested years %s not in the RACMO archive span %d-%d', ...
      mat2str(years), year(Time(1)), year(Time(end)))
   Time = Time(keep);

   % Read each available channel at the hyperslab and collapse to the
   % target (single cell or weighted polygon average).
   Data = timetable(Time);
   units = strings(1, size(channels, 1));
   for n = 1:size(channels, 1)
      if ~found(n)
         continue
      end
      [series, units(n)] = readChannel(files(n), channels{n, 1}, ...
         start, count, collapse, keep);
      Data.(channels{n, 2}) = series;
   end
   units = units(found);

   % Interpolate to hourly (legacy behavior) unless native requested.
   % The full-year hourly axis extends past the last 3-hourly posting
   % (21:00 on Dec 31); the trailing hours extrapolate linearly.
   if kwargs.dt == "1hr"
      t1 = dateshift(Time(1), 'start', 'year');
      t2 = dateshift(Time(end), 'start', 'year') + calyears(1) - hours(1);
      t_hourly = (t1:hours(1):t2)';
      t_hourly = t_hourly(ismember(year(t_hourly), years));
      Data = retime(Data, t_hourly, 'linear', 'EndValues', 'extrap');
   end

   % Site location (also needed for the optional MODIS channel below).
   site_lat = icemodel.forcing.helpers.slabMean(LAT, start, count, inslab);
   site_lon = icemodel.forcing.helpers.slabMean(LON, start, count, inslab);
   [site_x, site_y] = projfwd(proj, site_lat, site_lon);

   % Derived surface albedo. RACMO ships no albedo variable, so recover it
   % from downwelling and net shortwave: albedo = SWup/SWdown =
   % (SWdown - SWnet)/SWdown = 1 - swn/swd. The instantaneous ratio is
   % unreliable at low sun (small, noisy denominator near dawn/dusk and
   % through the polar night, where swd -> 0), so it is computed only where
   % SWdown exceeds a low-insolation floor and left NaN otherwise;
   % metchecks then linearly fills those gaps and clamps to [0.05, 0.98].
   % This refines the legacy saveRacmoData method, which formed the ratio
   % everywhere and relied on clamping alone. A daytime-only ratio is the
   % right hourly estimate here: where it is NaN, swd ~ 0, so the
   % processmet reconstruction swn = swd*(1-albedo) is insensitive to the
   % filled albedo. (Albedo is a forcing INPUT to icemodel.processmet, not
   % one of the radiation terms it derives, so it must be carried here.)
   swdown_floor = 10;   % [W m-2] below this, 1 - swn/swd is low-sun noise
   albedo = 1 - Data.swn ./ Data.swd;
   albedo(~(Data.swd >= swdown_floor)) = NaN;   % also catches swd == 0 / NaN
   Data.albedo = albedo;
   units(end+1) = "-";

   % Optional GEUS MODIS daily albedo at the site (nearest cell), matching
   % the legacy saveRacmoData MODIS channel.
   if kwargs.modis_dir ~= ""
      Data.modis = modisChannel(kwargs.modis_dir, years, site_lat, ...
         site_lon, Data.Time);
      units(end+1) = "reflectivity";
   end

   Data.Properties.VariableUnits = cellstr(units);

   % QA/QC: gap-fill + physical clamps. The legacy saveRacmoData ran
   % metchecks as its final step; here it also clamps the derived albedo to
   % [0.05, 0.98] and linearly fills the polar-night gap.
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data);

   % Userdata CustomProperties. Elevation from the topography channel
   % when present in the directory; otherwise NaN.
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = site_x;
   Data.Properties.CustomProperties.Y = site_y;
   Data.Properties.CustomProperties.Lat = site_lat;
   Data.Properties.CustomProperties.Lon = site_lon;
   Data.Properties.CustomProperties.Elev = ...
      readElevation(source_dir, start, count, inslab);
   Data.Properties.CustomProperties.Slope = NaN;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];

   metadata = struct( ...
      'source_files', files(found), ...
      'location_type', loctype, ...
      'method', kwargs.method, ...
      'remap', kwargs.remap, ...
      'grid_start', start, ...
      'grid_count', count, ...
      'n_cells', prod(count), ...
      'lat', site_lat, 'lon', site_lon, ...
      'dt', kwargs.dt, ...
      'mass_flux_units', "mWE/h (rate; cumulative sums need dt hours)", ...
      'checks', checks);
end

%% Local functions
function [files, found] = locateRacmoFiles(source_dir, prefixes)
   %LOCATERACMOFILES Resolve one file per RACMO variable prefix.
   n = numel(prefixes);
   files = strings(n, 1);
   found = false(n, 1);
   for k = 1:n
      match = dir(fullfile(source_dir, [prefixes{k} '.RACMO*.nc']));
      if isscalar(match)
         files(k) = string(fullfile(match.folder, match.name));
         found(k) = true;
      end
   end
end

function [series, unit] = readChannel(filename, prefix, start, count, ...
      collapse, keep)
   %READCHANNEL Read one RACMO variable (standard units), collapse, subset.
   % The hyperslab read + unit conversion is delegated to the shared reader
   % icemodel.forcing.readRacmo2p3; COLLAPSE (from gridLocation) then reduces
   % the cells-by-time block to the target series (nearest cell,
   % natural-neighbour point, or polygon mean), and KEEP subsets to the
   % requested years.
   [data, unit] = icemodel.forcing.readRacmo2p3(filename, prefix, ...
      start=start, count=count);
   series = collapse(data);
   series = series(keep);
   unit = string(unit);
end

function elev = readElevation(source_dir, start, count, inslab)
   %READELEVATION Mean terrain height from the FGRN11 topography file.
   %
   % The per-variable files carry only the model level height; terrain
   % comes from the companion FGRN11-topography file, expected in the
   % source directory or its parent.
   match = [dir(fullfile(source_dir, '*topography*.nc')); ...
      dir(fullfile(source_dir, '..', '*topography*.nc'))];
   if isempty(match)
      elev = NaN;
      return
   end
   topo = double(ncread( ...
      fullfile(match(1).folder, match(1).name), 'Topography'));
   elev = icemodel.forcing.helpers.slabMean(topo, start, count, inslab);
end

function [cellareas, validmask] = racmoGridStatics(source_dir)
   %RACMOGRIDSTATICS True cell areas [m^2] + ice mask from the topo file.
   %
   % Returns full grids aligned cell-for-cell with the data grid (they
   % share the FGRN11 rlon/rlat grid). gridarea is stored in km^2 (per the
   % FGRN11 topography file; the units attribute is uninformative), so it
   % is converted to m^2. Returns [] when the topography file is absent.
   match = [dir(fullfile(source_dir, '*topography*.nc')); ...
      dir(fullfile(source_dir, '..', '*topography*.nc'))];
   if isempty(match)
      cellareas = [];
      validmask = [];
      return
   end
   filename = fullfile(match(1).folder, match(1).name);
   cellareas = double(ncread(filename, 'gridarea')) * 1e6;   % km^2 -> m^2
   try
      validmask = logical(round(double(ncread(filename, 'IceMask'))));
   catch
      validmask = double(ncread(filename, 'Promicemask')) > 0;
   end
end

function [start, count, collapse, inslab, loctype] = ...
      resolveRacmoConservative(filename, source_dir, P, proj)
   %RESOLVERACMOCONSERVATIVE Rotated-pole conservative remap onto a polygon.
   %
   % RACMO FGRN11 is a CF rotated_latitude_longitude grid: regular in the
   % rotated rlon/rlat frame. The catchment polygon (EPSG:3413) is mapped to
   % true geographic coordinates for exactremap (which re-rotates it via the
   % grid mapping) and to rotated coordinates to bound the read hyperslab.
   % Weights are area-weighted by the shipped true cell areas (gridarea) with
   % off-ice cells (IceMask) inpainted. We use exactremap's 'weights' mode and
   % apply the weights ourselves because its rotated-pole 'areaavg' path is
   % currently broken (returns ~0; exactremap-0hv); the weights path is exact
   % (its rotated-pole geometry test passes to 1e-5) and matches the legacy
   % runoff approach of computing static remap weights once per polygon.
   loctype = "polygon";

   rlon = double(ncread(filename, 'rlon'));
   rlat = double(ncread(filename, 'rlat'));
   gm = racmoGridMapping(filename);
   [cellareas, validmask] = racmoGridStatics(source_dir);
   assert(~isempty(cellareas), ['conservative RACMO remap needs the FGRN11 ' ...
      'topography file (gridarea); none found near %s'], source_dir)

   % Polygon -> true geo (for exactremap) and rotated (to bound the slab).
   % Preserve NaN row separators so a multi-region polyshape (e.g. a catchment
   % with holes / disjoint parts) keeps its structure for exactremap; transform
   % only the finite vertices and scatter back, leaving NaN rows in place.
   vx = P.Vertices(:, 1);
   vy = P.Vertices(:, 2);
   fin = isfinite(vx) & isfinite(vy);
   vlat = nan(size(vx));
   vlon = nan(size(vx));
   [vlat(fin), vlon(fin)] = projinv(proj, vx(fin), vy(fin));
   Pgeo = [vlon, vlat];
   [vrlat, vrlon] = geo2rotated(vlat(fin), vlon(fin), ...
      gm.grid_north_pole_latitude, gm.grid_north_pole_longitude);

   % Slab covering the rotated bounding box plus a 2-cell pad.
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

   % Conservative weights over the slab. exactremap takes the gridvector axes
   % (rlon, rlat) in MESHGRID convention: its grid, the 2-D CellAreas/mask it
   % consumes, and the W it returns are all laid out [numel(rlat) x numel(rlon)]
   % (rlat down rows). The readRacmo2p3 data block, the slab, and slabMean are
   % NDGRID [numel(rlon) x numel(rlat)] (rlon down rows). So transpose the
   % cell-area/mask slabs going in, and reorient W back to ndgrid coming out,
   % so the weights align with the data block and slab indices.
   W = exactremap([], rlon(rows), rlat(cols), Pgeo, 'weights', ...
      'GridMapping', gm, 'CellAreas', cellareas(rows, cols).', ...
      'ValidCellsMask', validmask(rows, cols).', 'InfillMasked', true);
   w = reshape(W.W(:), [count(2), count(1)]).';   % meshgrid -> ndgrid
   w = w(:);
   w(~isfinite(w)) = 0;
   assert(sum(w) > 0, 'polygon does not overlap the RACMO grid')
   wn = w / sum(w);
   collapse = @(block) (wn.' * block).';

   inslab = find(w > 0);
   if isempty(inslab); inslab = 1; end
end

function gm = racmoGridMapping(filename)
   %RACMOGRIDMAPPING CF rotated-pole grid mapping struct from rotated_pole.
   gm = struct( ...
      'grid_mapping_name', 'rotated_latitude_longitude', ...
      'grid_north_pole_latitude', ...
      ncreadatt(filename, 'rotated_pole', 'grid_north_pole_latitude'), ...
      'grid_north_pole_longitude', ...
      ncreadatt(filename, 'rotated_pole', 'grid_north_pole_longitude'));
end

function modis = modisChannel(modis_dir, years, lat, lon, Time)
   %MODISCHANNEL GEUS MODIS daily albedo at the site, interpolated hourly.
   modis = nan(numel(Time), 1);
   for yyyy = years
      match = dir(fullfile(modis_dir, sprintf('*_%d_*.nc', yyyy)));
      if numel(match) ~= 1
         error('icemodel:forcing:buildRacmoData:modisNotFound', ...
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
