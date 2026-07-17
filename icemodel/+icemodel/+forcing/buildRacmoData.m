function [Data, metadata] = buildRacmoData(location, years, kwargs)
   %BUILDRACMODATA Build a Data timetable from RACMO2.3 NetCDF files.
   %
   %  [Data, metadata] = icemodel.forcing.buildRacmoData(location, years)
   %  [Data, metadata] = ... buildRacmoData(_, source_dir=..., dt="3hr")
   %
   % Extracts the RACMO2.3p3 surface energy/mass-balance channels at a
   % point or averaged over a polygon. The RACMO archive is organized
   % as one multi-year 3-hourly NetCDF per variable
   % (<var>.RACMO23p3_*_FGRN11_*.3H*.nc on the FGRN11 rotated-pole grid).
   % The AVAILABLE source files carry radiation (swsd/lwsd -> swd/lwd, with a
   % derived albedo = 1 - swn/swd), turbulent fluxes (senf/latf -> shf/lhf),
   % precip and the SMB components - but LACK the near-surface meteorological
   % STATE variables tair, wspd, rh, and psfc. So RACMO Data files serve
   % evaluation/reference and met-swap of the radiation/flux channels rather
   % than standalone met-file creation. (This is a property of the available
   % 2.3p3 files, not RACMO in general - the full state set must be obtained
   % from the RACMO developers, or borrowed from MAR/MERRA/PROMICE at the point,
   % for a RACMO-forced run.)
   %
   % Channels (file prefix -> output, standard units):
   %    swsd -> swd, lwsd -> lwd, swsn -> swn, lwsn -> lwn   [W m-2]
   %    senf -> shf, latf -> lhf                             [W m-2]
   %    precip -> ppt      [m s-1, the canonical water-equivalent precip rate]
   %    snowmelt -> melt, runoff -> runoff, smb -> smb,
   %    refreeze -> refreeze, -subl -> subl, sndiv -> sndiv,
   %    meltin -> meltin                                     [mWE/h]
   % RACMO stores sublimation as negative mass change and deposition as
   % positive. Canonical subl is positive surface loss, matching MAR SUH.
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
   % multiply by the timestep in hours (1 for dt="1hr", 3 for "3hr"). The
   % precipitation channel (ppt) is then carried in the canonical m s-1
   % water-equivalent rate; the diagnostic mass fluxes keep mWE/h.
   %
   % Legacy: reimplements runoff/functions/saveRacmoData.m (the original
   % retained, unchanged, as the legacy reference workflow).
   %
   % Inputs
   %  location - [lat lon] point, polyshape (EPSG:3413 m), or an Nx2 [lat lon]
   %             list of points. A point list returns a 1xN cell of Data
   %             timetables (metadata a 1xN struct array); the per-variable
   %             files and grid are read ONCE and every variable file opened
   %             ONCE, slicing each point's hyperslab from that single open
   %             (decisive for the multi-GB subsurface files). N=1 is the
   %             single-point path.
   %  years    - calendar years to keep (subset of the archive span)
   %
   % Name-value
   %  source_dir : directory with the per-variable RACMO files. Resolves
   %      ICEMODEL_RACMO_DIR, else the reference layout
   %      /Volumes/S03/DATA/greenland/racmo2p3/subsurface (the full RACMO2.3p3
   %      FGRN11 run, 2012-2018 - preferred for firn work over the shorter
   %      2012-2015 .../surface "no_subsurf_en" product). The per-variable
   %      filename pattern (<var>.RACMO*.nc) matches either product.
   %  modis_dir : directory with GEUS Greenland_Reflectivity_<YYYY>_5km_C6.nc
   %      files; when given, adds a daily MODIS albedo channel at the site.
   %      Reference layout: /Volumes/S03/DATA/greenland/geus/albedo/gris.
   %  dt : "1hr" (default; linear interpolation to hourly, the legacy
   %      behavior) or "3hr" (native posting)
   %  method : point sampling "nearest" (default) | "natural". When the
   %      FGRN11 topography is available, both methods exclude cells outside
   %      its rounded IceMask (fraction >= 0.5; Promicemask > 0 fallback) and
   %      reject a nearest valid cell farther than one native grid diagonal.
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
      source_dir = string(getenv("ICEMODEL_RACMO_DIR"));
      if source_dir == ""
         source_dir = "/Volumes/S03/DATA/greenland/racmo2p3/subsurface";
      end
   end
   if ~isfolder(source_dir)
      error('icemodel:forcing:buildRacmoData:sourceNotFound', ...
         ['RACMO source directory not found: %s. Pass source_dir or ' ...
         'stage the per-variable files (reference layout: ' ...
         '/Volumes/S03/DATA/greenland/racmo2p3/subsurface).'], source_dir)
   end

   % Channel table: file prefix -> output name.
   channels = {
      'swsd', 'swd'     ; 'lwsd', 'lwd'      ; 'swsn', 'swn'
      'lwsn', 'lwn'     ; 'senf', 'shf'      ; 'latf', 'lhf'
      'precip', 'ppt'   ; 'snowmelt', 'melt' ; 'runoff', 'runoff'
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
   [~, validmask] = racmoGridStatics(source_dir);
   if isnumeric(location) && isempty(validmask)
      error('icemodel:forcing:buildRacmoData:maskNotFound', ...
         ['RACMO point sampling requires the companion FGRN11 topography ' ...
         'IceMask; none was found in %s or its parent directory'], source_dir)
   end
   dx = hypot(diff(X, 1, 1), diff(Y, 1, 1));
   dy = hypot(diff(X, 1, 2), diff(Y, 1, 2));
   edge_lengths = [dx(:); dy(:)];
   edge_lengths = edge_lengths(isfinite(edge_lengths) & edge_lengths > 0);
   max_point_distance = sqrt(2) * median(edge_lengths);

   % Accept one location (1x2 point or polyshape, returns a single Data
   % timetable) OR a list of N points (Nx2 [lat lon], returns a 1xN cell of
   % Data timetables). N=1 is the single-point path. The per-variable files
   % and the grid above are read ONCE for the whole list; each point's grid
   % hyperslab + collapse rule is resolved here, then every variable file is
   % opened ONCE and each point's hyperslab sliced from that single open
   % (decisive for the multi-GB subsurface files).
   [locations, batch] = ...
      icemodel.forcing.helpers.normalizeLocations(location);
   npts = numel(locations);

   grid = struct('X', X, 'Y', Y, 'LAT', LAT, 'LON', LON, 'proj', proj, ...
      'first', first, 'source_dir', source_dir, ...
      'method', kwargs.method, 'remap', kwargs.remap, ...
      'validmask', validmask, 'max_point_distance', max_point_distance);
   [slabs, collapses, inslabs, sites] = deal( ...
      cell(1, npts), cell(1, npts), cell(1, npts), cell(1, npts));
   loctype = "";
   for p = 1:npts
      [slabs{p}, collapses{p}, inslabs{p}, sites{p}, loctype] = ...
         resolvePoint(grid, locations{p});
   end

   % Time axis (shared by all variables): days since 1950-01-01, native
   % 3-hourly posting.
   t = ncread(first, 'time');
   t_units = ncreadatt(first, 'time', 'units');
   assert(startsWith(t_units, 'days since 1950-01-01'), ...
      'unexpected RACMO time units: %s', t_units)
   Time = datetime(1950, 1, 1, 'TimeZone', 'UTC') + days(double(t));
   keep = ismember(year(Time), years);
   if ~any(keep)
      error('icemodel:forcing:buildRacmoData:yearNotInArchive', ...
         'requested years %s not in the RACMO archive span %d-%d', ...
         mat2str(years), year(Time(1)), year(Time(end)))
   end
   Time = Time(keep);

   % Read each available channel ONCE per file (every point's hyperslab from
   % a single open), then collapse + subset per point.
   series = cell(size(channels, 1), npts);
   for n = 1:size(channels, 1)
      if ~found(n)
         continue
      end
      blocks = icemodel.forcing.readRacmo2p3(files(n), channels{n, 1}, ...
         slabs=slabs);
      for p = 1:npts
         s = collapses{p}(blocks{p});
         series{n, p} = s(keep);
      end
   end

   data_out = cell(1, npts);
   meta_out = cell(1, npts);
   for p = 1:npts
      Data = timetable(Time);
      for n = 1:size(channels, 1)
         if found(n)
            Data.(channels{n, 2}) = series{n, p};
         end
      end
      [data_out{p}, meta_out{p}] = finalizeRacmoData(Data, sites{p}, ...
         slabs{p}, inslabs{p}, loctype, years, locations{p}, files, found, ...
         source_dir, grid, kwargs);
   end

   metadata = [meta_out{:}];
   if batch
      Data = data_out;
   else
      Data = data_out{1};
   end
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

function [slab, collapse, inslab, site, loctype] = resolvePoint(grid, location)
   %RESOLVEPOINT Map one point or polygon onto a RACMO grid hyperslab.
   %
   % Returns the hyperslab as a [start; count] 2x2 (for the batch reader),
   % the collapse handle, the slab-relative metadata index, and the site
   % lat/lon (slab mean). Conservative polygon remap runs in RACMO's NATIVE
   % rotated-pole frame (the FGRN11 rlon/rlat grid is regular there;
   % reprojecting to EPSG:3413 is curvilinear). exactremap's rotated-pole
   % support handles the rotation from the CF grid mapping and weights cells
   % by the shipped true cell areas (gridarea); off-ice cells (IceMask) are
   % inpainted. Point methods and equal-weight polygons stay in the projected
   % grid; point methods exclude the same masked cells when the mask is present.
   if isa(location, 'polyshape') && grid.remap == "conservative"
      [start, count, collapse, inslab, loctype] = ...
         resolveRacmoConservative(grid.first, grid.source_dir, ...
         location, grid.proj);
   else
      if isnumeric(location)
         assert(isequal(size(location), [1 2]), ...
            'point location must be [lat lon]')
         [xq, yq] = projfwd(grid.proj, location(1), location(2));
         location = [xq, yq];
      end
      [start, count, collapse, inslab, loctype] = ...
         icemodel.forcing.helpers.gridLocation(grid.X, grid.Y, location, ...
         grid.method, remap=grid.remap, validmask=grid.validmask, ...
         maxdistance=grid.max_point_distance);
   end
   slab = [start; count];
   site.lat = icemodel.forcing.helpers.slabMean(grid.LAT, start, count, inslab);
   site.lon = icemodel.forcing.helpers.slabMean(grid.LON, start, count, inslab);
end

function [Data, metadata] = finalizeRacmoData(Data, site, slab, inslab, ...
      loctype, years, location, files, found, source_dir, grid, kwargs)
   %FINALIZERACMODATA Post-process one point's assembled RACMO Data + metadata.
   % Identical to the legacy single-point tail: hourly interpolation, derived
   % albedo, optional MODIS channel, precip rate, units, metchecks, userdata
   % CustomProperties, and the provenance struct. start/count come from slab.
   start = slab(1, :);
   count = slab(2, :);

   % Interpolate to hourly (legacy behavior) unless native requested.
   % The full-year hourly axis extends past the last 3-hourly posting
   % (21:00 on Dec 31); the trailing hours extrapolate linearly.
   if kwargs.dt == "1hr"
      Time = Data.Time;
      t1 = dateshift(Time(1), 'start', 'year');
      t2 = dateshift(Time(end), 'start', 'year') + calyears(1) - hours(1);
      t_hourly = (t1:hours(1):t2)';
      t_hourly = t_hourly(ismember(year(t_hourly), years));
      Data = retime(Data, t_hourly, 'linear', 'EndValues', 'extrap');
   end

   % Site location (also needed for the optional MODIS channel below).
   [site_x, site_y] = projfwd(grid.proj, site.lat, site.lon);

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

   % Optional GEUS MODIS daily albedo, resolved at the SAME location and
   % aggregation as the RACMO channels: nearest/natural for a point,
   % conservative (or equal) catchment mean for a polygon. This generalizes
   % the legacy saveRacmoData MODIS channel (which took the nearest cell even
   % for catchments) to the area-weighted ROI mean for polygons.
   modis_metadata = struct();
   if kwargs.modis_dir ~= ""
      [modis, modis_metadata] = ...
         icemodel.forcing.helpers.modisAlbedoChannel( ...
         kwargs.modis_dir, years, location, kwargs.method, kwargs.remap, ...
         Data.Time);
      if ~isempty(modis_metadata.modis_coverage_years)
         Data.modis = modis;
      end
   end

   % Precipitation to the canonical water-equivalent rate m s-1. RACMO posts
   % precip as mWE/h, so dividing by 3600 s/h yields m s-1. The diagnostic
   % mass fluxes (melt/runoff/smb/refreeze/subl/sndiv/meltin) keep mWE/h.
   Data.ppt = Data.ppt / 3600;

   % Normalize RACMO's native mass-balance sign at the source boundary.
   % Native negative loss / positive deposition is the exact opposite of the
   % public positive-loss / negative-deposition sublimation contract.
   if ismember("subl", string(Data.Properties.VariableNames))
      Data.subl = -Data.subl;
   end

   % RACMO's physical precipitation flux contains small negative numerical
   % undershoots at the source grid cells. Enforce the nonnegative invariant at
   % the source-finalization boundary, after spatial sampling/remapping and
   % hourly interpolation, and retain the exact input minimum/count in metadata.
   [Data, ppt_qc] = ...
      icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl(Data);

   % Per-variable units from the shared canonical map.
   Data.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(Data.Properties.VariableNames));

   % Source-faithful QA/QC applies physical clamps without inventing values
   % across native outages. The legacy caller can still request gap filling at
   % its own boundary if needed.
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, fillgaps=false);

   % Attach the shared location schema. Elevation comes from the companion
   % topography channel when present; RACMO supplies no surface slope here.
   location_metadata = struct( ...
      'lat_wgs84', site.lat, 'lon_wgs84', site.lon, ...
      'x_epsg3413', site_x, 'y_epsg3413', site_y, ...
      'elev_m', readElevation(source_dir, start, count, inslab));
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, location_metadata);

   metadata = struct( ...
      'source_files', files(found), ...
      'location_type', loctype, ...
      'method', kwargs.method, ...
      'remap', kwargs.remap, ...
      'grid_start', start, ...
      'grid_count', count, ...
      'n_cells', prod(count), ...
      'lat', site.lat, 'lon', site.lon, ...
      'racmo_ice_mask_applied', ~isempty(grid.validmask), ...
      'racmo_point_max_distance_m', grid.max_point_distance, ...
      'dt', kwargs.dt, ...
      'racmo_subl_native_sign_convention', ...
      "negative_loss_positive_deposition", ...
      'racmo_subl_sign_convention', ...
      "positive_loss_negative_deposition", ...
      'mass_flux_units', ...
      "ppt m s-1; diagnostic fluxes mWE/h (rate; cumulative sums need dt hours)", ...
      'checks', checks);

   % Copy exact GEUS product/status/year provenance resolved by the channel
   % reader without any second inventory or NetCDF access.
   modis_fields = fieldnames(modis_metadata);
   for k = 1:numel(modis_fields)
      metadata.(modis_fields{k}) = modis_metadata.(modis_fields{k});
   end

   % Copy the flat precipitation-QC contract so writeuserdata carries it with
   % artifact metadata and the bounded repair can apply the identical policy.
   fields = fieldnames(ppt_qc);
   for k = 1:numel(fields)
      metadata.(fields{k}) = ppt_qc.(fields{k});
   end
   metadata = icemodel.forcing.helpers.columnizeMetadata(metadata);
   Data.Properties.UserData = metadata;
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
      ice_mask = double(ncread(filename, 'IceMask'));
      validmask = isfinite(ice_mask) & ice_mask >= 0.5;
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
   % apply the weights ourselves: the polygon geometry is solved ONCE and reused
   % across every year's data block (the legacy runoff static-weights pattern).
   % exactremap-0hv (the earlier rotated-pole 'areaavg' ~0) is fixed, so areaavg
   % is also correct now; weights mode is retained for the compute-once reuse.
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
   % NDGRID [numel(rlon) x numel(rlat)] (rlon down rows). For gridvector axes
   % this is a pure transpose (no flips), so transpose the cell-area/mask slabs
   % going in and the returned W back to ndgrid coming out; exactremap now
   % validates this orientation (rasterTransposed / ambiguousGridOrientation),
   % so a mismatch errors rather than misaligns. Weights mode is kept (over the
   % now-fixed exactremap-0hv 'areaavg' path) because it solves the polygon
   % geometry ONCE and reuses it across every year's data block via the collapse
   % closure below. (A coordinate-list 'weights' call computes the IDENTICAL
   % weights - verified: same sum, bit-identical sorted set - but returns them in
   % exactremap's reconstructed-grid order, not the input-list order, so it needs
   % an index remap-back and is no simpler than this transpose. The clean
   % migration is 'areaavg' on the multi-page V stack, which returns the
   % aggregate directly with no per-cell weight to reorder - see icemodel-1ps.18.)
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
