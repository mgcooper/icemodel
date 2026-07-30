function [albedo, Time, selection] = readGeusModis(filename, location, method, kwargs)
   %READGEUSMODIS Read the GEUS MODIS daily albedo at points or a polygon.
   %
   %  [albedo, Time] = icemodel.forcing.readGeusModis(filename, [lat lon])
   %  [albedo, Time] = icemodel.forcing.readGeusModis(filename, points)
   %  [albedo, Time] = icemodel.forcing.readGeusModis(filename, polygon)
   %  [albedo, Time, selection] = ... readGeusModis(_, method, remap=...)
   %
   % Reads the daily MODIS-derived surface albedo (fraction) from one
   % GEUS Greenland reflectivity file (Greenland_Reflectivity_<YYYY>_
   % 5km_C6.nc) at one or more points or averaged over a polygon
   % (catchment), using the same cell-selection / remap machinery as the
   % gridded-source builders (icemodel.forcing.helpers.gridLocation):
   %
   %  - location = [lat lon] rows (points, degrees): each row yields the
   %    grid cell nearest that point (method "nearest", default) or a
   %    natural-neighbour blend of the surrounding cells at the point
   %    (method "natural"). Passing every station at once builds the grid
   %    geometry once and reads ONE covering albedo hyperslab per file,
   %    which is how multi-station staging (stageModisAlbedo) avoids
   %    per-station re-reads of the same yearly NetCDF.
   %  - location = polyshape (vertices in EPSG:3413 metres): the MODIS
   %    cells are averaged over the polygon. remap="conservative" (default)
   %    is the exact overlap-area-weighted catchment mean via the exactremap
   %    toolbox (helpers.remapPolygon), matching the legacy readGeusModis
   %    ALBavgInPoly ROI mean; remap="equal" is a plain in-polygon
   %    cell-centre mean. The remap runs in the GEUS 5 km polar-stereographic
   %    frame (the grid is regular there).
   %
   % Inputs
   %  filename - GEUS reflectivity NetCDF for one year
   %  location - [lat lon] point rows (degrees) or polyshape (EPSG:3413
   %             metres)
   %  method   - point sampling "nearest" (default) | "natural"
   %
   % Name-value
   %  remap - polygon aggregation "conservative" (default) | "equal"
   %
   % Outputs
   %  albedo    - daily albedo series, one column per requested point (one
   %              column for a polygon) [-]; undocumented finite 999
   %              sentinels and other nonphysical samples are returned as
   %              NaN
   %  Time      - UTC daily datetime axis (Jan 1 of the file year onward)
   %  selection - 1 x ncolumns struct of cell-selection provenance: the
   %              sampling method, hyperslab start/count indices into the
   %              (x, y) grid, and, for single-cell point extractions, the
   %              selected cell centre (cell_lat/cell_lon, degrees) plus
   %              distance_m from the query point to that cell centre in
   %              the native GEUS frame (NaN for multi-cell selections).
   %              grid_size and grid_sha256 identify the complete raw
   %              latitude/longitude grid so multi-year staging cannot mix
   %              values from geometrically different source files.
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.stageModisAlbedo,
   %  icemodel.forcing.helpers.gridLocation,
   %  icemodel.forcing.helpers.dailyToHourly

   arguments
      filename (1, 1) string
      location
      method (1, 1) string {mustBeMember(method, ...
         ["nearest", "natural"])} = "nearest"
      kwargs.remap (1, 1) string {mustBeMember(kwargs.remap, ...
         ["equal", "conservative"])} = "conservative"
   end

   LON = double(ncread(filename, 'lon'));
   LAT = double(ncread(filename, 'lat'));
   grid_sha256 = coordinateGridSha256(LAT, LON);

   % Project the GEUS 5 km grid into its native polar-stereographic frame
   % (sphere, true scale 71N, central meridian 39W; see geusModisProjection),
   % where the 5 km posting is axis-aligned and uniform, then rebuild exactly
   % regular axes via linspace (a sub-metre correction over the 5 km cell) so
   % the conservative remap operates on a perfectly regular grid - the GEUS
   % analogue of the MAR Xnat/Ynat axes.
   geus = icemodel.forcing.helpers.geusModisProjection();
   [Xp, Yp] = projfwd(geus, LAT, LON);
   xax = linspace(mean(Xp(1, :)), mean(Xp(end, :)), size(Xp, 1)).';
   yax = linspace(mean(Yp(:, 1)), mean(Yp(:, end)), size(Yp, 2));
   [X, Y] = ndgrid(xax, yax);

   if isnumeric(location)
      % Point rows: one [lat lon] per requested series. A row list shares one
      % geometry build and one covering hyperslab read across all stations,
      % so multi-station staging never re-reads the same yearly file.
      assert(~isempty(location) && size(location, 2) == 2, ...
         'point locations must be [lat lon] rows')
      [xq, yq] = projfwd(geus, location(:, 1), location(:, 2));
      queries = num2cell([xq, yq], 2);
   else
      % Polygon vertices arrive in EPSG:3413 metres; map them into the GEUS
      % native frame (3413 metres -> lon/lat -> native metres) so the polygon
      % and the grid share the regular frame the remap operates in.
      proj = icemodel.forcing.helpers.psnProjection();
      [vlat, vlon] = projinv(proj, location.Vertices(:, 1), ...
         location.Vertices(:, 2));
      [qx, qy] = projfwd(geus, vlat, vlon);
      queries = {polyshape(qx, qy)};
   end

   % Map each target onto a grid hyperslab + collapse rule, reusing the
   % shared selection logic so the MODIS channel honours the same
   % point/polygon and nearest/natural/conservative/equal options as the
   % other gridded channels.
   nq = numel(queries);
   starts = zeros(nq, 2);
   counts = zeros(nq, 2);
   collapses = cell(nq, 1);
   for q = 1:nq
      [starts(q, :), counts(q, :), collapses{q}] = ...
         icemodel.forcing.helpers.gridLocation( ...
         X, Y, queries{q}, method, remap=kwargs.remap);
   end

   info = ncinfo(filename, 'albedo');
   ndays = info.Size(end);

   % Read ONE bounding hyperslab that covers every query over the two grid
   % dimensions and all days, then slice each query's block from memory and
   % flatten the cells column-major (cells x time) so the gridLocation
   % collapse handle reduces it to the target series exactly as it does for
   % the MAR / RACMO / MERRA channels. For a single query this is the same
   % read as the original per-target hyperslab.
   union_start = min(starts, [], 1);
   union_count = max(starts + counts, [], 1) - union_start;
   union_block = ncread(filename, 'albedo', ...
      [union_start 1], [union_count ndays]);

   % One output column per query, plus the selection provenance staging needs
   % to pin each extraction (indices, selected cell centre, query offset)
   % without re-deriving the grid geometry outside this reader.
   albedo = nan(ndays, nq);
   selection = repmat(struct('method', char(method), 'start', [0 0], ...
      'count', [0 0], 'cell_lat', NaN, 'cell_lon', NaN, ...
      'distance_m', NaN, 'grid_size', size(LAT), ...
      'grid_sha256', char(grid_sha256)), 1, nq);
   for q = 1:nq
      rows = starts(q, 1) - union_start(1) + (1:counts(q, 1));
      cols = starts(q, 2) - union_start(2) + (1:counts(q, 2));
      block = union_block(rows, cols, :);
      % Raw C6 files do not declare their finite 999 missing-data sentinel.
      % Mask nonphysical albedo before spatial collapse so no sentinel can
      % contaminate point interpolation or polygon aggregation.
      block = icemodel.forcing.helpers.normalizeGeusModisAlbedo(block);
      block = reshape(block, prod(counts(q, :)), ndays);
      % Collapse values and their validity weights separately. This preserves
      % a valid polygon/interpolated mean when only part of the selected slab
      % is missing, while an entirely missing target/day remains NaN.
      valid = isfinite(block);
      block(~valid) = 0;
      numerator = collapses{q}(block);
      denominator = collapses{q}(double(valid));
      series = numerator ./ denominator;
      series(denominator <= 0) = NaN;
      albedo(:, q) = ...
         icemodel.forcing.helpers.normalizeGeusModisAlbedo(series(:));

      selection(q).start = starts(q, :);
      selection(q).count = counts(q, :);
      % Cell-centre identity and query offset are only well defined when the
      % selection is one point mapped to one cell; blended or polygon
      % selections keep NaN placeholders.
      if isnumeric(queries{q}) && isequal(counts(q, :), [1 1])
         cell_x = X(starts(q, 1), starts(q, 2));
         cell_y = Y(starts(q, 1), starts(q, 2));
         [selection(q).cell_lat, selection(q).cell_lon] = ...
            projinv(geus, cell_x, cell_y);
         selection(q).distance_m = hypot( ...
            cell_x - queries{q}(1), cell_y - queries{q}(2));
      end
   end

   % Parse only the source basename: parent work directories can themselves
   % contain underscore-delimited four-digit tokens unrelated to the product.
   [~, source_name, source_ext] = fileparts(filename);
   if strcmpi(source_ext, '.gz')
      [~, source_name] = fileparts(source_name);
   end
   tok = regexp(source_name, '_(\d{4})_', 'tokens', 'once');
   assert(~isempty(tok), ...
      'cannot parse the file year from %s', filename)
   t0 = datetime(str2double(tok{1}), 1, 1, 'TimeZone', 'UTC');
   Time = (t0:days(1):(t0 + days(ndays - 1)))';
end

function digest = coordinateGridSha256(lat, lon)
   %COORDINATEGRIDSHA256 Identify one complete raw latitude/longitude grid.
   %
   % Include dimensions before both double arrays so equal byte sequences
   % with different shapes cannot identify as the same spatial grid.
   bytes = [typecast(uint64(size(lat)), 'uint8')'; ...
      typecast(lat(:), 'uint8'); typecast(lon(:), 'uint8')];
   md = java.security.MessageDigest.getInstance('SHA-256');
   raw = typecast(md.digest(bytes), 'uint8');
   digest = string(lower(reshape(dec2hex(raw, 2)', 1, [])));
end
