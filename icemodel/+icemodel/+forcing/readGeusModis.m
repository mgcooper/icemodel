function [albedo, Time] = readGeusModis(filename, location, method, kwargs)
   %READGEUSMODIS Read the GEUS MODIS daily albedo at a point or polygon.
   %
   %  [albedo, Time] = icemodel.forcing.readGeusModis(filename, [lat lon])
   %  [albedo, Time] = icemodel.forcing.readGeusModis(filename, polygon)
   %  [albedo, Time] = ... readGeusModis(_, method, remap=...)
   %
   % Reads the daily MODIS-derived surface albedo (fraction) from one
   % GEUS Greenland reflectivity file (Greenland_Reflectivity_<YYYY>_
   % 5km_C6.nc) at a point or averaged over a polygon (catchment), using
   % the same cell-selection / remap machinery as the gridded-source
   % builders (icemodel.forcing.helpers.gridLocation):
   %
   %  - location = [lat lon] (point, degrees): the grid cell nearest the
   %    point (method "nearest", default) or a natural-neighbour blend of
   %    the surrounding cells at the point (method "natural").
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
   %  location - [lat lon] point (degrees) or polyshape (EPSG:3413 metres)
   %  method   - point sampling "nearest" (default) | "natural"
   %
   % Name-value
   %  remap - polygon aggregation "conservative" (default) | "equal"
   %
   % Outputs
   %  albedo - daily albedo series at the target [-]
   %  Time   - UTC daily datetime axis (Jan 1 of the file year onward)
   %
   % See also: icemodel.forcing.buildMarData,
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
      assert(isequal(size(location), [1 2]), ...
         'point location must be [lat lon]')
      [xq, yq] = projfwd(geus, location(1), location(2));
      query = [xq, yq];
   else
      % Polygon vertices arrive in EPSG:3413 metres; map them into the GEUS
      % native frame (3413 metres -> lon/lat -> native metres) so the polygon
      % and the grid share the regular frame the remap operates in.
      proj = icemodel.forcing.helpers.psnProjection();
      [vlat, vlon] = projinv(proj, location.Vertices(:, 1), ...
         location.Vertices(:, 2));
      [qx, qy] = projfwd(geus, vlat, vlon);
      query = polyshape(qx, qy);
   end

   % Map the target onto a grid hyperslab + collapse rule, reusing the
   % shared selection logic so the MODIS channel honours the same
   % point/polygon and nearest/natural/conservative/equal options as the
   % other gridded channels.
   [start, count, collapse] = icemodel.forcing.helpers.gridLocation( ...
      X, Y, query, method, remap=kwargs.remap);

   info = ncinfo(filename, 'albedo');
   ndays = info.Size(end);

   % Read the bounding hyperslab over the two grid dimensions, all days, and
   % flatten the cells column-major (cells x time) so the gridLocation
   % collapse handle reduces it to the target series exactly as it does for
   % the MAR / RACMO / MERRA channels.
   block = double(ncread(filename, 'albedo', [start 1], [count ndays]));
   block = reshape(block, prod(count), ndays);
   albedo = collapse(block);

   tok = regexp(filename, '_(\d{4})_', 'tokens', 'once');
   assert(~isempty(tok), ...
      'cannot parse the file year from %s', filename)
   t0 = datetime(str2double(tok{1}), 1, 1, 'TimeZone', 'UTC');
   Time = (t0:days(1):(t0 + days(ndays - 1)))';
end
