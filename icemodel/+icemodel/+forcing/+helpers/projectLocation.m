function location = projectLocation(location)
   %PROJECTLOCATION Ensure a location struct carries EPSG:3413 coordinates.
   %
   %  location = icemodel.forcing.helpers.projectLocation(location)
   %
   % Source builders and manifest importers both need the same WGS84 +
   % polar-stereographic location metadata. This helper fills missing or NaN
   % x_epsg3413 / y_epsg3413 values from lat_wgs84 / lon_wgs84 and leaves
   % already-populated projected coordinates untouched.

   has_xy = isfield(location, 'x_epsg3413') ...
      && isfield(location, 'y_epsg3413') ...
      && all(isfinite([location.x_epsg3413, location.y_epsg3413]));
   has_ll = isfield(location, 'lat_wgs84') ...
      && isfield(location, 'lon_wgs84') ...
      && all(isfinite([location.lat_wgs84, location.lon_wgs84]));

   if ~has_xy && has_ll
      proj = icemodel.forcing.helpers.psnProjection();
      [x, y] = projfwd(proj, location.lat_wgs84, location.lon_wgs84);
      location.x_epsg3413 = x;
      location.y_epsg3413 = y;
   elseif ~isfield(location, 'x_epsg3413')
      location.x_epsg3413 = NaN;
      location.y_epsg3413 = NaN;
   elseif ~isfield(location, 'y_epsg3413')
      location.y_epsg3413 = NaN;
   end
end
