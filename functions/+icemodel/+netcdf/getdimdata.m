function dims = getdimdata(Z, dz)
   %GETDIMDATA Get dimension data
   %
   %  DIMS = GETDIMDATA(Z, DZ)
   %
   % See also:

   arguments
      Z  (1, 1) double {mustBeNumeric} = 0
      dz (1, 1) double {mustBeNumeric} = 0
   end

   % Load the grid coordinates and elevation.
   [dims.x_easting, dims.y_northing, dims.elevation, ...
      dims.latitude, dims.longitude] = loadIceMask("icemask", ...
      varnames=["X", "Y", "Elev", "Lat", "Lon"]);

   % Define the grid cell index
   dims.gridcell = 1:numel(dims.elevation);

   % Define the layer depth dimension
   dims.depth = dz/2:dz:Z;
end
