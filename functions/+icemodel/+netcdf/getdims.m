function dims = getdims(Z, dz)
   %GETDIMS Get dimensions
   %
   %  DIMS = GETDIMS(Z, DZ)
   %
   % See also:

   arguments
      Z = []
      dz = []
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
