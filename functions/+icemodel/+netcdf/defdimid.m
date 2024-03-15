function dimid = defdimid(ncid, ncells, nhrs, nlyrs)
   %DEFDIMID Define icemodel netcdf file dimensions
   %
   %  DIMID = DEFDIMID(NCID, NCELLS, NHRS, NLYRS)
   %
   % CF conventions suggest but do not require the dimensions appear in the
   % order: T, Z, Y, X. However, the first dimension in matlab is the last in
   % netcdf. Thus the data arrays are defined here as:
   %
   % ice1: gridcell x time
   % ice2: gridcell x depth x time
   %
   % And they appear in the netcdf files as:
   %
   % ice1: time x gridcell
   % ice2: time x depth x gridcell.
   %
   % Note that a "coordinate variable" is a "1-d variable with the same name as
   % its dimension, e.g., time(time), and it is defined as a numeric data type
   % with values that are ordered monotonically".
   %
   % See also:

   % Change dimid.ice1/2 to dimid.data for standard syntax in makencfile

   % Define the dimensions of the data arrays, in strict order.
   if nargin < 4 || nlyrs == 1

      dimid.gridcell = netcdf.defDim(ncid, 'gridcell', ncells);
      dimid.time = netcdf.defDim(ncid, 'time', nhrs);
      dimid.ice1 = [dimid.gridcell dimid.time];
   else

      dimid.gridcell = netcdf.defDim(ncid, 'gridcell', ncells);
      dimid.depth = netcdf.defDim(ncid, 'depth', nlyrs);
      dimid.time = netcdf.defDim(ncid, 'time', nhrs);
      dimid.ice2 = [dimid.gridcell dimid.depth dimid.time];
   end

   % For a lat lon or x y:
   %
   % dimid.ice1 = [dimid.lon dimid.lat dimid.time];
   % dimid.ice2 = [dimid.lon dimid.lat dimid.depth dimid.time];
   %
   % Notes from CF 1.11:
   %
   % If any or all of the dimensions of a variable have the interpretations of
   % "date or time" (T), "height or depth" (Z), "latitude" (Y), or "longitude"
   % (X) then we recommend, but do not require (see Section 1.5, "Relationship
   % to the COARDS Conventions"), those dimensions to appear in the relative
   % order T, then Z, then Y, then X in the CDL definition corresponding to the
   % file. All other dimensions should, whenever possible, be placed to the left
   % of the spatiotemporal dimensions.
   %
   % Also note:
   %
   % "we allow but do not require the units attribute of dimensionless vertical
   % coordinates to take the values "level", "layer", or "sigma_level.""
   %
   % But I was unable to find a similar "units" (or "standard_name") for grid
   % cell index
end
