function dimid = defdimid(ncid, dimdata, datasize, opts)
   %DEFDIMID Define icemodel netcdf file dimensions
   %
   %  DIMID = DEFDIMID(NCID, DIMDATA)
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

   % NOTE: With the current setup, GetSizeFromData should be identical to
   % GetSizeFromDims, because getdimsize sets the max(1, ...) check on depth.
   % The case where they would differ is if the datasize changes within years.

   arguments
      ncid  (1, 1) double {mustBeNumeric}
      dimdata (1, :) struct {mustBeStruct}
      datasize (1, :) double {mustBeNumeric}
      opts.GetSizeFromData (1, :) logical {mustBeNumericOrLogical} = true
      opts.GetSizeFromDims (1, :) logical {mustBeNumericOrLogical} = false
   end

   % Change dimid.ice1/2 to dimid.data for standard syntax in makencfile

   % Need to either use dimdata or get the sizes directly from the data.
   dimsize = icemodel.netcdf.getdimsize(dimdata);
   numcells = dimsize.gridcell;

   if opts.GetSizeFromDims
      numlayers = dimsize.depth;
      numtimesteps = dimsize.time;

   elseif opts.GetSizeFromData
      numlayers = datasize(1);
      numtimesteps = datasize(2);
   end

   % Define the dimensions of the data arrays, in strict order.
   if numlayers == 1

      dimid.gridcell = netcdf.defDim(ncid, 'gridcell', numcells);
      dimid.time = netcdf.defDim(ncid, 'time', numtimesteps);
      dimid.ice1 = [dimid.gridcell dimid.time];
   else

      dimid.gridcell = netcdf.defDim(ncid, 'gridcell', numcells);
      dimid.depth = netcdf.defDim(ncid, 'depth', numlayers);
      dimid.time = netcdf.defDim(ncid, 'time', numtimesteps);
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
