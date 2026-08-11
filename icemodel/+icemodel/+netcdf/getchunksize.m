function chunksize = getchunksize(whichdata, dimdata, datasize, opts)

   arguments
      whichdata
      dimdata (1, :) struct {mustBeStruct}
      datasize (1, :) double {mustBeNumeric}
      opts.GetSizeFromData (1, :) logical {mustBeNumericOrLogical} = true
      opts.GetSizeFromDims (1, :) logical {mustBeNumericOrLogical} = false
   end

   % The number of gridcells is set by the dimensions in both cases
   dimsizes = icemodel.netcdf.getdimsize(dimdata);
   numcells = dimsizes.gridcell;

   % Use GetSizeFromData to update the depth dimension values directly from the
   % size of the ice2 data set in getvarinfo. This does not catch the case
   % where the size changes from file to file within a year.

   if opts.GetSizeFromDims
      numlayers = dimsizes.depth;
      numtimesteps = dimsizes.time;

   elseif opts.GetSizeFromData
      numlayers = datasize(1);
      numtimesteps = datasize(2);
   end

   % Set the chunksizes
   switch whichdata
      case 'ice1'
         chunksize = [numcells, numtimesteps];   % all cells, annual chunks

      case 'ice2'
         chunksize = [1, numlayers, numtimesteps]; % one cell, all layers, annual
   end

   % Define chunkSize based on data access patterns. Larger chunk sizes
   % increase memory usage during read/write.
   %
   % When writing, the main concern is the number of writes. Writing the
   % entire array at once is usually best, and the netcdf software then
   % determines the chunk size.
   %
   % If the typical access patterns are known, chunking can improve
   % efficiency, because the data layout in the file then matches the later
   % access. If readers mostly access the data in large contiguous blocks,
   % align the data and the chunks to those blocks.
end
