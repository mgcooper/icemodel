function chunksize = getchunksize(whichdata, dimdata, datasize, opts)

   arguments
      whichdata
      dimdata
      datasize
      opts.GetSizeFromData (1, :) logical {mustBeNumericOrLogical} = true
      opts.GetSizeFromDims (1, :) logical {mustBeNumericOrLogical} = false
   end

   % The number of gridcells is set by the dimensions in both cases
   dimsizes = icemodel.netcdf.getdimsize(dimdata);
   numcells = dimsizes.gridcell;

   % Use GetSizeFromData to update the depth dimension values directly from the
   % size of the ice2 data set in getvarinfo. Note that this will not catch the
   % case where the size changes from file to file within a year.

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
   % When writing, the primary concern is reducing the number of writes.
   % Thus writing the entire array at once is typically ideal, and the
   % netcdf software will determine the chunk size.
   %
   % However, if typical access patterns are known, chunking can improve
   % efficiency by ensuring the data layout in the file matches how it is
   % accessed later. If data is predominantly accessed in large contiguous
   % blocks, having the data and chunks aligned to those patterns is ideal.
end
