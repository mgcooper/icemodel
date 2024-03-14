function writeice2(ncid, vars, dims, filepath)
   %WRITEICE2 Write ice2 data to icemodel nc file
   %
   % WRITEICE2(NCID, DATA, VARS)
   %
   %  access pattern: [n-1 0 0], [1 nlyrs nhrs]
   %
   %  Make sure the [start], [count] matches the dimid and the data
   %  orientation. Here [n-1 0 0], [1 nz nt]
   %  where nz is number of depth layers and nt number of time slices.
   %  The netcdf4 data model might allow the shape of the data sent in
   %  to differ, e.g., the putVar step works when ice2.(thisvar) is
   %  transposed, but the data is scrambled.
   %
   % See also:

   % Write the 2d variables
   for v = 1:numel(vars)
      thisvar = vars{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);

      % Writes one cell at a time
      for n = 1:numel(dims.gridcell)
         load(fullfile(filepath, ['ice2_' num2str(n) '.mat']), 'ice2');
         netcdf.putVar(ncid, thisvid, [n-1 0 0], [1 size(ice2.(thisvar))], ice2.(thisvar));
      end

      % This writes all data at once, which does not work with 16 GB ram
      % netcdf.putVar(ncid, thisvid, permute(data.(thisvar), [3 1 2]));

      % To preallocate in chunks:
      % chunksize = 10;
      % for n = 1:chunksize:numel(dims.gridcell)
      %    [s, e] = chunkLoopInds(n, 1, chunksize);
      %    tmp = preallocateDataArrays(vars, chunksize, 8760, 500, 'NC_FLOAT', 'ice2');
      %    for m = s:e
      %       tmp.ice2.(thisvar)(:, :, m) = ...
      %          load(fullfile(filepath, ['ice2_' num2str(m) '.mat'])).('ice2').(thisvar);
      %    end
      %    % The permute rearranges [depth time gridcell] to [gridcell depth time]
      %    netcdf.putVar(ncid, thisvid, [s-1 0 0], size(permute(tmp.ice2.(thisvar), [3 1 2])), ...
      %       permute(tmp.ice2.(thisvar), [3 1 2]));
      % end

   end

   % NOTE: The methods which use start, count should work if
   % size(ice2.(thisvar)) is used for count, but also recall this method
   % from the allocateData subfunction:
   %
   % data.(thisvar)(1:size(ice2.(thisvar), 1), :, n) = ice2.(thisvar); % nlyrs x nhrs x ncells

   % copied this out of allocate data b/c for ice2 we should loop over cells and
   % write all vars
   %
   % for n = 1:ncells
   %    ice2 = load(fullfile(filepath, ['ice2_' num2str(n) '.mat'])).('ice2');
   %    for v = 1:numel(vars.ice2)
   %       thisvar = vars.ice2{v};
   %       data.(thisvar)(1:size(ice2.(thisvar), 1), :, n) = ice2.(thisvar); % nlyrs x nhrs x ncells
   %    end
   % end
end
