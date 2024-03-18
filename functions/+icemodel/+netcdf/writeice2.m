function writeice2(ncid, datapath, varnames, dimdata, xtype, simmodel, maxGB)
   %WRITEICE2 Write ice2 data to icemodel nc file
   %
   % WRITEICE2(NCID, DATA, VARS)
   %
   %  The access pattern is: [gridcell, layer, timestep], thus:
   %  start = [n-1 0 0], count = [1 nlyrs nhrs]
   %
   %  Make sure the [start], [count] matches the dimid and the data
   %  orientation. Here [n-1 0 0], [1 nz nt]
   %  where nz is number of depth layers and nt number of time slices.
   %  The netcdf4 data model might allow the shape of the data sent in
   %  to differ, e.g., the putVar step works when ice2.(thisvar) is
   %  transposed, but the data is scrambled.
   %
   % See also:

   arguments
      ncid (1, 1) double {mustBeNumeric}
      datapath (1, :) char {mustBeFolder}
      varnames (1, :) cell {mustBeText}
      dimdata (1, :) struct {mustBeStruct}
      xtype (1, :) char {mustBeMember(xtype, ...
         {'NC_FLOAT', 'NC_DOUBLE', 'NC_INT64', 'NC_UINT64', 'NC_INT', ...
         'NC_UINT', 'NC_SHORT', 'NC_USHORT', 'NC_BYTE', 'NC_UBYTE', ...
         'NC_CHAR', 'NC_STRING'})}
      simmodel (1, :) char {mustBeTextScalar} = ""
      maxGB (1, 1) double {mustBeNumeric} = 8
   end


   dimsizes = icemodel.netcdf.getdimsize(dimdata);

   % Write in chunks of numcells, where numcells maximizes in-memory array size.
   numcells = icemodel.netcdf.maxcells(maxGB, dimsizes.depth, dimsizes.time, ...
      numel(varnames), nctype2mat(xtype));

   for startcell = 1:numcells:numel(dimdata.gridcell)

      countcell = min(numcells, numel(dimdata.gridcell) - startcell + 1);

      fprintf('working on cells %d:%d ... \n', startcell, startcell + countcell - 1)
      % [startcell startcell+countcell-1] % use to confirm start/count

      data = icemodel.netcdf.getvardata(datapath, varnames, dimdata, xtype, ...
         startcell, countcell); % layers x time x cells

      % Define the netcdf api start, count
      start = [startcell-1 0 0]; % thiscell x layer(1) x timestep(1)

      for v = 1:numel(varnames)
         thisvar = varnames{v};
         thisvid = netcdf.inqVarID(ncid, thisvar);

         thisdata = permute(data.(thisvar), [3 1 2]); % cells x layers x time
         count = size(thisdata);

         % The permute rearranges [depth time gridcell] to [gridcell depth time]
         netcdf.putVar(ncid, thisvid, start, count, thisdata);
      end
   end

   % % Write one cell at a time (use numcell = 1 in the chunk method instead)
   % for startcell = 1:numel(dims.gridcell)
   %
   %    data = load(fullfile(datapath, ['ice2_' num2str(startcell) '.mat'])).('ice2');
   %
   %    start = [startcell-1 0 0]; % thiscell x layer(1) x timestep(1)
   %
   %    % Write the 2d variables
   %    for v = 1:numel(varnames)
   %       thisvar = varnames{v};
   %       thisvid = netcdf.inqVarID(ncid, thisvar);
   %
   %       thisdata = data.(thisvar);
   %       count = [1 size(thisdata)]; % numcells x numlayers x numtimesteps
   %
   %       netcdf.putVar(ncid, thisvid, start, count, thisdata);
   %    end
   % end

   % If all data fit in memory, this writes it in one go:
   % for v = 1:numel(varnames)
   %    thisvar = varnames{v};
   %    thisvid = netcdf.inqVarID(ncid, thisvar);
   %    netcdf.putVar(ncid, thisvid, permute(data.(thisvar), [3 1 2]));
   % end
end
