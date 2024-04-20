function writeice1(ncid, datapath, varnames, dimdata, xtype, simmodel)
   %WRITEICE1 Write ice1 data to icemodel nc file
   %
   % WRITEICE1(NCID, VARS, DATA)
   %
   % The access pattern is: [gridcell, timestep], thus:
   % start = [n-1 0], count = [1 nhrs].
   %
   % See also: icemodel.netcdf.writeice2

   arguments
      ncid (1, 1) double {mustBeNumeric}
      datapath (1, :) char {mustBeFolder}
      varnames (1, :) cell {mustBeText}
      dimdata (1, :) struct {mustBeStruct}
      xtype (1, :) char {mustBeMember(xtype, ...
         {'NC_FLOAT', 'NC_DOUBLE', 'NC_INT64', 'NC_UINT64', 'NC_INT', ...
         'NC_UINT', 'NC_SHORT', 'NC_USHORT', 'NC_BYTE', 'NC_UBYTE', ...
         'NC_CHAR', 'NC_STRING'})}
      simmodel (1, :) char {mustBeTextScalar}
   end

   % Fill the data arrays in memory.
   data = icemodel.netcdf.getvardata( ...
      datapath, varnames, dimdata, xtype, 1, nan, simmodel);

   % Write the data arrays to netcdf.
   for v = 1:numel(varnames)
      thisvar = varnames{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);
      netcdf.putVar(ncid, thisvid, data.(thisvar).');
   end

   % To write one cell at a time:
   % for n = 1:numel(dimdata.gridcell)
   %    for v = 1:numel(varnames)
   %       thisvar = varnames{v};
   %       thisvid = netcdf.inqVarID(ncid, thisvar);
   %       netcdf.putVar(ncid, thisvid, [n-1 0], [1 nhrs], load(fullfile(...
   %          datapath, ['ice1_' num2str(n) '.mat']).('ice1').(thisvar)));
   %    end
   % end
end
