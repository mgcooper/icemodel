function defdatavars(ncid, dimid, vars, standardnames, longnames, ...
      units, chunksize, xtype, dochunking, ncprops)
   %DEFDATAVARS Define the icemodel data variables and attributes.
   %
   %  DEFDATAVARS(NCID, DIMID, VARS, STANDARDNAMES, LONGNAMES, UNITS, AXES, NCPROPS)
   %
   % See also:

   arguments
      ncid (1, 1) {mustBeNumeric}
      dimid
      vars
      standardnames
      longnames
      units
      chunksize
      xtype (1, :) char {mustBeMember(xtype, ...
         {'NC_FLOAT', 'NC_DOUBLE', 'NC_INT64', 'NC_UINT64', 'NC_INT', ...
         'NC_UINT', 'NC_SHORT', 'NC_USHORT', 'NC_BYTE', 'NC_UBYTE', ...
         'NC_CHAR', 'NC_STRING'})} = 'NC_DOUBLE'
      dochunking (1, 1) logical {mustBeNumericOrLogical} = false
      ncprops.shuffle (1, 1) logical {mustBeNumericOrLogical} = true
      ncprops.deflate (1, 1) logical {mustBeNumericOrLogical} = true
      ncprops.deflateLevel (1, 1) double {mustBeNumeric} = 1
   end

   % Define the data variables. NOTE: each variable must have the same dimid.
   for v = 1:numel(vars)
      thisvar = vars{v};
      thisvid = netcdf.defVar(ncid, thisvar, xtype, dimid);

      if dochunking
         netcdf.defVarChunking(ncid, thisvid, 'CHUNKED', chunksize);
      end

      netcdf.defVarDeflate(ncid, thisvid, ...
         ncprops.shuffle, ncprops.deflate, ncprops.deflateLevel);

      if ~isempty(standardnames{v})
         netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames{v});
      end

      if ~isempty(longnames{v})
         netcdf.putAtt(ncid, thisvid, 'long_name', longnames{v});
      end

      if ~isempty(units{v})
         netcdf.putAtt(ncid, thisvid, 'units', units{v});
      end

      netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');

      % Can also use the 'comments' field to clarify the standard / long name
      % netcdf.putAtt(ncid, thisvid, 'comments', '...');
   end
end
