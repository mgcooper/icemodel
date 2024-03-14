function defdatavars(ncid, dimid, vars, standardnames, longnames, ...
      units, chunksize, xtype, dochunking, ncprops)
   %DEFDATAVARS Define the icemodel data variables and attributes.
   %
   %  DEFDATAVARS(NCID, DIMID, VARS, STANDARDNAMES, LONGNAMES, UNITS, AXES, NCPROPS)
   %
   % See also:

   arguments
      ncid (1, 1)
      dimid
      vars
      standardnames
      longnames
      units
      chunksize
      xtype
      dochunking = false
      ncprops.shuffle (1, 1) logical = true
      ncprops.deflate (1, 1) logical = true
      ncprops.deflateLevel (1, 1) double = 1
   end

   % Define the data variables. NOTE: each variable must have the same dimid
   % i.e., the same dimensions. See defdimid.
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
