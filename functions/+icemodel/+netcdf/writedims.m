function writedims(ncid, dims, vars)
   %WRITEDIMS Write dimensions to icemodel nc file
   %
   % WRITEDIMS(NCID, DIMS, VARS) Writes the data in DIMS.(var)
   % for each var in VARS to the file with id NCID.
   %
   % See also: icemodel.netcdf.defdimid, icemodel.netcdf.defdimvars

   % Write the grid and time dimensions
   for v = 1:numel(vars)
      thisvar = vars{v};
      try
         thisvid = netcdf.inqVarID(ncid, thisvar);
      catch e
         if contains(e.message, 'Variable not found (NC_ENOTVAR)')
            % This catches the case where dims contains "depth"
            continue
         else
            rethrow(e)
         end
      end
      netcdf.putVar(ncid, thisvid, dims.(thisvar));
   end
end
