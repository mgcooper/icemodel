function writedims(ncid, dims, vars)
   %WRITEDIMS Write dimensions to icemodel nc file
   %
   % WRITEDIMS(NCID, DIMS, VARS)
   %
   % See also: icemodel.netcdf.defdimid, icemodel.netcdf.defdimvars

   % Write the grid and time dimensions
   for v = 1:numel(vars)
      thisvar = vars{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);
      netcdf.putVar(ncid, thisvid, dims.(thisvar));
   end
end
