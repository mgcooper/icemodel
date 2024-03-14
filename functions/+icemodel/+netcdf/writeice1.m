function writeice1(ncid, vars, data)
   %WRITEICE1 Write ice1 data to icemodel nc file
   %
   % WRITEICE1(NCID, VARS, DATA)
   %
   % See also:

   for v = 1:numel(vars)
      thisvar = vars{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);
      netcdf.putVar(ncid, thisvid, data.(thisvar).');

      % This writes one cell at a time. Keep this as a reminder of the access
      % pattern: [n-1 0], [1 nhrs]
      % for n = 1:ncells
      %    load(fullfile(infilepath, ['ice1_' num2str(n) '.mat']), 'ice1');
      %    netcdf.putVar(ncid, thisvid, [n-1 0], [1 nhrs], ice1.(thisvar));
      % end
   end
end
