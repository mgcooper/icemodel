function redefatt(filelist, varname, newatts)

   arguments
      filelist (1, :) cell
      varname (1, :) cell
      newatts (1, :) cell
   end

   for n = 1:numel(filelist)

      f = filelist{n};

      % Open the NetCDF file for writing
      fid = netcdf.open(f, 'WRITE');

      % Get the ID of the variable
      vid = netcdf.inqVarID(fid, varname);

      % Update the attribute of the variable
      netcdf.putAtt(fid, vid, 'units', newatts{n});

      % How I originally used this:
      % YYYY = num2str(allyears(n));
      % correctTimeUnits = ['seconds since ' YYYY '-01-01 00:00:00'];
      % netcdf.putAtt(fid, vid, 'units', correctTimeUnits);

      % Close the NetCDF file
      netcdf.close(fid);
   end
end
