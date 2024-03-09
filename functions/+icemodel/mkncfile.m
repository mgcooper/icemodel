function msg = mkncfile(filename, X, Y, ID, varnames, varargin)
   %MKNCFILE Make icemodel nc file

   % Input data is saved in varargin
   data = varargin;

   % Input validation
   if nargin < 6
      error( ...
         ['Invalid number of input arguments. ' ...
         'Expected 5 inputs: filename, X, Y, ID, and varnames.']);
   end
   % if numel(X) ~= numel(Y) || numel(X) ~= numel(ID) || numel(X) ~= size(TWD,1)
   %    error('The input arrays X, Y, ID, and TWD must have the same number of elements organized columnwise.');
   % end
   if ~(ischar(filename) || iscellstr(filename) || isstring(filename))
      error('Input filename must be a char, cellstr array, or string array.');
   end

   % Ensure filename is a cell array for iteration
   if ischar(filename) || isstring(filename)
      filename = cellstr(filename);
   end

   % Ensure the number of filenames matches the number of columns of TWD
   assert(numel(filename) == size(TWD, 2),'Number of filenames must equal the number of columns of TWD')

   % Determine the number of grid cells
   gridcell = numel(X);

   % Create the NetCDF files
   for n = 1:numel(filename)
      try
         ncid = netcdf.create(filename{n}, 'CLOBBER');

         % Define dimensions
         dimid_gridcell = netcdf.defDim(ncid, 'gridcell', gridcell);

         % Define variables
         varid_ID = netcdf.defVar(ncid, 'ID', 'double', dimid_gridcell);
         varid_lat = netcdf.defVar(ncid, 'lat', 'double', dimid_gridcell);
         varid_lon = netcdf.defVar(ncid, 'lon', 'double', dimid_gridcell);
         varid_x = netcdf.defVar(ncid, 'x', 'double', dimid_gridcell);
         varid_y = netcdf.defVar(ncid, 'y', 'double', dimid_gridcell);

         % Define variable attributes
         netcdf.putAtt(ncid, varid_lat, 'long_name', 'latitude');
         netcdf.putAtt(ncid, varid_lat, 'units', 'degrees_north');
         netcdf.putAtt(ncid, varid_lat, '_FillValue', -9999);
         netcdf.putAtt(ncid, varid_lat, 'standard_name', 'latitude');
         netcdf.putAtt(ncid, varid_lat, 'axis', 'Y');

         netcdf.putAtt(ncid, varid_lon, 'long_name', 'longitude');
         netcdf.putAtt(ncid, varid_lon, 'units', 'degrees_east');
         netcdf.putAtt(ncid, varid_lon, '_FillValue', -9999);
         netcdf.putAtt(ncid, varid_lon, 'standard_name', 'longitude');
         netcdf.putAtt(ncid, varid_lon, 'axis', 'X');

         % Attributes for X coordinate (eastings)
         netcdf.putAtt(ncid, varid_x, 'long_name', 'x coordinate of projection');
         netcdf.putAtt(ncid, varid_x, 'units', 'meters');
         netcdf.putAtt(ncid, varid_x, '_FillValue', -9999);
         netcdf.putAtt(ncid, varid_x, 'standard_name', 'projection_x_coordinate');
         netcdf.putAtt(ncid, varid_x, 'axis', 'X');

         % Attributes for Y coordinate (northings)
         netcdf.putAtt(ncid, varid_y, 'long_name', 'y coordinate of projection');
         netcdf.putAtt(ncid, varid_y, 'units', 'meters');
         netcdf.putAtt(ncid, varid_y, '_FillValue', -9999);
         netcdf.putAtt(ncid, varid_y, 'standard_name', 'projection_y_coordinate');
         netcdf.putAtt(ncid, varid_y, 'axis', 'Y');

         % Attributes for melt
         varid_melt = netcdf.defVar(ncid, 'melt', 'double', dimid_gridcell);
         netcdf.putAtt(ncid, varid_melt, 'units', 'm/hr');
         netcdf.putAtt(ncid, varid_melt, 'long_name', 'melt water production');
         netcdf.putAtt(ncid, varid_melt, '_FillValue', -9999.0);

         % Add user metadata
         varid = netcdf.getConstant('GLOBAL');
         netcdf.putAtt(ncid,varid,'Created_by' ,getenv('USER'));
         netcdf.putAtt(ncid,varid,'Created_on' ,char(datetime("now","Format","dd-MMM-uuuu hh:mm:ss")));

         % End "define" mode and enter "data" mode
         netcdf.endDef(ncid);

         % Write data to variables
         netcdf.putVar(ncid, varid_lat, Y);
         netcdf.putVar(ncid, varid_lon, X);
         netcdf.putVar(ncid, varid_melt, TWD(:,n));
         netcdf.putVar(ncid, varid_ID, ID);

         % Close the NetCDF file
         netcdf.close(ncid);

         % Return a success message
         msg = sprintf('NetCDF file successfully created: %s', filename{n});
      catch ME
         % If an error occurs, close the NetCDF file, delete it, and re-throw the error
         try
            netcdf.close(ncid);
            delete(filename{n});
         catch
            % Do nothing if an error occurs while closing the NetCDF file and deleting it
         end
         rethrow(ME);
      end
   end

end
