function info = makencfile(pathdata, pathsave, simyears, opts, ncprops)

   arguments
      pathdata (1, :) char
      pathsave (1, :) char
      simyears (:, :)

      % Custom options
      opts.dz (1, 1) = 0.04
      opts.Z (1, 1) = 20;
      opts.test_write (1, 1) logical = true
      opts.test_numcells (1, 1) = 10
      opts.make_backups (1, 1) logical = true
      opts.time_units (1, :) char {mustBeMember(opts.time_units, ...
         {'hours', 'seconds'})} = 'seconds'
      opts.dochunking = false

      % Netcdf api options
      ncprops.format (1, :) char {mustBeMember(ncprops.format, ...
         {'NC_FORMAT_CLASSIC', 'NC_FORMAT_64BIT', ...
         'NC_FORMAT_NETCDF4', 'NC_FORMAT_NETCDF4_CLASSIC'})} ...
         = getenv('ICEMODEL_NC_DEFAULT_FORMAT')

      ncprops.xtype (1, :) char {mustBeMember(ncprops.xtype, ...
         {'NC_FLOAT', 'NC_DOUBLE', 'NC_INT64', 'NC_UINT64', 'NC_INT', ...
         'NC_UINT', 'NC_SHORT', 'NC_USHORT', 'NC_BYTE', 'NC_UBYTE', ...
         'NC_CHAR', 'NC_STRING'})} ...
         = 'NC_DOUBLE'
      ncprops.shuffle (1, 1) logical = true
      ncprops.deflate (1, 1) logical = true
      ncprops.deflateLevel (1, 1) double = 1
   end

   % Note:
   % NC_FLOAT = single
   % NC_INT64, NC_UINT64, NC_UINT, NC_USHORT, NC_UBYTE, NC_STRING only for nc4

   % Set the file format
   oldformat = netcdf.setDefaultFormat(ncprops.format);

   % Create the output folder if it does not exist
   if ~isfolder(pathsave)
      mkdir(pathsave);
   end

   % Pull out the netcdf api options
   [xtype, shuffle, deflate, deflateLevel] = deal( ...
      ncprops.xtype, ncprops.shuffle, ncprops.deflate, ncprops.deflateLevel);

   % Pull out the function options
   [Z, dz, time_units, make_backups, dochunking] = deal( ...
      opts.Z, opts.dz, opts.time_units, opts.make_backups, opts.dochunking);

   if opts.test_write
      numyears = 1;
   else
      numyears = numel(simyears);
   end

   % Define variables and attributes.
   for v = ["dims", "ice1", "ice2"]
      [vars.(v), longnames.(v), units.(v), axes.(v), standardnames.(v)] ...
         = icemodel.netcdf.getdefaults(v, ...
         {'varnames', 'longnames', 'units', 'axes', 'standardnames'});
   end

   % Load the grid coordinates and elevation.
   [dims.x_easting, dims.y_northing, dims.elevation, ...
      dims.latitude, dims.longitude] = loadIceMask("oldmask", ...
      varnames=["X", "Y", "Elev", "Lat", "Lon"]);

   % Set the timestep
   switch time_units
      case 'seconds'
         dt = 3600;
      case 'hours'
         dt = 1;
   end

   % FOR TESTING
   if opts.test_write
      ncells = 3;
      [dims.x_easting, dims.y_northing, dims.elevation, ...
         dims.latitude, dims.longitude] = deal( ...
         dims.x_easting(1:ncells), dims.y_northing(1:ncells), ...
         dims.elevation(1:ncells), dims.latitude(1:ncells), ...
         dims.longitude(1:ncells));
   end

   % Create the files year by year.
   for m = 1:numyears

      thisyear = num2str(simyears(m));

      % Get the input data dimensions
      [ncells, nlyrs, nhrs, chunksizes, ...
         vars, units, longnames, standardnames] = getDataDims( ...
         pathdata, thisyear, dims, ...
         vars, units, longnames, standardnames);

      % Update the depth
      % dims.depth = dz/2:dz:(nlyrs*opts.dz);
      dims.depth = dz/2:dz:Z;

      % Overrule nlyrs
      nlyrs = numel(dims.depth);

      % Read in the ice1 and ice2 data
      data = allocateDataArrays(pathdata, thisyear, ncells, nhrs, nlyrs, xtype, vars);

      % Update the time dimension
      dims.time = 0:dt:dt*(nhrs-1);
      units.dims = strrep(units.dims, 'time', ...
         sprintf('%s since %s-01-01 00:00:00', ...
         time_units, thisyear));

      % Create the file
      ncid = createNcFile(pathsave, thisyear, make_backups);

      % Define dimensions and variables
      defineNcDimsAndVars(ncid, ncells, nlyrs, nhrs, vars, units, axes, ...
         longnames, standardnames, chunksizes, xtype, shuffle, deflate, ...
         deflateLevel, dochunking)

      % Close definition mode
      netcdf.endDef(ncid)

      % Write the data
      writeNcData(ncid, data, dims, vars)

      % Close the file
      netcdf.close(ncid);
   end

   % Reset the file format
   netcdf.setDefaultFormat(oldformat);

   % Parse outputs
   if nargout > 1
      info = ncinfo(outfilename);
   end
end

%%
% function [ncells, nlyrs, nhrs, chunksizes] = getDataDims(pathdata, thisyear, dims)
function [ncells, nlyrs, nhrs, chunksizes, ...
      vars, units, longnames, standardnames] = getDataDims( ...
      pathdata, thisyear, dims, ...
      vars, units, longnames, standardnames)

   % Load one ice2 file to get the dimensions
   ice1 = load(fullfile(pathdata, thisyear, 'ice1_1.mat')).('ice1');
   ice2 = load(fullfile(pathdata, thisyear, 'ice2_1.mat')).('ice2');

   % Remove variables from the dictionaries that are not present in ice1/2
   ivars1 = ismember(vars.ice1, ice1.Properties.VariableNames);
   ivars2 = ismember(vars.ice2, fieldnames(ice2));

   vars.ice1 = vars.ice1(ivars1);
   vars.ice2 = vars.ice2(ivars2);
   units.ice1 = units.ice1(ivars1);
   units.ice2 = units.ice2(ivars2);
   longnames.ice1 = longnames.ice1(ivars1);
   longnames.ice2 = longnames.ice2(ivars2);
   standardnames.ice1 = standardnames.ice1(ivars1);
   standardnames.ice2 = standardnames.ice2(ivars2);

   % Set the dimensions
   ncells = numel(dims.x_easting);        % number of grid cells
   [nlyrs, ...                            % number of vertical layers
      nhrs] = size(ice2.Tice);            % number of hours per year

   % Define chunkSize based on data access patterns. Larger chunk sizes
   % increase memory usage during read/write.
   chunksizes.ice1 = ceil([ncells/2, nhrs/2]);  % all cells, annual chunks
   chunksizes.ice2 = ceil([1, nlyrs, nhrs/12]); % one cell, all layers, annual chunks

   % chunksizes.ice1 = [ncells, nhrs];      % all cells, annual chunks
   % chunksizes.ice2 = [1, nlyrs, nhrs];    % one cell, all layers, annual chunks
   % chunksizes.ice1 = [ncells, 24];      % all cells, daily chunks
   % chunksizes.ice2 = [1, nlyrs, 24];    % one cell, all layers, daily chunks

   % default chunking was set to:
   % [744,4392] for [1487,8784] data (ncells x nhrs)
   % [124,42,732] for [1487,500,8784] data (ncells x nlayers x ndays)

   % For the 3D data consider [X, 500, Y] where X represents a typical subset
   % size of grid cells (could be a smaller number if you access grid cells
   % individually or in small groups), and Y is a large fraction of nhrs that
   % represents typical time period access patterns (e.g., several months).

   % Filling the arrays must be done column-major, but what about writing them?
   % When writing, the primary concern isn't the order of operations but
   % ensuring the data layout in the file matches how it will be accessed later,
   % by other programs which assume row-major.
   %
   % When the entire array is written at once, as it is here, the efficiency
   % concern shifts from the order of writes to ensuring that the chunking in
   % the NetCDF file matches your typical access patterns. If data is
   % predominantly accessed in large contiguous blocks, having the data and
   % chunks aligned to those access patterns is ideal.
end

%%
function data = allocateDataArrays(pathdata, thisyear, ncells, nhrs, ...
      nlyrs, xtype, vars)

   % Note: allocate the data using matlab column-major format. Transpose the
   % data to row-major when writing to netcdf.

   % Preallocate data arrays based on specified dimensions and data type
   data = preallocateDataArrays(vars, ncells, nhrs, nlyrs, xtype);

   % Read in the ice1 and ice2 data and fill the arrays
   filepath = fullfile(pathdata, thisyear);
   for n = 1:ncells
      ice1 = load(fullfile(filepath, ['ice1_' num2str(n) '.mat'])).('ice1');
      ice2 = load(fullfile(filepath, ['ice2_' num2str(n) '.mat'])).('ice2');

      for v = 1:numel(vars.ice1)
         thisvar = vars.ice1{v};
         data.ice1.(thisvar)(:, n) = ice1.(thisvar); % nhrs x ncells
      end

      for v = 1:numel(vars.ice2)
         thisvar = vars.ice2{v};
         data.ice2.(thisvar)(1:size(ice2.(thisvar), 1), :, n) = ice2.(thisvar); % nlyrs x nhrs x ncells
      end
   end
end
function data = preallocateDataArrays(vars, ncells, nhrs, nlyrs, xtype)

   matlabType = nctype2mat(xtype);

   % Preallocate based on the data type
   switch matlabType
      case {'char', 'string', 'cell'}
         error( ...
            'Preallocation for %s is not supported in this function.', xtype);
      otherwise

         for v = 1:numel(vars.ice1)
            thisvar = vars.ice1{v};
            data.ice1.(thisvar) = nan(nhrs, ncells, matlabType);
         end

         for v = 1:numel(vars.ice2)
            thisvar = vars.ice2{v};
            data.ice2.(thisvar) = nan(nlyrs, nhrs, ncells, matlabType);
         end
   end
end
%%
function ncid = createNcFile(pathsave, thisyear, make_backups)

   % Set the output netcdf file name
   filename = fullfile(pathsave, ['icemodel_' thisyear '.nc4']);

   % Back up the file if it exists
   if isfile(filename)
      backupfile(filename, make_backups);
      fprintf('Moving to recycle bin: %s \n', filename)
      status = recycle;
      recycle("on");
      delete(filename)
      recycle(status);
   end

   % Create the netcdf file.
   ncid = netcdf.create(filename, 'NETCDF4');

   % Set NOFILL to avoid writing fill values that are later replaced by data.
   netcdf.setFill(ncid, 'NC_NOFILL');

   % Add user metadata
   varid = netcdf.getConstant('GLOBAL');
   netcdf.putAtt(ncid, varid, 'model', ['IceModel ' getenv('ICEMODEL_VERSION')]);
   netcdf.putAtt(ncid, varid, 'created_by', getenv('USER'));
   netcdf.putAtt(ncid, varid, 'date_created', ...
      char(datetime("now", "Format", "dd-MMM-uuuu hh:mm:ss", "TimeZone", "UTC")));
   netcdf.putAtt(ncid, varid, 'contact', getenv('ICEMODEL_CONTACT'));
end

%% Define netcdf file dimensions
function dimid = defineNcDimensions(ncid, ncells, nlyrs, nhrs)

   % Define the dimensions of the data arrays, in order. UPDATE: a "coordinate
   % variable" is a "1-d variable with the same name as its dimension [e.g.,
   % time(time) ], and it is defined as a numeric data type with values that are
   % ordered monotonically". Thus, each dimension needs a dimid.
   dimid.gridcell = netcdf.defDim(ncid, 'gridcell', ncells);
   dimid.depth = netcdf.defDim(ncid, 'depth', nlyrs);
   dimid.time = netcdf.defDim(ncid, 'time', nhrs);

   % The first dimension in matlab is the last in netcdf. Thus ice1 will be:
   % time x gridcell, and ice2 will be time x depth x gridcell.
   dimid.ice1 = [dimid.gridcell dimid.time];
   dimid.ice2 = [dimid.gridcell dimid.depth dimid.time];

   % Note: for a lat lon or x y, it would be:
   % dimid.ice1 = [dimid.lon dimid.lat dimid.time];
   % dimid.ice2 = [dimid.lon dimid.lat dimid.depth dimid.time];
   %
   % As explained here:
   %
   % If any or all of the dimensions of a variable have the interpretations of
   % "date or time" (T), "height or depth" (Z), "latitude" (Y), or "longitude"
   % (X) then we recommend, but do not require (see Section 1.5, "Relationship
   % to the COARDS Conventions"), those dimensions to appear in the relative
   % order T, then Z, then Y, then X in the CDL definition corresponding to the
   % file. All other dimensions should, whenever possible, be placed to the left
   % of the spatiotemporal dimensions.

   % NOTE: This suggests I should reverse the ordering above. However, the
   % actual data is reversed, see ncdump -h, ordering is (time, gridcell) and
   % (time, depth, gridcell)
   %
   % Also note:
   %
   % "we allow but do not require the units attribute of dimensionless vertical
   % coordinates to take the values "level", "layer", or "sigma_level.""
   %
   % But I was unable to find a similar "units" (or "standard_name") for grid
   % cell index

end

%% Define netcdf file variables and attributes
function defineNcDimsAndVars(ncid, ncells, nlyrs, nhrs, vars, units, axes, ...
      longnames, standardnames, chunksize, xtype, shuffle, deflate, ...
      deflateLevel, dochunking)

   % Define the dimensions IDs.
   dimid = defineNcDimensions(ncid, ncells, nlyrs, nhrs);

   % Define the grid and time dimensions and attributes.
   for v = 1:numel(vars.dims)
      thisvar = vars.dims{v};

      % Explicitly assign the fundamental dimensions. Their
      switch thisvar

         case "gridcell"
            thisvid = netcdf.defVar(ncid, thisvar, 'NC_INT', dimid.gridcell);

         case "depth"
            thisvid = netcdf.defVar(ncid, thisvar, 'NC_FLOAT', dimid.depth);

         case "time"
            thisvid = netcdf.defVar(ncid, thisvar, 'NC_DOUBLE', dimid.time);

         otherwise % lat, lon, x, y, elev
            thisvid = netcdf.defVar(ncid, thisvar, 'NC_FLOAT', dimid.gridcell);
      end

      netcdf.defVarDeflate(ncid, thisvid, shuffle, deflate, deflateLevel);
      netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.dims{v});
      netcdf.putAtt(ncid, thisvid, 'long_name', longnames.dims{v});
      netcdf.putAtt(ncid, thisvid, 'units', units.dims{v});
      % netcdf.putAtt(ncid, thisvid, '_FillValue', -9999);

      % Axis
      if ~isempty(axes.dims{v})
         netcdf.putAtt(ncid, thisvid, 'axis', axes.dims{v});
      end

      % Additional attributes e.g. coordinates (linked to auxiliary coordinate
      % variables), positive direction for vertical coordinates, calendar, etc.
      switch thisvar

         case "gridcell"
            netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');
            netcdf.putAtt(ncid, thisvid, 'comment', 'A unique identifier for each grid cell');

         case ["x_easting", "y_northing"]

            netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');

         case "elevation"
            netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');
            netcdf.putAtt(ncid, thisvid, 'positive', 'up');

         case "depth"
            netcdf.putAtt(ncid, thisvid, 'positive', 'down');

         case "time"
            netcdf.putAtt(ncid, thisvid, 'calendar', 'noleap');

         otherwise
            % lat, lon - do not add the coordinates
      end
   end

   % Define the 1d x time variables
   for v = 1:numel(vars.ice1)
      thisvar = vars.ice1{v};
      thisvid = netcdf.defVar(ncid, thisvar, xtype, dimid.ice1);

      if dochunking
         netcdf.defVarChunking(ncid, thisvid, 'CHUNKED', chunksize.ice1);
      end
      netcdf.defVarDeflate(ncid, thisvid, shuffle, deflate, deflateLevel);
      netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.ice1{v});
      netcdf.putAtt(ncid, thisvid, 'long_name', longnames.ice1{v});
      netcdf.putAtt(ncid, thisvid, 'units', units.ice1{v});

      netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');

      % Can also use the 'comments' field to clarify the standard / long name
      % netcdf.putAtt(ncid, thisvid, 'comments', '...');
   end

   % Define the 2d variables
   for v = 1:numel(vars.ice2)
      thisvar = vars.ice2{v};
      thisvid = netcdf.defVar(ncid, thisvar, xtype, dimid.ice2);

      if dochunking
         netcdf.defVarChunking(ncid, thisvid, 'CHUNKED', chunksize.ice2);
      end
      netcdf.defVarDeflate(ncid, thisvid, shuffle, deflate, deflateLevel);
      netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.ice2{v});
      netcdf.putAtt(ncid, thisvid, 'long_name', longnames.ice2{v});
      netcdf.putAtt(ncid, thisvid, 'units', units.ice2{v});

      netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');
   end
end
%%
function writeNcData(ncid, data, dims, vars)

   % Write the grid and time dimensions
   for v = 1:numel(vars.dims)
      thisvar = vars.dims{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);
      netcdf.putVar(ncid, thisvid, dims.(thisvar));
   end

   % Write the 1d variables
   for v = 1:numel(vars.ice1)
      thisvar = vars.ice1{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);
      netcdf.putVar(ncid, thisvid, data.ice1.(thisvar).');

      % This writes one cell at a time. Keep this as a reminder of the access
      % pattern: [n-1 0], [1 nhrs]
      % for n = 1:ncells
      %    load(fullfile(infilepath, ['ice1_' num2str(n) '.mat']), 'ice1');
      %    netcdf.putVar(ncid, thisvid, [n-1 0], [1 nhrs], ice1.(thisvar));
      % end
   end

   % Write the 2d variables
   for v = 1:numel(vars.ice2)
      thisvar = vars.ice2{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);
      netcdf.putVar(ncid, thisvid, permute(data.ice2.(thisvar), [3 1 2]));

      % This writes one cell at a time. Keep this as a reminder of the access
      % pattern: [n-1 0 0], [1 nlyrs nhrs]
      % for n = 1:ncells
      %    load(fullfile(infilepath, ['ice2_' num2str(n) '.mat']), 'ice2');
      %    netcdf.putVar(ncid, thisvid, [n-1 0 0], [1 nlyrs nhrs], ice2.(thisvar)');
      % end
   end
end
