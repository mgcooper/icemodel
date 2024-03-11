function [info, data] = makencfile(pathdata, pathsave, simyears, ncopts)

   arguments
      pathdata (1, :) char
      pathsave (1, :) char
      simyears (:, :)
      ncopts.test_write (1, 1) logical = true
      ncopts.make_backups (1, 1) logical = true
      ncopts.time_units (1, :) char {mustBeMember(ncopts.time_units, ...
         {'hours', 'seconds'})} = 'seconds'
      ncopts.xtype (1, :) char {mustBeMember(ncopts.xtype, ...
         {'NC_FLOAT', 'NC_DOUBLE', 'NC_INT64', 'NC_UINT64', 'NC_INT', ...
         'NC_UINT', 'NC_SHORT', 'NC_USHORT', 'NC_BYTE', 'NC_UBYTE', ...
         'NC_CHAR', 'NC_STRING'})} = 'NC_FLOAT'
      ncopts.shuffle (1, 1) logical = true
      ncopts.deflate (1, 1) logical = true
      ncopts.deflateLevel (1, 1) double = 1
      ncopts.dochunking = false
   end

   % Note:
   % NC_FLOAT = single
   % NC_INT64, NC_UINT64, NC_UINT, NC_USHORT, NC_UBYTE, NC_STRING only for nc4

   % Set the file format
   netcdf.setDefaultFormat('NC_FORMAT_NETCDF4');

   % Create the output folder if it does not exist
   if ~isfolder(pathsave)
      mkdir(pathsave);
   end

   % Pull out the function options
   [xtype, shuffle, deflate, deflateLevel, dochunking] = deal(ncopts.xtype, ...
      ncopts.shuffle, ncopts.deflate, ncopts.deflateLevel, ncopts.dochunking);

   if ncopts.test_write
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
   switch ncopts.time_units
      case 'seconds'
         dt = 3600;
      case 'hours'
         dt = 1;
   end

   % FOR TESTING
   if ncopts.test_write
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
      [ncells, nlyrs, nhrs, chunksizes] = getDataDims(pathdata, thisyear, dims);

      % Read in the ice1 and ice2 data
      data = allocateDataArrays(pathdata, thisyear, ncells, nhrs, ...
         nlyrs, xtype, vars);

      % Update the time dimension
      dims.time = 0:dt:dt*(nhrs-1);
      units.dims = strrep(units.dims, 'time', ...
         sprintf('%s since %s-01-01 00:00:00 -2:00', ...
         ncopts.time_units, thisyear));

      % Create the file
      ncid = createNcFile(pathsave, thisyear, ncopts.make_backups);

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

   % Parse outputs
   % varargout = cell(nargout, 1);
   % switch nargout
   %    case 1
   %       varargout{1} = ncinfo(outfilename);
   %    case 2
   %       varargout{1} = ncinfo(outfilename);
   %       varargout{2} = ncreaddata(outfilename);
   % end
end

%%
function [ncells, nlyrs, nhrs, chunksizes] = getDataDims(pathdata, thisyear, dims)

   % Load one ice2 file to get the dimensions
   load(fullfile(pathdata, thisyear, 'ice2_1.mat'), 'ice2')

   % Set the dimensions
   ncells = numel(dims.x_easting);        % number of grid cells
   [nlyrs, ...                            % number of vertical layers
      nhrs] = size(ice2.f_liq);           % number of hours per year

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
         data.ice2.(thisvar)(:, :, n) = ice2.(thisvar); % nlyrs x nhrs x ncells
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
            data.ice1.(thisvar) = zeros(nhrs, ncells, matlabType);
         end

         for v = 1:numel(vars.ice2)
            thisvar = vars.ice2{v};
            data.ice2.(thisvar) = zeros(nlyrs, nhrs, ncells, matlabType);
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
end

%% Define netcdf file dimensions
function dimid = defineNcDimensions(ncid, ncells, nlyrs, nhrs)

   % Define the dimensions of the data arrays, in order
   dimid.cell = netcdf.defDim(ncid, 'gridcell', ncells);
   dimid_layer = netcdf.defDim(ncid, 'layer', nlyrs);
   dimid.time = netcdf.defDim(ncid, 'time', nhrs);

   dimid.ice1 = [dimid.cell dimid.time];
   dimid.ice2 = [dimid.cell dimid_layer dimid.time];
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
      if thisvar == "time"
         % Use double for time.
         thisvid = netcdf.defVar(ncid, thisvar, 'NC_DOUBLE', dimid.time);
      else
         thisvid = netcdf.defVar(ncid, thisvar, 'NC_FLOAT', dimid.cell);
      end
      netcdf.defVarDeflate(ncid, thisvid, shuffle, deflate, deflateLevel);
      netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.dims{v});
      netcdf.putAtt(ncid, thisvid, 'long_name', longnames.dims{v});
      netcdf.putAtt(ncid, thisvid, 'units', units.dims{v});
      netcdf.putAtt(ncid, thisvid, 'axis', axes.dims{v});
      % netcdf.putAtt(ncid, thisvid, '_FillValue', -9999);

      if thisvar == "elevation"
         netcdf.putAtt(ncid, thisvid, 'positive', 'up');
      end
   end

   % Define the 1d variables
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
      netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.ice1{v});
      netcdf.putAtt(ncid, thisvid, 'long_name', longnames.ice2{v});
      netcdf.putAtt(ncid, thisvid, 'units', units.ice2{v});
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
