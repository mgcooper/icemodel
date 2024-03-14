function info = makencfile(pathdata, pathsave, simmodel, forcings, ...
      userdata, simyears, opts, ncprops)

   arguments
      pathdata (1, :) char
      pathsave (1, :) char
      simmodel (1, :) char
      forcings (1, :) char
      userdata (1, :) char
      simyears (:, :)

      % Custom options
      opts.whichdata (1, :) string = ["ice1", "ice2"]
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

   % Update the API configuration
   icemodel.netcdf.config();

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
   [Z, dz, time_units, make_backups, dochunking, whichdata] = deal( ...
      opts.Z, opts.dz, opts.time_units, opts.make_backups, opts.dochunking, ...
      opts.whichdata);

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

   % Get the spatial dimension data
   dims = icemodel.netcdf.getdimensions(Z, dz);

   % Set the timestep
   switch time_units
      case 'seconds'
         dt = 3600;
      case 'hours'
         dt = 1;
   end

   % FOR TESTING
   if opts.test_write
      ncells = opts.test_numcells;
      fields = fieldnames(dims);
      for n = 1:numel(fields)
         if ~strcmp(fields{n}, {'depth', 'time'})
            dims.(fields{n}) = dims.(fields{n})(1:ncells);
         end
      end
   end

   % Create the files year by year.
   for m = 1:numyears

      thisyear = num2str(simyears(m));

      % Get the input data dimensions
      [ncells, nlyrs, nhrs, chunksizes, ...
         vars, units, longnames, standardnames] = getDataDims( ...
         pathdata, thisyear, dims, ...
         vars, units, longnames, standardnames, simmodel);

      % Update the depth
      % dims.depth = dz/2:dz:(nlyrs*opts.dz);

      % Overrule nlyrs
      nlyrs = numel(dims.depth);

      % Read in the ice1 and ice2 data
      filepath = fullfile(pathdata, thisyear);
      data = allocateDataArrays(filepath, ncells, nhrs, nlyrs, xtype, vars, ...
         whichdata, simmodel);

      % Update the time dimension
      dims.time = 0:dt:dt*(nhrs-1);
      units.dims = strrep(units.dims, 'time', ...
         sprintf('%s since %s-01-01 00:00:00', ...
         time_units, thisyear));

      % Create the file
      % <model>.<forcings>.<albedo>.<sitename>.<yyyy>.nc4
      filename = [simmodel '.' char(strjoin(whichdata, '.')) '.' ...
         forcings '.' userdata '.sw.' thisyear '.nc4'];
      ncid = icemodel.netcdf.create(fullfile(pathsave, filename), ...
         'make_backups', make_backups);
      % job = onCleanup(@() netcdf.close(ncid));

      % Define dimensions and variables
      defineNcDimsAndVars(ncid, ncells, nlyrs, nhrs, vars, units, axes, ...
         longnames, standardnames, chunksizes, xtype, shuffle, deflate, ...
         deflateLevel, dochunking, whichdata)

      % Close definition mode
      netcdf.endDef(ncid)

      % Write the data
      writeNcData(ncid, data, dims, vars, whichdata, filepath)

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
      vars, units, longnames, standardnames, simmodel)

   % Load one ice1 file to get the dimensions
   ice1 = load(fullfile(pathdata, thisyear, 'ice1_1.mat')).('ice1');
   ice2 = load(fullfile(pathdata, thisyear, 'ice2_1.mat')).('ice2');

   % Remove variables from the dictionaries that are not present in ice1/2
   ivars1 = ismember(vars.ice1, ice1.Properties.VariableNames);
   vars.ice1 = vars.ice1(ivars1);
   units.ice1 = units.ice1(ivars1);
   longnames.ice1 = longnames.ice1(ivars1);
   standardnames.ice1 = standardnames.ice1(ivars1);

   try
      ivars2 = ismember(vars.ice2, fieldnames(ice2));
      vars.ice2 = vars.ice2(ivars2);
      units.ice2 = units.ice2(ivars2);
      longnames.ice2 = longnames.ice2(ivars2);
      standardnames.ice2 = standardnames.ice2(ivars2);
   catch
   end

   % Set the dimensions
   ncells = numel(dims.x_easting);        % number of grid cells
   [nlyrs, ...                            % number of vertical layers
      nhrs] = size(ice2.Tice);            % number of hours per year

   % Define chunkSize based on data access patterns. Larger chunk sizes
   % increase memory usage during read/write.
   % chunksizes.ice1 = ceil([ncells/2, nhrs/2]);  % all cells, annual chunks
   % chunksizes.ice2 = ceil([1, nlyrs, nhrs/12]); % one cell, all layers, monthly chunks

   chunksizes.ice1 = [ncells, nhrs];      % all cells, annual chunks
   chunksizes.ice2 = [1, nlyrs, nhrs];    % one cell, all layers, annual chunks
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
function data = allocateDataArrays(filepath, ncells, nhrs, ...
      nlyrs, xtype, vars, whichdata, simmodel)

   % Note: allocate the data using matlab column-major format. Transpose the
   % data to row-major when writing to netcdf.

   data = preallocateDataArrays(vars, ncells, nhrs, nlyrs, xtype, 'ice1');

   % Read in the ice1 data and fill the arrays
   if ismember('ice1', whichdata)
      for n = 1:ncells
         ice1 = load(fullfile(filepath, ['ice1_' num2str(n) '.mat'])).('ice1');

         % Set freeze zero for skinmodel (freeze is Qf, not refreeze)
         if strcmp('skinmodel', simmodel) && isvariable('freeze', ice1)
            ice1.freeze = 0 * ice1.freeze;
         end

         for v = 1:numel(vars.ice1)
            thisvar = vars.ice1{v};
            data.ice1.(thisvar)(:, n) = ice1.(thisvar); % nhrs x ncells
         end
      end
   end

   % Read in the ice2 data and fill the arrays
   if ismember('ice2', whichdata)
      for n = 1:ncells
         ice2 = load(fullfile(filepath, ['ice2_' num2str(n) '.mat'])).('ice2');
         for v = 1:numel(vars.ice2)
            thisvar = vars.ice2{v};
            data.ice2.(thisvar)(1:size(ice2.(thisvar), 1), :, n) = ice2.(thisvar); % nlyrs x nhrs x ncells
         end
      end
   end
end

function data = preallocateDataArrays(vars, ncells, nhrs, nlyrs, xtype, whichdata)

   matlabType = nctype2mat(xtype);

   % Preallocate based on the data type
   switch matlabType
      case {'char', 'string', 'cell'}
         error( ...
            'Preallocation for %s is not supported in this function.', xtype);
      otherwise
         % Preallocate data arrays based on specified dimensions and data type

         if ismember('ice1', whichdata)
            for v = 1:numel(vars.ice1)
               thisvar = vars.ice1{v};
               data.ice1.(thisvar) = nan(nhrs, ncells, matlabType);
            end
         end

         if ismember('ice2', whichdata)
            for v = 1:numel(vars.ice2)
               thisvar = vars.ice2{v};
               data.ice2.(thisvar) = nan(nlyrs, nhrs, ncells, matlabType);
            end
         end
   end

   % The purpose of this is to return both data structs in one to later access
   % it as data.('ice1').(thisvar) and data.('ice2').(thisvar) while also being
   % able to just get ice1/ice2 out of this function. To use this, remove the
   % data. prefix on ice1 and ice2 above. However, this might lead to confusion,
   % so I reverted to the original output where it's data.ice1, or data.ice2, or
   % data.ice1 and data.ice2.
   % switch nargout
   %    case 1
   %       if all(ismember({'ice1', 'ice2'}, whichdata))
   %          data.ice1 = ice1;
   %          data.ice2 = ice2;
   %          varargout{1} = data;
   %
   %       elseif ismember('ice1', whichdata)
   %          varargout{1} = ice1;
   %
   %       elseif ismember('ice2', whichdata)
   %          varargout{1} = ice2;
   %       end
   %
   %    case 2
   %       varargout{1} = ice1;
   %       varargout{2} = ice2;
   % end
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
      deflateLevel, dochunking, whichdata)

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

      % standard names
      if ~isempty(standardnames.dims{v})
         netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.dims{v});
      end

      % long names
      if ~isempty(longnames.dims{v})
         netcdf.putAtt(ncid, thisvid, 'long_name', longnames.dims{v});
      end

      % Units
      if ~isempty(units.dims{v})
         netcdf.putAtt(ncid, thisvid, 'units', units.dims{v});
      end

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
   if ismember('ice1', whichdata)
      for v = 1:numel(vars.ice1)
         thisvar = vars.ice1{v};
         thisvid = netcdf.defVar(ncid, thisvar, xtype, dimid.ice1);

         if dochunking
            netcdf.defVarChunking(ncid, thisvid, 'CHUNKED', chunksize.ice1);
         end
         netcdf.defVarDeflate(ncid, thisvid, shuffle, deflate, deflateLevel);

         if ~isempty(standardnames.ice1{v})
            netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.ice1{v});
         end
         if ~isempty(longnames.ice1{v})
            netcdf.putAtt(ncid, thisvid, 'long_name', longnames.ice1{v});
         end
         if ~isempty(units.ice1{v})
            netcdf.putAtt(ncid, thisvid, 'units', units.ice1{v});
         end

         netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');

         % Can also use the 'comments' field to clarify the standard / long name
         % netcdf.putAtt(ncid, thisvid, 'comments', '...');
      end
   end

   % Define the 2d x time variables
   if ismember('ice2', whichdata)
      for v = 1:numel(vars.ice2)
         thisvar = vars.ice2{v};
         thisvid = netcdf.defVar(ncid, thisvar, xtype, dimid.ice2);

         if dochunking
            netcdf.defVarChunking(ncid, thisvid, 'CHUNKED', chunksize.ice2);
         end
         netcdf.defVarDeflate(ncid, thisvid, shuffle, deflate, deflateLevel);

         if ~isempty(standardnames.ice2{v})
            netcdf.putAtt(ncid, thisvid, 'standard_name', standardnames.ice2{v});
         end
         if ~isempty(longnames.ice2{v})
            netcdf.putAtt(ncid, thisvid, 'long_name', longnames.ice2{v});
         end
         if ~isempty(units.ice2{v})
            netcdf.putAtt(ncid, thisvid, 'units', units.ice2{v});
         end

         netcdf.putAtt(ncid, thisvid, 'coordinates', 'longitude latitude');
      end
   end
end
%%
function writeNcData(ncid, data, dims, vars, whichdata, filepath)

   % Write the grid and time dimensions
   for v = 1:numel(vars.dims)
      thisvar = vars.dims{v};
      thisvid = netcdf.inqVarID(ncid, thisvar);
      netcdf.putVar(ncid, thisvid, dims.(thisvar));
   end

   % Write the 1d variables
   if ismember('ice1', whichdata)
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
   end

   % Write the 2d variables
   if ismember('ice2', whichdata)
      for v = 1:numel(vars.ice2)
         thisvar = vars.ice2{v};
         thisvid = netcdf.inqVarID(ncid, thisvar);

         % NOTE: The methods which use start, count should work if
         % size(ice2.(thisvar)) is used for count, but also recall this method
         % from the allocateData subfunction:
         %
         % data.ice2.(thisvar)(1:size(ice2.(thisvar), 1), :, n) = ice2.(thisvar); % nlyrs x nhrs x ncells

         % This writes all data at once, which does not work with 16 GB ram
         % netcdf.putVar(ncid, thisvid, permute(data.ice2.(thisvar), [3 1 2]));

         % To preallocate in chunks:
         % chunksize = 10;
         % for n = 1:chunksize:numel(dims.gridcell)
         %    [s, e] = chunkLoopInds(n, 1, chunksize);
         %    tmp = preallocateDataArrays(vars, chunksize, 8760, 500, 'NC_FLOAT', 'ice2');
         %    for m = s:e
         %       tmp.ice2.(thisvar)(:, :, m) = ...
         %          load(fullfile(filepath, ['ice2_' num2str(m) '.mat'])).('ice2').(thisvar);
         %    end
         %    % The permute rearranges [depth time gridcell] to [gridcell depth time]
         %    netcdf.putVar(ncid, thisvid, [s-1 0 0], size(permute(tmp.ice2.(thisvar), [3 1 2])), ...
         %       permute(tmp.ice2.(thisvar), [3 1 2]));
         % end

         % This writes one cell at a time. Keep this as a reminder of the access
         % pattern: [n-1 0 0], [1 nlyrs nhrs]. Make sure the [start], [count]
         % matches the dimid and the data orientation. Here [n-1 0 0], [1 nz nt]
         % where nz is number of depth layers and nt number of time slices. The
         % netcdf4 data model might allow the shape of the data sent in to
         % differ, e.g., the putVar step works when ice2.(thisvar) is
         % transposed, but the data is scrambled.
         for n = 1:numel(dims.gridcell)
            load(fullfile(filepath, ['ice2_' num2str(n) '.mat']), 'ice2');
            netcdf.putVar(ncid, thisvid, [n-1 0 0], [1 size(ice2.(thisvar))], ice2.(thisvar));
         end
      end
   end

   % This is wrong. I thought it explained why the transpose was necessary, but
   % in fact the transpose was wrong.
   % recall: [gridcell x depth x time] becomes [time x depth x gridcell]
   % so [n-1 0 0], [1 500 8760] becomes [0 0 n-1], [8760 500 1], thus the data,
   % which is [500 x 8760], must be transposed to become [8760 x 500], which,
   % when appended to the n-1, 0 start count, yields [0 0 n-1], [8760 500 1],
   % which is an ideal access pattern in a column-major layout, but worst
   % pattern in row-major, but it depends on whether the right way to think of
   % it is in terms of how the dims are defined: [gridcell x depth x time] =
   % [1 500 8760] (ideal for row major) or how it goes in when the data is
   % transposed / written to disk: [8760 500 1] (ideal for column major).
   %
   % But this seems to contradict the full array case:
   % the data is [500 8760 2479], which is permuted to become [2479 500 8760]
   % the data is [nlyrs nhrs ncell] which is permuted to become [ncell nlyrs nhrs]
   % thus permuation makes the data match the dimid: [gridcell depth time]
   %
   % if it's written directly to the dimensions defined in dimid:
end
