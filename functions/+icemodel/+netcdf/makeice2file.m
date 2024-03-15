function info = makencfile(pathdata, pathsave, simmodel, forcings, ...
      userdata, sitename, simyears, opts, ncprops)

   arguments
      pathdata (1, :) char
      pathsave (1, :) char
      simmodel (1, :) char
      forcings (1, :) char
      userdata (1, :) char
      sitename (1, :) char
      simyears (:, :)

      % Custom options
      opts.dz (1, 1) = 0.04
      opts.Z (1, 1) = 20;
      opts.whichdata (1, :) char {mustBeMember(opts.whichdata, ...
         {'ice1', 'ice2', 'met'})} ...
         = 'ice1'

      opts.timeunits (1, :) char {mustBeMember(opts.timeunits, ...
         {'hours', 'seconds'})} ...
         = 'seconds'

      opts.testwrite (1, 1) logical = true
      opts.numcells (1, 1) = 10
      opts.makebackups (1, 1) logical = true
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
   [whichdata, Z, dz, timeunits, makebackups, dochunking] = deal( ...
      opts.whichdata, opts.Z, opts.dz, opts.timeunits, opts.makebackups, ...
      opts.dochunking);

   if opts.testwrite
      numyears = 1;
   else
      numyears = numel(simyears);
   end

   % Define variables and attributes.
   for v = ["dims", whichdata]
      [vars.(v), longnames.(v), units.(v), axes.(v), standardnames.(v)] ...
         = icemodel.netcdf.getdefaults(v, ...
         {'varnames', 'longnames', 'units', 'axes', 'standardnames'});
   end

   % Remove 'depth'
   if strcmp(whichdata, 'ice1')
      drop = ismember(vars.dims, 'depth');
      vars.dims(drop) = [];
      axes.dims(drop) = [];
      units.dims(drop) = [];
      longnames.dims(drop) = [];
      standardnames.dims(drop) = [];
   end

   % Get the spatial dimension data
   dims = icemodel.netcdf.getdimdata(Z, dz);

   % Set the timestep
   switch timeunits
      case 'seconds'
         dt = 3600;
      case 'hours'
         dt = 1;
   end

   % FOR TESTING
   if opts.testwrite
      for f = fieldnames(dims)'
         if ~strcmp(f{:}, {'depth', 'time'})
            dims.(f{:}) = dims.(f{:})(1:opts.numcells);
         end
      end
   end

   % Create the files year by year.
   for m = 1:numyears

      thisyear = num2str(simyears(m));
      filepath = fullfile(pathdata, thisyear);

      % Define the output filename
      filename = strjoin( ...
         {simmodel, whichdata, forcings, userdata, sitename, thisyear, 'nc4'}, '.');

      % Get the input data dimensions
      [ncells, nlyrs, nhrs, chunksize, ...
         vars, units, longnames, standardnames] = getDataDims( ...
         pathdata, thisyear, dims, ...
         vars, units, longnames, standardnames, whichdata);

      % Update the time dimension
      dims.time = 0:dt:dt*(nhrs-1);
      units.dims = strrep(units.dims, 'time', ...
         sprintf('%s since %s-01-01 00:00:00', ...
         timeunits, thisyear));

      % Create the file
      ncid = icemodel.netcdf.create(fullfile(pathsave, filename), ...
         'makebackups', makebackups);
      % job = onCleanup(@() netcdf.close(ncid));

      % Define the dimensions IDs.
      dimid = icemodel.netcdf.defdimid(ncid, ncells, nhrs, nlyrs);

      % Note: dimid is passed to defdimvars b/c each dim has a named dimid i.e.
      % dimid.gridcell, dimid.time, but defdatavars only needs the data dimid so
      % can pass dimid.ice1/2 (or dimid.data if it is redefined in defdimid).

      % Define the grid and time dimensions and attributes.
      icemodel.netcdf.defdimvars(ncid, dimid, vars.dims, ...
         standardnames.dims, longnames.dims, units.dims, axes.dims, ...
         "shuffle", shuffle, "deflate", deflate, "deflateLevel", deflateLevel)

      % Define the data variables.
      icemodel.netcdf.defdatavars(ncid, dimid.(whichdata), vars.(whichdata), ...
         standardnames.(whichdata), longnames.(whichdata), units.(whichdata), ...
         chunksize, xtype, dochunking, "shuffle", shuffle, "deflate", deflate, ...
         "deflateLevel", deflateLevel)

      % Close definition mode
      netcdf.endDef(ncid)

      % Write the grid and time dimensions
      icemodel.netcdf.writedims(ncid, dims, vars.dims);

      % Write the data
      if strcmp(whichdata, 'ice1')

         data = icemodel.netcdf.getvardata(filepath, ...
            ncells, nhrs, xtype, vars.ice1, simmodel);

         icemodel.netcdf.writeice1(ncid, vars.ice1, data);

      elseif strcmp(whichdata, 'ice2')

         icemodel.netcdf.writeice2(ncid, vars.ice2, dims, filepath);
      end

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
function [ncells, nlyrs, nhrs, chunksize, ...
      vars, units, longnames, standardnames] = getDataDims( ...
      pathdata, thisyear, dims, ...
      vars, units, longnames, standardnames, whichdata)

   % Load one file to get the dimensions
   switch whichdata
      case 'ice1'

         ice1 = load(fullfile(pathdata, thisyear, 'ice1_1.mat')).('ice1');

         % Set the dimensions
         ncells = numel(dims.x_easting);        % number of grid cells
         nhrs = height(ice1);                   % number of hours per year
         nlyrs = [];                            % set empty for ice1

         allvars = ice1.Properties.VariableNames;

         chunksize = [ncells, nhrs];            % all cells, annual chunks
         % chunksize = [ncells, 24];            % all cells, daily chunks

      case 'ice2'

         ice2 = load(fullfile(pathdata, thisyear, 'ice2_1.mat')).('ice2');

         % Set the dimensions
         ncells = numel(dims.x_easting);        % number of grid cells
         [nlyrs, ...                            % number of vertical layers
            nhrs] = size(ice2.Tice);            % number of hours per year

         allvars = fieldnames(ice2);

         chunksize = [1, nlyrs, nhrs];             % one cell, all layers, annual
         % chunksize = ceil([1, nlyrs, nhrs/12]);  % one cell, all layers, monthly
         % chunksize = [1, nlyrs, 24];             % one cell, all layers, daily

         % This update the depth dimension values for the case where ice2 files
         % have different Z, but uses the test file loaded above and thus won't
         % catch when it changes. If ice2 files are written to individual nc
         % files then this would need to be placed in writeice2. If multiple
         % ice2 files are written to one nc file, then the supplied opts.dz/Z
         % should be used and possibly NOFILL removed to account for different
         % sized arrays.
         % dims.depth = dz/2:dz:(nlyrs*opts.dz);

         % Overrule nlyrs
         nlyrs = numel(dims.depth);
   end

   % Remove variables from the dictionaries that are not present in ice1/2
   ivars = ismember(vars.(whichdata), allvars);
   vars.(whichdata) = vars.(whichdata)(ivars);
   units.(whichdata) = units.(whichdata)(ivars);
   longnames.(whichdata) = longnames.(whichdata)(ivars);
   standardnames.(whichdata) = standardnames.(whichdata)(ivars);

   % Define chunkSize based on data access patterns. Larger chunk sizes
   % increase memory usage during read/write.
   %
   % When writing, the primary concern is reducing the number of writes.
   % Thus writing the entire array at once is typically ideal, and the
   % netcdf software will determine the chunk size.
   %
   % However, if typical access patterns are known, chunking can improve
   % efficiency by ensuring the data layout in the file matches how it is
   % accessed later. If data is predominantly accessed in large contiguous
   % blocks, having the data and chunks aligned to those patterns is ideal.
end
