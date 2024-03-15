function info = makencfile(pathdata, pathsave, simmodel, forcings, ...
      userdata, sitename, simyears, opts, ncprops)

   arguments
      % required arguments
      pathdata (1, :) char {mustBeFolder}
      pathsave (1, :) char {mustBeFolder}
      simmodel (1, :) char {mustBeTextScalar}
      forcings (1, :) char {mustBeTextScalar}
      userdata (1, :) char {mustBeTextScalar}
      sitename (1, :) char {mustBeTextScalar}
      simyears (:, :) double {mustBeNumeric}

      % optional arguments
      opts.dz  (1, 1) double {mustBeNumeric} = 0
      opts.Z   (1, 1) double {mustBeNumeric} = 0
      opts.whichdata (1, :) char {mustBeMember(opts.whichdata, ...
         {'ice1', 'ice2', 'met'})} = 'ice1'
      opts.timeunits (1, :) char {mustBeMember(opts.timeunits, ...
         {'hours', 'seconds'})} = 'seconds'
      opts.testwrite (1, 1) logical {mustBeNumericOrLogical} = true
      opts.numcells (1, 1) {mustBeNumeric} = 10
      opts.makebackups (1, 1) logical {mustBeNumericOrLogical} = true
      opts.setchunks (1, 1) logical {mustBeNumericOrLogical} = false

      % netcdf api options
      ncprops.format (1, :) char {mustBeMember(ncprops.format, ...
         {'NC_FORMAT_CLASSIC', 'NC_FORMAT_64BIT', ...
         'NC_FORMAT_NETCDF4', 'NC_FORMAT_NETCDF4_CLASSIC'})} ...
         = getenv('ICEMODEL_NC_DEFAULT_FORMAT')

      ncprops.xtype (1, :) char {mustBeMember(ncprops.xtype, ...
         {'NC_FLOAT', 'NC_DOUBLE', 'NC_INT64', 'NC_UINT64', 'NC_INT', ...
         'NC_UINT', 'NC_SHORT', 'NC_USHORT', 'NC_BYTE', 'NC_UBYTE', ...
         'NC_CHAR', 'NC_STRING'})} = 'NC_DOUBLE'
      ncprops.shuffle (1, 1) logical {mustBeNumericOrLogical} = true
      ncprops.deflate (1, 1) logical {mustBeNumericOrLogical} = true
      ncprops.deflateLevel (1, 1) double {mustBeNumeric} = 1
   end

   % NC_INT64, NC_UINT64, NC_UINT, NC_USHORT, NC_UBYTE, NC_STRING only for nc4

   % Update the API configuration
   icemodel.netcdf.config();

   % Set the file format
   oldformat = netcdf.setDefaultFormat(ncprops.format);

   % Pull out the netcdf api options
   [xtype, shuffle, deflate, deflateLevel] = deal( ...
      ncprops.xtype, ncprops.shuffle, ncprops.deflate, ncprops.deflateLevel);

   % Pull out the function options
   [whichdata, Z, dz, timeunits, makebackups, setchunks] = deal( ...
      opts.whichdata, opts.Z, opts.dz, opts.timeunits, opts.makebackups, ...
      opts.setchunks);

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

   % Get the spatial dimension data. Note: this is the only place Z,dz are used.
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
      units.dims{end} = sprintf('%s since %s-01-01 00:00:00', timeunits, thisyear);

      % Create the file
      ncid = icemodel.netcdf.create(fullfile(pathsave, filename), ...
         'makebackups', makebackups);

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
         chunksize, xtype, setchunks, "shuffle", shuffle, "deflate", deflate, ...
         "deflateLevel", deflateLevel)

      % Close definition mode
      netcdf.endDef(ncid)

      % Write the grid and time dimensions
      icemodel.netcdf.writedims(ncid, dims, vars.dims);

      % Get the
      data = icemodel.netcdf.getvardata(filepath, ...
         ncells, nhrs, nlyrs, xtype, vars.(whichdata), simmodel);

      % Write the data
      if strcmp(whichdata, 'ice1')

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

   % Load one file to get a sample of the data
   tmp = load(fullfile(pathdata, thisyear, [whichdata '_1.mat'])).(whichdata);

   switch whichdata
      case 'ice1'
         data = tmp.Tsfc.'; % transpose so it is nlyrs x nhrs
         allvars = tmp.Properties.VariableNames;

      case 'ice2'
         data = tmp.Tice;
         allvars = fieldnames(tmp);
   end

   % Remove variables from the dictionaries that are not present in the data
   [vars, units, longnames, standardnames] = trimvars( ...
      vars, units, longnames, standardnames, allvars, whichdata);

   % Get the dimensions of the data
   ncells = numel(dims.x_easting);  % number of grid cells
   [nlyrs, ...                      % number of hours per year
      nhrs] = size(data);           % number of vertical layers

   % Set the chunksizes
   switch whichdata
      case 'ice1'
         chunksize = [ncells, nhrs];   % all cells, annual chunks

      case 'ice2'

         % This updates the depth dimension values for the case where ice2 files
         % have different Z. However, this uses the test file loaded above, thus
         % if each individual ice2 file is written to an nc file, this would
         % need to be updated in writeice2. If multiple ice2 files are written
         % to one nc file, then the supplied opts.dz/Z should be used and
         % possibly NOFILL removed to account for different sized arrays.

         % dims.depth = dz/2:dz:(nlyrs*opts.dz);

         % Overrule nlyrs with the value set by opts.Z/dz
         nlyrs = numel(dims.depth);

         chunksize = [1, nlyrs, nhrs]; % one cell, all layers, annual
   end

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

%% set time units
function units = settimeunits(units, timeunits, thisyear)
   itime = contains(units.dims, '00:00:00');
   units.dims{itime} = sprintf('%s since %s-01-01 00:00:00', timeunits, thisyear);
end

   ivars = ismember(vars.(whichdata), allvars);
   vars.(whichdata) = vars.(whichdata)(ivars);
   units.(whichdata) = units.(whichdata)(ivars);
   longnames.(whichdata) = longnames.(whichdata)(ivars);
   standardnames.(whichdata) = standardnames.(whichdata)(ivars);
end
