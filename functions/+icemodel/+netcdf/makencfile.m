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

      opts.GetSizeFromData (1, :) logical {mustBeNumericOrLogical} = true
      opts.GetSizeFromDims (1, :) logical {mustBeNumericOrLogical} = false

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

   % ---------------------------------------------------------------------
   % ----------------------------  Parse inputs and set parameters

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

   % Set the timestep
   T = 8760; % update inside the loop to support leap years
   switch timeunits
      case 'seconds'
         dt = 3600;
      case 'hours'
         dt = 1;
   end
   numyears = numel(simyears);

   % ---------------------------------------------------------------------
   % ----------------------------  Define variables and attributes.

   for v = ["dims", whichdata]
      [varnames.(v), longnames.(v), units.(v), axes.(v), standardnames.(v)] ...
         = icemodel.netcdf.getdefaults(v, ...
         {'varnames', 'longnames', 'units', 'axes', 'standardnames'});
   end

   % I commented this out to test if I could get ice1 to work w/o it and it
   % does, the reason it works is because of the try-catch in writedims and the
   % if-else in defdimvars. Once getdefaults is refactored to use explcit
   % fieldnames, those should be able to be removed.
   %
   % % Remove 'depth' dim from ice1 - this is done for defdimvars
   % [varnames, axes, units, longnames, standardnames] = dropdims(whichdata, ...
   %    varnames, axes, units, longnames, standardnames);

   % ---------------------------------------------------------------------
   % ----------------------------  Create the files year by year
   for m = 1:numyears

      thisyear = num2str(simyears(m));
      filepath = fullfile(pathdata, thisyear);

      % Define the output filename
      filename = fullfile(pathsave, strjoin( ...
         {simmodel, whichdata, forcings, userdata, sitename, thisyear, 'nc4'}, '.'));

      % Get the spatial and time dimension data. This is the only place Z,dz,T,dt are used.
      dimdata = icemodel.netcdf.getdimdata(Z, dz, T, dt);

      % Need to update dimdata directly from data here before passing dimdata to
      % other functions which use getdimsize.
      [datavars, datasize] = icemodel.netcdf.getvarinfo(filepath, whichdata, 1);

      % NOTE: I think I need to actually reset dimdata with the datasize
      if opts.GetSizeFromData

         % The problem is when Z, dz are default 0, then dimdata.depth is an
         % empty double, so it needs to be getdimsfromdata ... so this is just a
         % test here to reset dimdata.depth but either need to enforce explicit
         % Z, dz or refactor

         % dimdata.depth = ... yeah so this is the problem, we don't know Z, dz
         % so just enforce it

      elseif opts.GetSizeFromDims
      end



      % ------------------------------------------------------------------
      % Content in this section is to be removed or only for testing

      % For testing on fewer grid cells
      dimdata = resetDimData(dimdata, opts);

      % Remove variables from the dictionaries that are not present in the data

      % Note: The problem here is if the datavars change, then varnames, units,
      % longnames, and standardnames are trimmed and cannot be reset

      [varnames, units, longnames, standardnames] = trimvars(whichdata, ...
         varnames, units, longnames, standardnames, datavars);

      % If setting the dimensions directly from the data, then getchunksize and
      % defdimid need to get the same values.

      % ------------------------------------------------------------------

      chunksizes = icemodel.netcdf.getchunksize(whichdata, dimdata, datasize, ...
         "GetSizeFromData", true, "GetSizeFromDims", false);

      % Update the time units for this year
      units = settimeunits(units, timeunits, thisyear);

      % Create the file
      ncid = icemodel.netcdf.create(filename, 'makebackups', makebackups);

      % Define the dimensions IDs.
      dimid = icemodel.netcdf.defdimid(ncid, dimdata, datasize, ...
         "GetSizeFromData", true, "GetSizeFromDims", false);

      % Define the grid and time dimensions and attributes.
      icemodel.netcdf.defdimvars(ncid, dimid, varnames.dims, ...
         standardnames.dims, longnames.dims, units.dims, axes.dims, ...
         "shuffle", shuffle, "deflate", deflate, "deflateLevel", deflateLevel)

      % Define the data variables.
      icemodel.netcdf.defdatavars(ncid, dimid.(whichdata), varnames.(whichdata), ...
         standardnames.(whichdata), longnames.(whichdata), units.(whichdata), ...
         chunksizes, xtype, setchunks, "shuffle", shuffle, "deflate", deflate, ...
         "deflateLevel", deflateLevel)

      % Close definition mode
      netcdf.endDef(ncid)

      % Write the grid and time dimensions
      icemodel.netcdf.writedims(ncid, dimdata, varnames.dims);

      % Get the
      data = icemodel.netcdf.getvardata(filepath, varnames.(whichdata), ...
         dimdata, xtype, simmodel);

      % Write the data
      if strcmp(whichdata, 'ice1')

         icemodel.netcdf.writeice1(ncid, varnames.ice1, data);

      elseif strcmp(whichdata, 'ice2')

         icemodel.netcdf.writeice2(ncid, varnames.ice2, dimdata, filepath);
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

%% trim varnames
function [varnames, units, longnames, standardnames] = trimvars(whichdata, ...
      varnames, units, longnames, standardnames, datavars)

   ivars = ismember(varnames.(whichdata), datavars);
   varnames.(whichdata) = varnames.(whichdata)(ivars);
   units.(whichdata) = units.(whichdata)(ivars);
   longnames.(whichdata) = longnames.(whichdata)(ivars);
   standardnames.(whichdata) = standardnames.(whichdata)(ivars);
end

%% drop dims
function [varnames, axes, units, longnames, standardnames] = dropdims( ...
      whichdata, varnames, axes, units, longnames, standardnames)

   if strcmp(whichdata, 'ice1')
      drop = ismember(varnames.dims, 'depth');
      varnames.dims(drop) = [];
      axes.dims(drop) = [];
      units.dims(drop) = [];
      longnames.dims(drop) = [];
      standardnames.dims(drop) = [];
   end
end

%% set time units
function units = settimeunits(units, timeunits, thisyear)
   itime = contains(units.dims, '00:00:00');
   units.dims{itime} = sprintf('%s since %s-01-01 00:00:00', timeunits, thisyear);
end

%% reset testdims
function dims = trimGridCells(dims, testwrite, numcells)
   % FOR TESTING
   if testwrite
      for f = fieldnames(dims)'
         if ~strcmp(f{:}, {'depth', 'time'})
            dims.(f{:}) = dims.(f{:})(1:numcells);
         end
      end
   end
end
