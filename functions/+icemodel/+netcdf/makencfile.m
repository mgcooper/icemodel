function info = makencfile(datafile, datapath, savepath, simmodel, forcings, ...
      userdata, sitename, simyears, opts, ncprops)

   arguments
      % required arguments
      datafile (1, :) char {mustBeMember(datafile, ...
         {'ice1', 'ice2', 'met'})}
      datapath (1, :) char {mustBeFolder}
      savepath (1, :) char {mustBeFolder}
      simmodel (1, :) char {mustBeTextScalar}
      forcings (1, :) char {mustBeTextScalar}
      userdata (1, :) char {mustBeTextScalar}
      sitename (1, :) char {mustBeTextScalar}
      simyears (:, :) double {mustBeNumeric}

      % optional arguments
      opts.Z   (1, 1) double {mustBeNumeric} = 0
      opts.dz  (1, 1) double {mustBeNumeric} = 0
      opts.T   (1, 1) double {mustBeNumeric} = 8760
      opts.dt  (1, 1) double {mustBeNumeric} = 3600
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
         = 'NC_FORMAT_NETCDF4' % getenv('ICEMODEL_NC_DEFAULT_FORMAT')

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

   % Validate the grid spacing.
   if datafile == "ice2" && opts.Z == 0
      error('Set Z and dz if saving ice2 data')
   end

   % Validate the timestep.
   timeunits = opts.timeunits;
   switch timeunits
      case 'seconds'
         assert(opts.dt == 3600)
      case 'hours'
         assert(opts.dt == 1)
   end

   % Set the file format.
   oldformat = netcdf.setDefaultFormat(ncprops.format);
   cleanupfn = onCleanup(@() netcdf.setDefaultFormat(oldformat));

   % ---------------------------------------------------------------------
   % ----------------------------  Create the files year by year

   for thisyear = string(simyears)

      % Define the output filename.
      filename = fullfile(savepath, strjoin( ...
         {simmodel, datafile, forcings, userdata, sitename, char(thisyear), ...
         'nc4'}, '.'));

      % Update the time units using thisyear as the datum.
      opts.timeunits = gettimeunits(timeunits, thisyear);

      % Create the file.
      processOneYear(fullfile(datapath, thisyear), ...
         datafile, filename, ncprops, opts, simmodel);
   end

   % Parse outputs.
   if nargout > 1
      info = ncinfo(filename);
   end
end
%%

% Note: The call to icemodel.netcdf.getdefaults returns varnames, units, axes,
% etc, and then trimvars removes variables which are not present in the actual
% ice1/2 data files. Rather than call getdefaults once in the main function, it
% is included in this subfunction to account for the case where the data files
% change year by year, otherwise trimvars could remove variables in one year
% which are present in a later year.
%
% A similar issue occurs if the depth dimension (or time) changes year by year.
% The changing depth dimension is nominally handled by the call to getvarinfo,
% and the "GetSizeFromData" vs "GetSizeFromDims" options.
%
% However, neither trimvars nor the "GetSize" options account for the case where
% the vars or dims change within a year from file to file.
%
% For dims, specifically depth, the important thing is setting the dimsizes in
% the files to the maximum depth so if some files have 300 layers and other 500,
% the nc files are defined to have 500 layers and when a file with 300 layers is
% encountered the data is written to the first 300 layers. Thus Z and dz are
% used to set the file-wise dims, which means using GetSizeFromDims == true, and
% GetSizeFromData could be removed altogether. In general they should be
% interchangeable and the Z,dz inputs could be removed to simplify the
% interface.
%
% There may be a use case for the two separate "GetSize" paths - Z, dz could
% be used to set the file-wise dims, where getdimsize returns the size of the
% depth grid defined by Z, dz. But GetSizeFromData is used in getchunksize so
% the chunks match the actual data ... but actually that's not right either, the
% chunksize is file-wise. So there may not be any use case for GetSizeFromData
% unless we want to eliminated Z, dz and rely entirely on the data.

% The use case could just be for validating consistent datasize and dimsize

function processOneYear(datapath, datafile, filename, ncprops, opts, simmodel)

   % Note: simmodel is only added as an input to patch the skinmodel
   % ice1.freeze data. Once those files are written, remove simmodel.

   % Pull out the netcdf api options
   [xtype, shuffle, deflate, deflateLevel] = deal( ...
      ncprops.xtype, ncprops.shuffle, ncprops.deflate, ncprops.deflateLevel);

   % Pull out the function options
   [Z, dz, T, dt, timeunits, makebackups, setchunks] = deal( ...
      opts.Z, opts.dz, opts.T, opts.dt, opts.timeunits, ...
      opts.makebackups, opts.setchunks);

   % ---------------------------------------------------------------------
   % ----------------------------  Define variables and attributes.
   [varnames, longnames, units, axes, standardnames] ...
      = icemodel.netcdf.getdefaults(["dims", datafile], ...
      {'varnames', 'longnames', 'units', 'axes', 'standardnames'});

   % Get the spatial and time dimension data.
   dimdata = icemodel.netcdf.getdimdata(Z, dz, T, dt, whichmask="icemask");

   % Reduce the number of grid cells if testwrite == true
   dimdata = trimGridCells(dimdata, opts.testwrite, opts.numcells);

   % Update dimdata from the data before passing to getdimsize.
   [datavars, datasize] = icemodel.netcdf.getvarinfo(datapath, datafile, 1);

   % Remove varnames (and associated units etc) that are not present in datavars
   [varnames, units, longnames, standardnames] = trimvars(datafile, ...
      varnames, units, longnames, standardnames, datavars);

   % Define the chunksizes. If setting the dimensions directly from the data,
   % then getchunksize and defdimid need to get the same values.
   chunksizes = icemodel.netcdf.getchunksize(datafile, dimdata, datasize, ...
      "GetSizeFromData", true, "GetSizeFromDims", false);

   % Update the time units for this year
   units = settimeunits(units, timeunits);

   % Create the file
   ncid = icemodel.netcdf.create(filename, 'makebackups', makebackups);
   job = onCleanup(@() netcdf.close(ncid));

   % Define the dimensions IDs.
   dimid = icemodel.netcdf.defdimid(ncid, dimdata, datasize);

   % Define the grid and time dimensions and attributes.
   icemodel.netcdf.defdimvars(ncid, dimid, varnames.dims, ...
      standardnames.dims, longnames.dims, units.dims, axes.dims, ...
      "shuffle", shuffle, "deflate", deflate, "deflateLevel", deflateLevel)

   % Define the data variables.
   icemodel.netcdf.defdatavars(ncid, dimid.(datafile), varnames.(datafile), ...
      standardnames.(datafile), longnames.(datafile), units.(datafile), ...
      chunksizes, xtype, setchunks, "shuffle", shuffle, "deflate", deflate, ...
      "deflateLevel", deflateLevel)

   % Close definition mode
   netcdf.endDef(ncid)

   % Write the grid and time dimensions
   icemodel.netcdf.writedims(ncid, dimdata, varnames.dims);

   % Write the data
   if strcmp(datafile, 'ice1')

      icemodel.netcdf.writeice1(ncid, datapath, varnames.(datafile), ...
         dimdata, xtype, simmodel);

   elseif strcmp(datafile, 'ice2')

      icemodel.netcdf.writeice2(ncid, datapath, varnames.(datafile), ...
         dimdata, xtype, simmodel);
   end
end

%% trim varnames
function [varnames, units, longnames, standardnames] = trimvars(datafile, ...
      varnames, units, longnames, standardnames, datavars)

   ivars = ismember(varnames.(datafile), datavars);
   varnames.(datafile) = varnames.(datafile)(ivars);
   units.(datafile) = units.(datafile)(ivars);
   longnames.(datafile) = longnames.(datafile)(ivars);
   standardnames.(datafile) = standardnames.(datafile)(ivars);
end

%% set time units
function timeunits = gettimeunits(baseunits, thisyear)
   timeunits = sprintf('%s since %s-01-01 00:00:00', baseunits, thisyear);
end
function units = settimeunits(units, timeunits)
   itime = contains(units.dims, '00:00:00');
   units.dims{itime} = timeunits;
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
