function [data, units, Time] = readRacmo2p3(filename, varname, kwargs)
   %READRACMO2P3 Read one RACMO 2.3 (FGRN11) variable in standard units.
   %
   %  [data, units, Time] = icemodel.forcing.readRacmo2p3(filename, varname)
   %  [data, units, Time] = ... readRacmo2p3(_, start=[i j], count=[ni nj])
   %  [blocks, units, Time] = ... readRacmo2p3(_, slabs={[i j;ni nj], ...})
   %
   % Reads one variable from a per-variable RACMO 2.3p3 FGRN11 NetCDF file
   % (optionally a spatial hyperslab) and converts the native units to
   % icemodel-standard ones. Mirrors icemodel.forcing.readMar3p11 / the
   % legacy runoff readRacmo2p3, so the gridded-source readers share a
   % contract: a cells-by-time block (cells flattened in native grid order,
   % matching the X, Y grids from the same file and gridLocation's column-
   % major slab), the unit string, and the UTC time axis.
   %
   % RACMO per-variable files are dimensioned [rlon rlat height(=1) time]
   % on the rotated-pole FGRN11 grid, posted 3-hourly. The singleton height
   % level is squeezed out.
   %
   % Unit conversions (native -> standard):
   %    kg m-2 s-1 -> mWE/h   (mass fluxes: runoff, melt, precip, ...)
   %    mmWE/h     -> mWE/h
   %    C          -> K
   %    g/kg       -> kg/kg
   %    hPa        -> Pa
   %    W m-2      -> W/m2    (label only)
   %
   % NOTE on mass fluxes: kg m-2 s-1 is converted to meters water equivalent
   % per HOUR (x3600/1000), i.e. it represents the 3-hourly-mean rate. To
   % accumulate, either keep 3-hourly posting and multiply by 3 before
   % cumsum, or interpolate to hourly first (buildRacmoData does the latter).
   %
   % Inputs
   %  filename - RACMO per-variable NetCDF (e.g. runoff.RACMO23p3_..._3H.nc)
   %  varname  - NetCDF variable name (e.g. 'runoff', 'latf', 'swsd')
   %
   % Name-value
   %  start, count - optional grid hyperslab: start cell [i j] (1-based) and
   %                 extent [ni nj] over [rlon rlat]. Default reads the full grid.
   %  slabs        - optional cell array of [start; count] 2x2 hyperslab specs
   %                 ({[i j; ni nj], ...}). When given, the (multi-GB) file is
   %                 OPENED ONCE and every listed hyperslab is read from the same
   %                 open file, returning a cell array of blocks (one per slab)
   %                 instead of a single matrix - the batch path that extracts
   %                 many points from one variable file without re-opening it per
   %                 point. start/count are ignored when slabs is given.
   %
   % Outputs
   %  data  - (ncells x ntime) double in standard units, native grid order;
   %          OR, when slabs is given, a cell array {(ncells x ntime), ...}
   %  units - unit string after conversion
   %  Time  - UTC datetime axis (computed only when requested)
   %
   % See also: icemodel.forcing.readMar3p11, icemodel.forcing.readMerra2,
   %  icemodel.forcing.buildRacmoData

   arguments
      filename (1, 1) string
      varname (1, 1) string
      kwargs.start (1, :) double = []
      kwargs.count (1, :) double = []
      kwargs.slabs cell = {}
   end

   [dims, native_units] = ncVarInfo(filename, varname);

   if isempty(kwargs.slabs)
      [start, count] = slabWindow(dims, kwargs.start, kwargs.count);
      data = convertSlab(double(squeeze( ...
         ncread(filename, varname, start, count))), native_units, count);
   else
      % Batch path: open once, read every listed hyperslab from the same ncid.
      % RACMO files are the multi-GB subsurface product, so re-opening per
      % point dominates a staging job - this opens the file exactly once.
      ncid = netcdf.open(filename, 'NOWRITE');
      cleanup = onCleanup(@() netcdf.close(ncid));
      varid = netcdf.inqVarID(ncid, char(varname));
      data = cell(size(kwargs.slabs));
      for k = 1:numel(kwargs.slabs)
         [start, count] = slabWindow(dims, ...
            kwargs.slabs{k}(1, :), kwargs.slabs{k}(2, :));
         data{k} = convertSlab( ...
            double(squeeze(netcdf.getVar(ncid, varid, start - 1, count))), ...
            native_units, count);
      end
   end

   units = convertUnits(native_units);

   if nargout >= 3
      Time = racmoTime(filename);
   end
end

%% Local functions
function [dims, units] = ncVarInfo(filename, varname)
   %NCVARINFO Variable dimensions and unit string from ncinfo.
   info = ncinfo(filename, varname);
   dims = info.Size;
   units = '';
   has_units = strcmp({info.Attributes.Name}, 'units');
   if any(has_units)
      units = info.Attributes(has_units).Value;
   end
end

function [start, count] = slabWindow(dims, kstart, kcount)
   %SLABWINDOW Read window: requested spatial hyperslab, singleton level, all
   % times. RACMO variables are [rlon rlat height(=1) time].
   if isempty(kstart)
      start = ones(1, numel(dims));
      count = dims;
   else
      start = [kstart, ones(1, numel(dims) - 2)];
      count = [kcount, dims(3:end)];
   end
end

function data = convertSlab(data, units, count)
   %CONVERTSLAB Reshape to cells x time and convert to standard units.
   % Shared by the ncread (single) and netcdf.getVar (batch) read paths.

   % Collapse to cells x time (cells flattened in native [rlon rlat] order).
   ncells = prod(count(1:2));
   data = reshape(data, ncells, []);

   % Standard unit conversions (shared with readMar3p11; the RACMO archive
   % posts mass fluxes as kg m-2 s-1 rather than MAR's mmWE/h). The unit
   % STRING is relabelled once at the top level by convertUnits.
   switch units
      case 'kg m-2 s-1'
         data = data * 3600 / 1000;   % -> meters water equivalent per hour
      case 'mmWE/h'
         data = data / 1000;
      case 'C'
         data = data + 273.15;
      case 'g/kg'
         data = data / 1000;
      case 'hPa'
         data = data * 100;
   end
end

function units = convertUnits(units)
   %CONVERTUNITS Relabel a native RACMO unit string to its standard form.
   % Mirrors the value scaling in convertSlab (data and label stay in sync).
   switch units
      case 'kg m-2 s-1'
         units = 'mWE/h';
      case 'mmWE/h'
         units = 'mWE/h';
      case 'C'
         units = 'K';
      case 'g/kg'
         units = 'kg/kg';
      case 'hPa'
         units = 'Pa';
      case 'W m-2'
         units = 'W/m2';
   end
end

function Time = racmoTime(filename)
   %RACMOTIME UTC datetime axis from the RACMO 'time' variable.
   t = double(ncread(filename, 'time'));
   t_units = ncreadatt(filename, 'time', 'units');
   assert(startsWith(t_units, 'days since 1950-01-01'), ...
      'unexpected RACMO time units: %s', t_units)
   Time = datetime(1950, 1, 1, 'TimeZone', 'UTC') + days(t);
end
