function [data, units, Time] = readRacmo2p3(filename, varname, kwargs)
   %READRACMO2P3 Read one RACMO 2.3 (FGRN11) variable in standard units.
   %
   %  [data, units, Time] = icemodel.forcing.readRacmo2p3(filename, varname)
   %  [data, units, Time] = ... readRacmo2p3(_, start=[i j], count=[ni nj])
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
   %
   % Outputs
   %  data  - (ncells x ntime) double in standard units, native grid order
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
   end

   info = ncinfo(filename, varname);
   dims = info.Size;
   units = '';
   has_units = strcmp({info.Attributes.Name}, 'units');
   if any(has_units)
      units = info.Attributes(has_units).Value;
   end

   % Assemble the read window: requested spatial hyperslab, singleton level,
   % all times. RACMO variables are [rlon rlat height(=1) time].
   if isempty(kwargs.start)
      start = ones(1, numel(dims));
      count = dims;
   else
      start = [kwargs.start, ones(1, numel(dims) - 2)];
      count = [kwargs.count, dims(3:end)];
   end
   data = double(squeeze(ncread(filename, varname, start, count)));

   % Collapse to cells x time (cells flattened in native [rlon rlat] order).
   ncells = prod(count(1:2));
   data = reshape(data, ncells, []);

   % Standard unit conversions (shared with readMar3p11; the RACMO archive
   % posts mass fluxes as kg m-2 s-1 rather than MAR's mmWE/h).
   switch units
      case 'kg m-2 s-1'
         data = data * 3600 / 1000;   % -> meters water equivalent per hour
         units = 'mWE/h';
      case 'mmWE/h'
         data = data / 1000;
         units = 'mWE/h';
      case 'C'
         data = data + 273.15;
         units = 'K';
      case 'g/kg'
         data = data / 1000;
         units = 'kg/kg';
      case 'hPa'
         data = data * 100;
         units = 'Pa';
      case 'W m-2'
         units = 'W/m2';
   end

   if nargout >= 3
      Time = racmoTime(filename);
   end
end

%% Local functions
function Time = racmoTime(filename)
   %RACMOTIME UTC datetime axis from the RACMO 'time' variable.
   t = double(ncread(filename, 'time'));
   t_units = ncreadatt(filename, 'time', 'units');
   assert(startsWith(t_units, 'days since 1950-01-01'), ...
      'unexpected RACMO time units: %s', t_units)
   Time = datetime(1950, 1, 1, 'TimeZone', 'UTC') + days(t);
end
