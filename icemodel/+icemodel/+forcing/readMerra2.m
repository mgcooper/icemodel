function [data, units, Time] = readMerra2(filename, varname, kwargs)
   %READMERRA2 Read one MERRA-2 variable in standard units.
   %
   %  [data, units, Time] = icemodel.forcing.readMerra2(filename, varname)
   %  [data, units, Time] = ... readMerra2(_, start=[i j], count=[ni nj])
   %
   % Reads one variable from a single MERRA-2 daily collection file
   % (tavg1_2d_{slv,rad,flx}_Nx or tavg3_2d_glc_Nx; optionally a spatial
   % hyperslab) and converts the native units to icemodel-standard ones.
   % Mirrors icemodel.forcing.readMar3p11 / readRacmo2p3 and the legacy
   % merra.readMerra2, so the gridded-source readers share a contract: a
   % cells-by-time block (cells flattened in native [lon lat] order, matching
   % the X, Y grids and gridLocation's column-major slab), the unit string,
   % and the UTC time axis.
   %
   % MERRA-2 surface collections are dimensioned [lon lat time], with 24
   % hourly bins (tavg1) or 8 three-hourly bins (tavg3/glc) per daily file.
   %
   % Unit conversions (native -> standard):
   %    kg m-2 s-1 -> mWE/h   (mass-flux RATE: PRECTOTCORR, PRECSNO, EVAP, RUNOFF)
   %    C          -> K
   %    g/kg       -> kg/kg
   %    hPa        -> Pa
   %    W m-2      -> W/m2    (label only)
   %    m s-1      -> m/s     (label only)
   % kg m-2 (a STORE, e.g. SNOMAS_GL/swe) is left untouched - only the
   % per-second flux RATE is scaled to mWE/h.
   %
   % Inputs
   %  filename - MERRA-2 daily NetCDF (e.g. MERRA2_400.tavg1_2d_slv_Nx....nc4)
   %  varname  - NetCDF variable name (e.g. 'T2M', 'PRECTOTCORR', 'RUNOFF')
   %
   % Name-value
   %  start, count - optional grid hyperslab: start cell [i j] (1-based) and
   %                 extent [ni nj] over [lon lat]. Default reads the full grid.
   %
   % Outputs
   %  data  - (ncells x ntime) double in standard units, native grid order
   %  units - unit string after conversion
   %  Time  - UTC datetime axis from the file (computed only when requested)
   %
   % See also: icemodel.forcing.readMar3p11, icemodel.forcing.readRacmo2p3,
   %  icemodel.forcing.buildMerraData

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

   % Assemble the read window: requested spatial hyperslab, all times.
   if isempty(kwargs.start)
      start = ones(1, numel(dims));
      count = dims;
   else
      start = [kwargs.start, ones(1, numel(dims) - 2)];
      count = [kwargs.count, dims(3:end)];
   end
   data = double(squeeze(ncread(filename, varname, start, count)));

   % Collapse to cells x time (cells flattened in native [lon lat] order).
   ncells = prod(count(1:2));
   data = reshape(data, ncells, []);

   % MERRA-2 _FillValue (~1e15) over masked tiles (e.g. glacier-tile
   % variables off-ice) -> NaN, so it never corrupts a downstream mean/sum.
   data(data >= 1e14) = NaN;

   % Standard unit conversions (shared reader family). Note kg m-2 s-1 (a
   % flux rate) converts to mWE/h, but kg m-2 (a store, e.g. swe) does not.
   switch units
      case 'kg m-2 s-1'
         data = data * 3600 / 1000;   % -> meters water equivalent per hour
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
      case 'm s-1'
         units = 'm/s';
   end

   if nargout >= 3
      Time = merraTime(filename);
   end
end

%% Local functions
function Time = merraTime(filename)
   %MERRATIME UTC bin-center datetime axis from the MERRA-2 'time' variable.
   % MERRA-2 stores time as "minutes since <yyyy-mm-dd hh:mm:ss>".
   t = double(ncread(filename, 'time'));
   t_units = ncreadatt(filename, 'time', 'units');
   tok = regexp(t_units, 'minutes since (\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2})', ...
      'tokens', 'once');
   assert(~isempty(tok), 'unexpected MERRA time units: %s', t_units)
   t0 = datetime(tok{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
   Time = t0 + minutes(t);
end
