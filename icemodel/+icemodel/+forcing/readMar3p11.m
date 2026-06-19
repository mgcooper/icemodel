function [data, units, Time] = readMar3p11(filename, varname, kwargs)
   %READMAR3P11 Read one MAR v3.11 variable in standard units.
   %
   %  [data, units, Time] = icemodel.forcing.readMar3p11(filename, varname)
   %  [data, units, Time] = ... readMar3p11(_, start=[i j], count=[ni nj])
   %
   % Reads a MAR v3.11 NetCDF variable (optionally a spatial hyperslab)
   % and converts the legacy MAR units to icemodel-standard ones:
   %
   %    C      -> K        (air/surface temperature)
   %    g/kg   -> kg/kg    (specific humidity)
   %    mmWE/h -> mWE/h    (precipitation / melt / runoff / SMB)
   %    hPa    -> Pa       (surface pressure)
   %
   % MAR yearly files store hourly variables as [nx ny 24 ndays] and
   % daily variables as [nx ny ndays] (or [nx ny 1 ndays]); both collapse
   % to a cells-by-time list with cells flattened in native grid order
   % (matching marGridInfo linear indices) and time across columns.
   %
   % Inputs
   %  filename - MAR yearly NetCDF (e.g. MARv3.11-ERA5-15km-2009.nc)
   %  varname  - MAR variable name (TTH, SWDH, SHSN2, ...)
   %
   % Name-value
   %  start, count - optional grid hyperslab: start cell [i j] (1-based)
   %                 and extent [ni nj]. Default reads the full grid.
   %
   % Outputs
   %  data  - (ncells x ntime) double in standard units
   %  units - unit string after conversion
   %  Time  - UTC datetime axis (hourly or daily, matching the variable).
   %          MAR yearly files begin Jan 1 00:00 UTC of the file year.
   %
   % See also: icemodel.forcing.marGridInfo, icemodel.forcing.buildMarMet

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

   % Collapse to cells x time. ncread returns [ni nj (24) ndays]; squeeze
   % drops singleton dims, so reshape from the requested cell count.
   ncells = prod(count(1:2));
   data = reshape(data, ncells, []);

   % No-data handling (before unit conversion, so sentinels are caught at
   % their native magnitude). MAR uses two no-data conventions: a ~1e36
   % _FillValue (unambiguous for any variable) and a 999 sentinel that
   % appears in the mass-flux fields (notably RUH/SMBH outside the melt
   % season; see runoff/NEW/test/BAD_MAR_DATA.m). A blanket 999 threshold
   % is unsafe for e.g. SWDH (~1000 W m-2), so the 999 cut is applied only
   % to the mmWE/h mass fluxes, where a legitimate value never reaches it.
   % Left as NaN here; downstream metchecks / conservative-remap inpainting
   % handle the gaps. Without this filter a single 1e36 cell would corrupt
   % any catchment mean/sum that includes it.
   data(data >= 1e30) = NaN;
   if strcmp(units, 'mmWE/h')
      data(data >= 999) = NaN;
   end

   % Standard unit conversions (legacy readMar3p11 table).
   switch units
      case 'C'
         data = data + 273.15;
         units = 'K';
      case 'g/kg'
         data = data / 1000;
         units = 'kg/kg';
      case 'mmWE/h'
         data = data / 1000;
         units = 'mWE/h';
      case 'hPa'
         data = data * 100;
         units = 'Pa';
   end

   % Time axis: yearly files start Jan 1 00:00 UTC; hourly variables have
   % 24 samples per day, daily variables one.
   yyyy = marFileYear(filename);
   ndays = dims(end);
   ntime = size(data, 2);
   t0 = datetime(yyyy, 1, 1, 'TimeZone', 'UTC');
   if ntime == 24 * ndays
      Time = (t0:hours(1):(t0 + days(ndays) - hours(1)))';
   else
      Time = (t0:days(1):(t0 + days(ndays - 1)))';
   end
end

function yyyy = marFileYear(filename)
   %MARFILEYEAR Parse the calendar year from a MAR yearly filename.
   tok = regexp(filename, '(\d{4})\.nc$', 'tokens', 'once');
   assert(~isempty(tok), ...
      'cannot parse the file year from %s (expected ...-YYYY.nc)', filename)
   yyyy = str2double(tok{1});
end
