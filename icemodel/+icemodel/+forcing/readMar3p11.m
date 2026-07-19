function [data, units, Time] = readMar3p11(filename, varname, kwargs)
   %READMAR3P11 Read one MAR v3.11 variable in standard units.
   %
   %  [data, units, Time] = icemodel.forcing.readMar3p11(filename, varname)
   %  [data, units, Time] = ... readMar3p11(_, start=[i j], count=[ni nj])
   %  [blocks, units, Time] = ... readMar3p11(_, slabs={[i j;ni nj], ...})
   %  [blocks, units, Time] = ... readMar3p11(_, slabs=..., sector=[1 2])
   %
   % Reads a MAR v3.11 NetCDF variable (optionally a spatial hyperslab)
   % and converts the legacy MAR units to icemodel-standard ones:
   %
   %    C      -> K        (air/surface temperature)
   %    g/kg   -> kg/kg    (specific humidity)
   %    mmWE/h -> mWE/h    (hourly precipitation / melt / runoff / SMB /
   %                         sublimation)
   %    mmWE/day -> mWE/day (daily mass-balance diagnostics)
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
   %  slabs        - optional cell array of [start; count] 2x2 hyperslab
   %                 specs ({[i j; ni nj], ...}). When given, the file is
   %                 OPENED ONCE and every listed hyperslab is read from
   %                 the same open file, returning a cell array of blocks
   %                 (one per slab) instead of a single matrix. This is the
   %                 batch path used to extract many points from one source
   %                 file without re-opening it per point. start/count are
   %                 ignored when slabs is given.
   %  sector       - optional surface-sector index for variables whose third
   %                 dimension is SECTOR (1 = permanent ice, 2 = tundra).
   %                 Supply one value for a single slab or one value per slab.
   %
   % Outputs
   %  data  - (ncells x ntime) double in standard units; OR, when slabs is
   %          given, a cell array {(ncells x ntime), ...} aligned to slabs
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
      kwargs.slabs cell = {}
      kwargs.sector (1, :) double {mustBeInteger, mustBePositive} = []
   end

   [dims, native_units] = ncVarInfo(filename, varname);

   if isempty(kwargs.slabs)
      % Single-hyperslab path (the original contract, returns one matrix).
      assert(numel(kwargs.sector) <= 1, ...
         'a single MAR slab accepts at most one surface sector')
      [start, count] = slabWindow( ...
         dims, kwargs.start, kwargs.count, kwargs.sector);
      data = readSlab(filename, varname, start, count, native_units);
   else
      % Batch path: open the file ONCE and read each requested hyperslab
      % from the same open file (one ncid for all listed points), so a
      % staging job does not re-open the file per point.
      ncid = netcdf.open(filename, 'NOWRITE');
      cleanup = onCleanup(@() netcdf.close(ncid));
      varid = netcdf.inqVarID(ncid, char(varname));
      data = cell(size(kwargs.slabs));
      sectors = normalizeSectors(kwargs.sector, numel(kwargs.slabs));
      for k = 1:numel(kwargs.slabs)
         [start, count] = slabWindow(dims, ...
            kwargs.slabs{k}(1, :), kwargs.slabs{k}(2, :), sectors(k));
         data{k} = convertSlab( ...
            ncGetSlab(ncid, varid, start, count), native_units, count);
      end
   end

   units = convertUnits(native_units);

   if nargout >= 3
      Time = marTime(filename, dims);
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

function [start, count] = slabWindow(dims, kstart, kcount, sector)
   %SLABWINDOW Assemble the read window: requested spatial hyperslab, all times.
   if isempty(kstart)
      start = ones(1, numel(dims));
      count = dims;
   else
      start = [kstart, ones(1, numel(dims) - 2)];
      count = [kcount, dims(3:end)];
   end

   % Sector selection keeps the native daily time dimension while removing the
   % permanent-ice/tundra dimension. It is valid only for 4-D daily products.
   if ~isempty(sector) && isfinite(sector)
      assert(numel(dims) >= 4 && sector <= dims(3), ...
         'requested MAR surface sector is not present in this variable')
      start(3) = sector;
      count(3) = 1;
   end
end

function sectors = normalizeSectors(sector, nslab)
   %NORMALIZESECTORS Broadcast or validate one sector per requested slab.
   if isempty(sector)
      sectors = nan(1, nslab);
   elseif isscalar(sector)
      sectors = repmat(sector, 1, nslab);
   else
      assert(numel(sector) == nslab, ...
         'MAR sector must be scalar or have one value per slab')
      sectors = reshape(sector, 1, []);
   end
end

function data = readSlab(filename, varname, start, count, units)
   %READSLAB ncread one hyperslab and standardize (single-open ncread path).
   data = double(squeeze(ncread(filename, varname, start, count)));
   data = convertSlab(data, units, count);
end

function data = ncGetSlab(ncid, varid, start, count)
   %NCGETSLAB Read one hyperslab from an already-open file (0-based netcdf API).
   % netcdf.getVar uses 0-based start; ncread uses 1-based. The squeeze/double
   % matches the ncread path so the two readers return identical arrays.
   data = double(squeeze(netcdf.getVar(ncid, varid, start - 1, count)));
end

function data = convertSlab(data, units, count)
   %CONVERTSLAB Reshape to cells x time, mask no-data, convert to std units.
   % Shared by the ncread (single) and netcdf.getVar (batch) read paths so
   % both produce byte-identical blocks.

   % Collapse to cells x time. ncread/getVar return [ni nj (24) ndays];
   % squeeze drops singleton dims, so reshape from the requested cell count.
   ncells = prod(count(1:2));
   data = reshape(data, ncells, []);

   % No-data handling precedes unit conversion. MAR uses signed ~1e34/1e36
   % fill values; both signs are unambiguous. Do not apply the legacy >=999
   % mass-flux cut: raw RUH values above 999 sum to native no-delay RU2, and
   % paired positive/negative SMBH pulses sum to native daily SMB. Treating
   % only the positive pulse as missing corrupts the daily mass balance.
   data(abs(data) >= 1e30) = NaN;

   % Standard unit conversions (legacy readMar3p11 table). The unit STRING is
   % relabelled once at the top level by convertUnits (kept in lockstep here).
   switch units
      case 'C'
         data = data + 273.15;
      case 'g/kg'
         data = data / 1000;
      case 'mmWE/h'
         data = data / 1000;
      case 'mmWE/day'
         data = data / 1000;
      case 'hPa'
         data = data * 100;
   end
end

function units = convertUnits(units)
   %CONVERTUNITS Relabel a native MAR unit string to its standard form.
   % Mirrors the value scaling in convertSlab (data and label stay in sync).
   switch units
      case 'C'
         units = 'K';
      case 'g/kg'
         units = 'kg/kg';
      case 'mmWE/h'
         units = 'mWE/h';
      case 'mmWE/day'
         units = 'mWE/day';
      case 'hPa'
         units = 'Pa';
   end
end

function Time = marTime(filename, dims)
   %MARTIME UTC datetime axis from a MAR yearly file (interval-start stamped).
   %
   % Yearly files start Jan 1 00:00 UTC; hourly variables have 24 samples per
   % day, daily variables one. The hourly axis t0+(0:23)h is INTERVAL-START
   % stamped, matching the icemodel forcing time convention (the [t, t+dt)
   % label = interval start; see icemodel/+icemodel/+forcing/README.md "Time
   % convention"). No shift is applied. The hourly/daily branch is decided
   % from the time-dimension layout, not the (collapsed) data, so it is
   % independent of the single/batch read path.
   yyyy = marFileYear(filename);
   ndays = dims(end);
   t0 = datetime(yyyy, 1, 1, 'TimeZone', 'UTC');
   if numel(dims) >= 4 && dims(3) == 24
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
