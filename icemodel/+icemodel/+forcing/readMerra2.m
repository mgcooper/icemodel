function [data, units, Time] = readMerra2(filename, varname, kwargs)
   %READMERRA2 Read one MERRA-2 variable in standard units.
   %
   %  [data, units, Time] = icemodel.forcing.readMerra2(filename, varname)
   %  [data, units, Time] = ... readMerra2(_, start=[i j], count=[ni nj])
   %  [blocks, units, Time] = ... readMerra2(_, slabs={[i j;ni nj], ...})
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
   %  slabs        - optional cell array of [start; count] 2x2 hyperslab specs
   %                 ({[i j; ni nj], ...}). When given, the file is OPENED ONCE
   %                 and every listed hyperslab is read from the same open file,
   %                 returning a cell array of blocks (one per slab) instead of a
   %                 single matrix - the batch path that extracts many points
   %                 from one daily file without re-opening it per point.
   %                 start/count are ignored when slabs is given.
   %
   % Outputs
   %  data  - (ncells x ntime) double in standard units, native grid order;
   %          OR, when slabs is given, a cell array {(ncells x ntime), ...}
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
      kwargs.slabs cell = {}
   end

   [dims, native_units] = ncVarInfo(filename, varname);

   if isempty(kwargs.slabs)
      [start, count] = slabWindow(dims, kwargs.start, kwargs.count);
      data = convertSlab(double(squeeze( ...
         ncread(filename, varname, start, count))), native_units, count);
   else
      % Batch path: open once, read every listed hyperslab from the same ncid.
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
      Time = merraTime(filename);
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
   %SLABWINDOW Assemble the read window: requested spatial hyperslab, all times.
   if isempty(kstart)
      start = ones(1, numel(dims));
      count = dims;
   else
      start = [kstart, ones(1, numel(dims) - 2)];
      count = [kcount, dims(3:end)];
   end
end

function data = convertSlab(data, units, count)
   %CONVERTSLAB Reshape to cells x time, mask fill, convert to standard units.
   % Shared by the ncread (single) and netcdf.getVar (batch) read paths.

   % Collapse to cells x time (cells flattened in native [lon lat] order).
   ncells = prod(count(1:2));
   data = reshape(data, ncells, []);

   % MERRA-2 _FillValue (~1e15) over masked tiles (e.g. glacier-tile
   % variables off-ice) -> NaN, so it never corrupts a downstream mean/sum.
   data(data >= 1e14) = NaN;

   % Standard unit conversions (shared reader family). Note kg m-2 s-1 (a
   % flux rate) converts to mWE/h, but kg m-2 (a store, e.g. swe) does not.
   % The unit STRING is relabelled once at the top level by convertUnits.
   switch units
      case 'kg m-2 s-1'
         data = data * 3600 / 1000;   % -> meters water equivalent per hour
      case 'C'
         data = data + 273.15;
      case 'g/kg'
         data = data / 1000;
      case 'hPa'
         data = data * 100;
   end
end

function units = convertUnits(units)
   %CONVERTUNITS Relabel a native MERRA-2 unit string to its standard form.
   % Mirrors the value scaling in convertSlab (data and label stay in sync).
   switch units
      case 'kg m-2 s-1'
         units = 'mWE/h';
      case 'C'
         units = 'K';
      case 'g/kg'
         units = 'kg/kg';
      case 'hPa'
         units = 'Pa';
      case 'W m-2'
         units = 'W/m2';
      case 'm s-1'
         units = 'm/s';
   end
end

function Time = merraTime(filename)
   %MERRATIME UTC interval-START datetime axis from the MERRA-2 'time' variable.
   % MERRA-2 stores time as "minutes since <yyyy-mm-dd hh:mm:ss>" and posts tavg
   % samples at the bin CENTER (:30 hourly, 1:30 three-hourly). The icemodel
   % forcing convention labels an averaged interval at its START t (the
   % [t, t+dt) integration; see icemodel/+icemodel/+forcing/README.md "Time
   % convention"), matching readMar3p11. We therefore shift the native bin
   % centers back by half a step so MERRA aligns with the other sources.
   t = double(ncread(filename, 'time'));
   t_units = ncreadatt(filename, 'time', 'units');
   tok = regexp(t_units, 'minutes since (\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2})', ...
      'tokens', 'once');
   assert(~isempty(tok), 'unexpected MERRA time units: %s', t_units)
   t0 = datetime(tok{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
   Time = t0 + minutes(t);
   if numel(t) > 1
      Time = Time - minutes(mode(diff(t))) / 2;   % bin center -> interval start
   end
end
