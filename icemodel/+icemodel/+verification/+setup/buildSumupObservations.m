function [observations, metadata] = buildSumupObservations(point, kwargs)
   %BUILDSUMUPOBSERVATIONS Convert SUMup firn records to verification targets.
   %
   %  [observations, metadata] = ...
   %     icemodel.verification.setup.buildSumupObservations([lat lon])
   %  [observations, metadata] = ...
   %     icemodel.verification.setup.buildSumupObservations([lat lon], ...
   %        source_dir=..., radius_km=..., startdate=..., enddate=...)
   %
   %  Reads the cached SUMup density / accumulation / subsurface-temperature
   %  files and returns the firn observations nearest one [lat lon] point as a
   %  verification-target struct. SUMup is the firn observation source; it
   %  subsumes Humphrey 2012 subsurface temperature and GreenTRACS
   %  accumulation. Mirrors buildEsmSnowmipObservations: a reusable per-point
   %  observation builder used by importSumup during staging.
   %
   %  Output target schema (verification timeseries / profile bundle):
   %    observations.format   = "subsurface_profile_bundle"
   %    observations.density            depth-indexed profile TABLE
   %                                    (depth, density, error, datetime)
   %    observations.subsurface_temperature  datetime-indexed TIMETABLE T(z,t)
   %    observations.accumulation       period-indexed TABLE with start_date /
   %                                    end_date datetimes
   %  Each sub-bundle is present only when the corresponding SUMup variable
   %  file is in the cache and has a record within radius_km of the point.
   %
   %  The bundle struct is a WRAPPER for three heterogeneous profile tables
   %  (density rho(z), subsurface temperature T(z,t), accumulation/SMB), not a
   %  storage choice; the three tables carry different indexing axes (depth,
   %  time, period) and cannot share one table. The generic name
   %  "subsurface_profile_bundle" (not "firn_profile_bundle") covers ablation
   %  sites too, where the bare-ice/seasonal-snow column is not firn.
   %
   %  TIME AXIS. SUMup stores a numeric `timestamp` (days since 1900-01-01).
   %  This builder converts it to real UTC datetimes:
   %    - subsurface_temperature is a time series, returned as a TIMETABLE
   %      indexed by the measurement datetime (row time `Time`);
   %    - density is a depth profile, returned as a TABLE with an added
   %      `datetime` column carrying the profile timestamp;
   %    - accumulation is a period observation, returned as a TABLE with
   %      `start_date` / `end_date` converted to datetimes (the integer
   %      days-since-1900 columns are replaced in place).
   %
   %  Inputs
   %    point : [lat lon]  WGS84 site coordinates.
   %
   %  Name-value
   %    source_dir : string (default data/verification/sumup)
   %    radius_km  : double (default 7.5)  point-selection radius
   %    startdate  : datetime or "" (default "")
   %    enddate    : datetime or "" (default "")
   %
   %  Returns
   %    observations : struct  firn-observation target bundle
   %    metadata     : struct  provenance + which variables were found
   %
   %  Role
   %    Reusable per-point SUMup observation builder, symmetric with
   %    buildSumupForcing. Used by importSumup. The low-level SUMup file
   %    parsing is intentionally isolated here so importSumup stays a staging
   %    orchestrator.
   %
   %  The concrete SUMup parsing targets the real 2025 release Greenland files
   %  (grouped NetCDF: /DATA + /METADATA). When the cache is missing, fetchSumup
   %  (strict=true) errors with the retrieval banner rather than fabricating
   %  records.
   %
   %  icemodel.verification.setup.importSumup,
   %  icemodel.verification.setup.fetchSumup

   arguments
      point (1, 2) double
      kwargs.source_dir (1, 1) string = ""
      kwargs.radius_km (1, 1) double {mustBePositive} = 7.5
      kwargs.startdate = ""
      kwargs.enddate   = ""
   end

   % Resolve and verify the cache (fetch is the single source of truth for
   % "are the SUMup files present?"). strict=true errors with the retrieval
   % banner when the cache is empty, so this builder never fabricates records.
   source_dir = icemodel.verification.setup.fetchSumup( ...
      cache_dir=resolveCacheDir(kwargs.source_dir), strict=true);

   % Resolve the optional comparison window.
   has_start = ~strcmp(string(kwargs.startdate), "");
   has_end = ~strcmp(string(kwargs.enddate), "");
   if has_start ~= has_end
      error('icemodel:verification:buildSumupObservations:halfWindow', ...
         'startdate and enddate must be provided together')
   end

   % Read each SUMup variable group nearest the point. Each reader returns a
   % table and a per-variable provenance note, or empty when no record falls
   % within radius_km. The 2025 release stores the measured channel as `SMB`
   % (the firn accumulation/SMB observation source).
   [density, density_note] = readSumupVariable(source_dir, "density", ...
      point, kwargs.radius_km, kwargs.startdate, kwargs.enddate);
   [temperature, temp_note] = readSumupVariable(source_dir, "temperature", ...
      point, kwargs.radius_km, kwargs.startdate, kwargs.enddate);
   [accumulation, accum_note] = readSumupVariable(source_dir, "SMB", ...
      point, kwargs.radius_km, kwargs.startdate, kwargs.enddate);

   observations = struct( ...
      'format', 'subsurface_profile_bundle', ...
      'density', density, ...
      'subsurface_temperature', temperature, ...
      'accumulation', accumulation);

   metadata = icemodel.verification.setup.metadataStruct({ ...
      'observation_source', 'SUMup 2025 release (NSIDC G02288)'
      'point_lat_wgs84', point(1)
      'point_lon_wgs84', point(2)
      'selection_radius_km', kwargs.radius_km
      'density_note', density_note
      'subsurface_temperature_note', temp_note
      'accumulation_note', accum_note});
end

%% Local helpers
function cache_dir = resolveCacheDir(source_dir)
   %RESOLVECACHEDIR Resolve the SUMup cache dir, defaulting to the standard one.
   if strlength(source_dir) > 0
      cache_dir = source_dir;
   else
      cache_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'sumup'));
   end
end

function [record, note] = readSumupVariable(source_dir, variable, point, ...
      radius_km, startdate, enddate)
   %READSUMUPVARIABLE Read the nearest SUMup records for one variable group.
   %
   % Locates the SUMup 2025 Greenland file for the requested variable group
   % (density / temperature / SMB), reads the gridded /DATA group into a flat
   % table, selects the records within radius_km of the point, optionally
   % windows them by timestamp, resolves the name_key into the human-readable
   % core/site name, and returns the selection plus a provenance note. Returns
   % [] with an explanatory note when the variable file is absent or no record
   % falls within range.
   %
   % SUMup 2025 NetCDF layout: each release file holds two groups, /DATA and
   % /METADATA. /DATA exposes 1-D variables along measurement_id (latitude,
   % longitude, timestamp/start_date, the measured channel, depth columns, the
   % per-measurement error, and the name_key/reference_key/method_key links).
   % /METADATA maps name_key->name, method_key->method, reference_key->
   % reference_short. The reader is scoped to the Greenland files via the glob.

   record = [];
   nc = dir(fullfile(source_dir, sprintf('*%s*greenland*.nc', variable)));

   if isempty(nc)
      note = sprintf('no SUMup %s greenland file in cache', variable);
      return
   end

   file = fullfile(nc(1).folder, nc(1).name);
   tbl = readSumupNetcdf(file, variable);
   [record, note] = selectNearest( ...
      tbl, point, radius_km, startdate, enddate, variable);
end

function tbl = readSumupNetcdf(file, variable)
   %READSUMUPNETCDF Read the /DATA group of a SUMup 2025 file into a table.
   %
   % Reads the coordinate, time, value, depth, error and key channels that
   % apply to the requested variable group, then attaches the resolved
   % name_key -> name string from /METADATA. Depth/value column sets differ per
   % group (density carries start/stop/midpoint depths; temperature a single
   % depth; SMB a start/end date span), so the channel list is per-variable.
   variable = char(variable);

   % Channels common to all SUMup variable groups.
   common = ["latitude", "longitude", "elevation", ...
      "name_key", "reference_key", "method_key", "error"];

   switch lower(variable)
      case "density"
         channels = [common, "density", ...
            "start_depth", "stop_depth", "midpoint", "timestamp"];
      case "temperature"
         channels = [common, "temperature", ...
            "depth", "duration", "timestamp"];
      case "smb"
         channels = [common, "smb", ...
            "start_date", "end_date", "start_year", "end_year"];
      otherwise
         channels = common;
   end

   tbl = table();
   for k = 1:numel(channels)
      name = char(channels(k));
      v = ncread(file, ['/DATA/' name]);
      tbl.(name) = double(v(:));
   end

   % Resolve name_key -> name (core / location label) from /METADATA. The
   % /METADATA/name char matrix is fixed-width: each name is padded to the
   % column width with trailing NUL (char 0) and/or space. Strip both so the
   % name is the bare label (e.g. "s5", not "s5" + 31 pad chars). deblank
   % removes trailing whitespace; char(0) is removed explicitly first.
   meta_keys = double(ncread(file, '/METADATA/name_key'));
   meta_names = string(ncread(file, '/METADATA/name')');
   meta_names = strtrim(erase(meta_names, char(0)));
   [tf, loc] = ismember(tbl.name_key, meta_keys);
   names = strings(height(tbl), 1);
   names(tf) = meta_names(loc(tf));
   tbl.name = names;
end

function [record, note] = selectNearest(tbl, point, radius_km, ...
      startdate, enddate, variable)
   %SELECTNEAREST Keep table rows within radius_km of the point (optional window).
   record = [];

   proj = icemodel.forcing.helpers.psnProjection();
   [px, py] = projfwd(proj, point(1), point(2));
   [rx, ry] = projfwd(proj, tbl.latitude, tbl.longitude);
   d_km = hypot(rx - px, ry - py) / 1000;
   keep = d_km <= radius_km;

   % Optional temporal window. Density/temperature carry `timestamp`; SMB
   % carries a `start_date` span. Both are int64 days since 1900-01-01.
   time_name = "";
   if ismember("timestamp", string(tbl.Properties.VariableNames))
      time_name = "timestamp";
   elseif ismember("start_date", string(tbl.Properties.VariableNames))
      time_name = "start_date";
   end
   if time_name ~= "" && ~strcmp(string(startdate), "")
      t = datetime(1900, 1, 1) + days(tbl.(time_name));
      t.TimeZone = 'UTC';
      t1 = icemodel.verification.setup.ensureUtc(startdate);
      t2 = icemodel.verification.setup.ensureUtc(enddate);
      keep = keep & t >= t1 & t <= t2;
   end

   if ~any(keep)
      note = sprintf('no SUMup %s record within %.1f km', variable, radius_km);
      return
   end
   record = tbl(keep, :);

   % Attach real UTC datetimes from the SUMup days-since-1900 integers, and
   % shape each variable to its natural index axis (see header). Temperature is
   % a time series -> datetime-indexed timetable; density is a depth profile ->
   % table + profile datetime column; accumulation is a period observation ->
   % table with start_date/end_date as datetimes.
   record = attachDatetimes(record, variable);

   note = sprintf('%d SUMup %s record(s) within %.1f km of %d profile(s)', ...
      nnz(keep), variable, radius_km, numel(unique(record.name_key)));
end

function record = attachDatetimes(record, variable)
   %ATTACHDATETIMES Convert SUMup days-since-1900 integers to UTC datetimes.
   epoch = datetime(1900, 1, 1, 'TimeZone', 'UTC');
   switch lower(char(variable))
      case "temperature"
         % Subsurface temperature is a time series: index by measurement time.
         t = epoch + days(record.timestamp);
         record.timestamp = [];
         record = table2timetable(record, 'RowTimes', t);
      case "density"
         % Depth profile: keep the table, add the profile datetime column.
         record.datetime = epoch + days(record.timestamp);
      case "smb"
         % Period observation: replace the integer span with datetimes.
         record.start_date = epoch + days(record.start_date);
         record.end_date = epoch + days(record.end_date);
   end
end
