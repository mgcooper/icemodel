function [observations, metadata] = buildSumupObservations(point, kwargs)
   %BUILDSUMUPOBSERVATIONS Convert SUMup firn records to verification targets.
   %
   %  [observations, metadata] = ...
   %     icemodel.verification.setup.buildSumupObservations([lat lon])
   %  [observations, metadata] = ...
   %     icemodel.verification.setup.buildSumupObservations([lat lon], ...
   %        source_dir=..., radius_km=..., startdate=..., enddate=...)
   %
   %  Reads the cached SUMup density / SMB / subsurface-temperature files and
   %  returns the firn observations nearest one [lat lon] point as a
   %  verification-target struct. Mirrors buildEsmSnowmipObservations: a reusable
   %  per-point observation builder used by importSumup during staging.
   %
   %  Output target schema (verification timeseries / profile bundle):
   %    observations.format   = "subsurface_profile_bundle"
   %    observations.density            depth-indexed profile TABLE
   %                                    (depth, density, error, datetime)
   %    observations.subsurface_temperature  datetime-indexed TIMETABLE T(z,t)
   %    observations.smb                period-indexed TABLE with start_date /
   %                                    end_date datetimes (signed surface mass
   %                                    balance: accumulation OR ablation)
   %    metadata                        per-bundle provenance including exact
   %                                    raw/unique/duplicate-removed row counts
   %  Each sub-bundle is present only when the corresponding SUMup variable
   %  file is in the cache and has a record within radius_km of the point.
   %
   %  The bundle struct is a WRAPPER for three heterogeneous profile tables
   %  (density rho(z), subsurface temperature T(z,t), smb), not a storage choice;
   %  the three tables carry different indexing axes (depth, time, period) and
   %  cannot share one table. The generic name "subsurface_profile_bundle" (not
   %  "firn_profile_bundle") covers ablation sites too, where the bare-ice/
   %  seasonal-snow column is not firn.
   %
   %  TIME AXIS. SUMup stores a numeric `timestamp` (days since 1900-01-01).
   %  This builder converts it to real UTC datetimes:
   %    - subsurface_temperature is a time series, returned as a TIMETABLE
   %      indexed by the measurement datetime (row time `Time`);
   %    - density is a depth profile, returned as a TABLE with an added
   %      `datetime` column carrying the profile timestamp;
   %    - smb is a period observation, returned as a TABLE with
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
   %    Reusable per-point SUMup observation builder used by importSumup. The
   %    low-level SUMup file parsing is intentionally isolated here so
   %    importSumup stays a staging orchestrator.
   %
   %  The concrete SUMup parsing targets the real 2025 release Greenland files
   %  (grouped NetCDF: /DATA + /METADATA). When the cache is missing, fetchSumup
   %  (strict=true) errors with the retrieval banner rather than fabricating
   %  records. Exact scientific duplicates are removed after spatial/window
   %  selection and before this datetime shaping; aliases/elevation retain the
   %  first selected source row as provenance.
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

   % Reject malformed public windows before resolving or probing the SUMup cache.
   [window_start, window_end, has_window] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);

   % Resolve and verify the cache (fetch is the single source of truth for
   % "are the SUMup files present?"). strict=true errors with the retrieval
   % banner when the cache is empty, so this builder never fabricates records.
   source_dir = icemodel.verification.setup.fetchSumup( ...
      cache_dir=icemodel.verification.setup.sumupCacheDir(kwargs.source_dir), ...
      strict=true);

   % Use the normalized optional comparison window for all three source readers.
   if ~has_window
      window_start = "";
      window_end = "";
   end

   % Read each SUMup variable group nearest the point. Each reader returns a
   % table and a per-variable provenance note, or empty when no record falls
   % within radius_km. The 2025 release stores the measured channel as `SMB`
   % (the firn accumulation/SMB observation source).
   [density, density_note, density_counts] = ...
      readSumupVariable(source_dir, "density", ...
      point, kwargs.radius_km, window_start, window_end);
   [temperature, temp_note, temperature_counts] = ...
      readSumupVariable(source_dir, "temperature", ...
      point, kwargs.radius_km, window_start, window_end);
   % SMB (surface mass balance) - the third obs axis. NOT "accumulation": SUMup
   % spans accumulation AND ablation zones, so this quantity is signed SMB
   % (positive net accumulation / negative net ablation), not accumulation per se.
   [smb, smb_note, smb_counts] = readSumupVariable(source_dir, "SMB", ...
      point, kwargs.radius_km, window_start, window_end);

   observations = struct( ...
      'format', 'subsurface_profile_bundle', ...
      'density', density, ...
      'subsurface_temperature', temperature, ...
      'smb', smb);
   observations = icemodel.verification.setup.stampArtifactMetadata( ...
      observations);
   % Reported uncertainties have the same physical units as their measured
   % channel; the generic `error` name cannot supply that context by itself.
   observations.density = stampObservedUncertaintyUnits( ...
      observations.density, "density");
   observations.subsurface_temperature = stampObservedTemperatureUnits( ...
      observations.subsurface_temperature);
   observations.subsurface_temperature = stampObservedUncertaintyUnits( ...
      observations.subsurface_temperature, "subsurface_temperature");
   observations.smb = stampObservedSmbUnits(observations.smb);
   observations.smb = stampObservedUncertaintyUnits( ...
      observations.smb, "smb");

   % Record the actual observation coverage so unbounded imports can still
   % write inspectable manifest periods without imposing a hidden fixed window.
   [period_start, period_end] = observationPeriod(observations);

   metadata = icemodel.verification.setup.metadataStruct({ ...
      'observation_source', 'SUMup 2025 release (NSIDC G02288)'
      'point_lat_wgs84', point(1)
      'point_lon_wgs84', point(2)
      'selection_radius_km', kwargs.radius_km
      'observation_period_start', period_start
      'observation_period_end', period_end
      'density_note', density_note
      'density_raw_rows', density_counts.raw_rows
      'density_unique_rows', density_counts.unique_rows
      'density_duplicate_rows_removed', density_counts.removed_rows
      'subsurface_temperature_note', temp_note
      'subsurface_temperature_raw_rows', temperature_counts.raw_rows
      'subsurface_temperature_unique_rows', temperature_counts.unique_rows
      'subsurface_temperature_duplicate_rows_removed', ...
      temperature_counts.removed_rows
      'smb_note', smb_note
      'smb_raw_rows', smb_counts.raw_rows
      'smb_unique_rows', smb_counts.unique_rows
      'smb_duplicate_rows_removed', smb_counts.removed_rows});
end

function [t1, t2] = observationPeriod(observations)
   %OBSERVATIONPERIOD Return min/max datetimes across SUMup observation tables.
   t = NaT(0, 1, 'TimeZone', 'UTC');
   fields = string(fieldnames(observations));
   for k = 1:numel(fields)
      value = observations.(fields(k));
      if istimetable(value)
         t = appendTimes(t, value.Time);
      elseif istable(value)
         if ismember("datetime", string(value.Properties.VariableNames))
            t = appendTimes(t, value.datetime);
         end
         if all(ismember(["start_date", "end_date"], ...
               string(value.Properties.VariableNames)))
            t = appendTimes(t, value.start_date);
            t = appendTimes(t, value.end_date);
         end
      end
   end

   if isempty(t)
      t1 = NaT(1, 1, 'TimeZone', 'UTC');
      t2 = NaT(1, 1, 'TimeZone', 'UTC');
   else
      t1 = min(t);
      t2 = max(t);
   end
end

function out = appendTimes(out, values)
   %APPENDTIMES Add finite datetimes to an accumulated UTC datetime vector.
   values = values(:);
   values = values(~isnat(values));
   if isempty(values)
      return
   end
   if isempty(values.TimeZone)
      values.TimeZone = 'UTC';
   end
   values.TimeZone = 'UTC';
   out = [out; values];
end

function smb = stampObservedSmbUnits(smb)
   %STAMPOBSERVEDSMBUNITS Override forcing-rate metadata for SUMup SMB.
   if ~istable(smb) || ~ismember("smb", string(smb.Properties.VariableNames))
      return
   end

   idx = find(string(smb.Properties.VariableNames) == "smb", 1);
   smb.Properties.VariableUnits{idx} = 'm w.e.';
   smb.Properties.VariableDescriptions{idx} = ...
      'accumulated surface mass balance over the observation interval';
   if isprop(smb.Properties.CustomProperties, "StandardNames")
      smb.Properties.CustomProperties.StandardNames(idx) = "";
   end
end

function temperature = stampObservedTemperatureUnits(temperature)
   %STAMPOBSERVEDTEMPERATUREUNITS Keep SUMup profile temperatures in Celsius.
   if ~(istable(temperature) || istimetable(temperature)) ...
         || ~ismember("subsurface_temperature", ...
         string(temperature.Properties.VariableNames))
      return
   end

   idx = find(string(temperature.Properties.VariableNames) ...
      == "subsurface_temperature", 1);
   temperature.Properties.VariableUnits{idx} = 'degC';
   temperature.Properties.VariableDescriptions{idx} = ...
      'observed subsurface temperature';
   if isprop(temperature.Properties.CustomProperties, "StandardNames")
      temperature.Properties.CustomProperties.StandardNames(idx) = "";
   end
end

function artifact = stampObservedUncertaintyUnits(artifact, measured_name)
   %STAMPOBSERVEDUNCERTAINTYUNITS Apply the measured channel's units to error.
   names = string.empty(1, 0);
   if istable(artifact) || istimetable(artifact)
      names = string(artifact.Properties.VariableNames);
   end
   if ~all(ismember([measured_name, "error"], names))
      return
   end

   measured_idx = find(names == measured_name, 1);
   error_idx = find(names == "error", 1);
   artifact.Properties.VariableUnits{error_idx} = ...
      artifact.Properties.VariableUnits{measured_idx};
   artifact.Properties.VariableDescriptions{error_idx} = ...
      sprintf('reported uncertainty in %s', measured_name);
end

%% Local helpers
function [record, note, counts] = readSumupVariable(source_dir, variable, ...
      point, ...
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
   counts = struct('raw_rows', 0, 'unique_rows', 0, 'removed_rows', 0);
   nc = dir(fullfile(source_dir, sprintf('*%s*greenland*.nc', variable)));

   if isempty(nc)
      note = sprintf([ ...
         'no SUMup %s greenland file in cache; rows raw=0 unique=0 ' ...
         'duplicates_removed=0'], variable);
      return
   end

   file = fullfile(nc(1).folder, nc(1).name);
   tbl = readSumupNetcdf(file, variable);
   [record, note, counts] = selectNearest( ...
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

   % measurement_id is the native NetCDF dimension rather than a coordinate
   % variable. Preserve its one-based source row index as provenance, but keep it
   % outside scientific identity because duplicate groups have distinct ids.
   tbl.measurement_id = (1:height(tbl))';

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

function [record, note, counts] = selectNearest(tbl, point, radius_km, ...
      startdate, enddate, variable)
   %SELECTNEAREST Keep table rows within radius_km of the point (optional window).
   record = [];
   counts = struct('raw_rows', 0, 'unique_rows', 0, 'removed_rows', 0);

   proj = icemodel.forcing.helpers.psnProjection();
   [px, py] = projfwd(proj, point(1), point(2));
   [rx, ry] = projfwd(proj, tbl.latitude, tbl.longitude);
   d_km = hypot(rx - px, ry - py) / 1000;
   keep = d_km <= radius_km;

   % Optional temporal window. Point observations use timestamp containment.
   % SMB period observations are accumulated over their full start/end interval,
   % so explicit windows retain only records fully contained in that window.
   names = string(tbl.Properties.VariableNames);
   if ismember("timestamp", names) && ~strcmp(string(startdate), "")
      t = sumupDatetime(tbl.timestamp);
      t1 = icemodel.verification.setup.ensureUtc(startdate);
      t2 = icemodel.verification.setup.ensureUtc(enddate);
      keep = keep & t >= t1 & t <= t2;
   elseif all(ismember(["start_date", "end_date"], names)) ...
         && ~strcmp(string(startdate), "")
      t_start = sumupDatetime(tbl.start_date);
      t_end = sumupDatetime(tbl.end_date);
      t1 = icemodel.verification.setup.ensureUtc(startdate);
      t2 = icemodel.verification.setup.ensureUtc(enddate);
      keep = keep & t_start >= t1 & t_end <= t2;
   end

   if ~any(keep)
      note = sprintf([ ...
         'no SUMup %s record within %.1f km; rows raw=0 unique=0 ' ...
         'duplicates_removed=0'], variable, radius_km);
      return
   end
   record = tbl(keep, :);

   % Scientific identity is resolved while every source field is still numeric.
   % Keep the pre-dedup profile count for the existing selection provenance.
   n_profiles = numel(unique(record.name_key));
   [record, counts.raw_rows, counts.unique_rows, counts.removed_rows] = ...
      icemodel.verification.setup.deduplicateSumupRecords(record, variable);

   % Attach real UTC datetimes from the SUMup days-since-1900 integers, and
   % shape each variable to its natural index axis (see header). Temperature is
   % a time series -> datetime-indexed timetable; density is a depth profile ->
   % table + profile datetime column; SMB is a period observation ->
   % table with start_date/end_date as datetimes.
   record = attachDatetimes(record, variable);

   note = sprintf([ ...
      '%d SUMup %s record(s) within %.1f km of %d profile(s); ' ...
      'rows raw=%d unique=%d duplicates_removed=%d'], ...
      counts.unique_rows, variable, radius_km, n_profiles, ...
      counts.raw_rows, counts.unique_rows, counts.removed_rows);
end

function record = attachDatetimes(record, variable)
   %ATTACHDATETIMES Convert SUMup days-since-1900 integers to UTC datetimes.
   switch lower(char(variable))
      case "temperature"
         % Subsurface temperature is a time series: index by measurement time.
         t = sumupDatetime(record.timestamp);
         record.timestamp = [];
         record = table2timetable(record, 'RowTimes', t);
         record = renamevars(record, "temperature", ...
            "subsurface_temperature");
         record.Properties.DimensionNames{1} = 'Time';
      case "density"
         % Depth profile: keep the table, add the profile datetime column.
         record.datetime = sumupDatetime(record.timestamp);
         record.timestamp = [];
      case "smb"
         % Period observation: replace the integer span with datetimes.
         record.start_date = sumupDatetime(record.start_date);
         record.end_date = sumupDatetime(record.end_date);
   end
end

function t = sumupDatetime(days_since_1900)
   %SUMUPDATETIME Convert SUMup day offsets to UTC datetime values.
   epoch = datetime(1900, 1, 1, 'TimeZone', 'UTC');
   t = epoch + days(days_since_1900);
end
