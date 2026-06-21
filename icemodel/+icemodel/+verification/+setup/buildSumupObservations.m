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
   %    observations.format   = "firn_profile_bundle"
   %    observations.density            profile table (depth, density, error)
   %    observations.subsurface_temperature  T(z,t) records
   %    observations.accumulation       SMB / accumulation records
   %  Each sub-bundle is present only when the corresponding SUMup variable
   %  file is in the cache and has a record within radius_km of the point.
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
   %  NOTE (data-gated): the concrete SUMup column/variable parsing is
   %    finalized against the real 2024 release files; until the gitignored
   %    cache is populated (see fetchSumup), this builder errors with the
   %    fetch retrieval banner rather than fabricating records.
   %
   % See also: icemodel.verification.setup.buildSumupForcing,
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
   % timetable/profile and a per-variable provenance note, or empty when no
   % record falls within radius_km.
   [density, density_note] = readSumupVariable(source_dir, "density", ...
      point, kwargs.radius_km, kwargs.startdate, kwargs.enddate);
   [temperature, temp_note] = readSumupVariable(source_dir, "temperature", ...
      point, kwargs.radius_km, kwargs.startdate, kwargs.enddate);
   [accumulation, accum_note] = readSumupVariable(source_dir, "accumulation", ...
      point, kwargs.radius_km, kwargs.startdate, kwargs.enddate);

   observations = struct( ...
      'format', 'firn_profile_bundle', ...
      'density', density, ...
      'subsurface_temperature', temperature, ...
      'accumulation', accumulation);

   metadata = icemodel.verification.setup.metadataStruct({ ...
      'observation_source', 'SUMup 2024 release (NSIDC G02288)'
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
   %READSUMUPVARIABLE Read the nearest SUMup record for one variable group.
   %
   % Locates the SUMup file for the requested variable (NetCDF preferred, CSV
   % fallback), selects the records within radius_km of the point, optionally
   % windows them, and returns the selection plus a provenance note. Returns
   % [] with an explanatory note when the variable file is absent or no
   % record falls within range.
   %
   % The SUMup release ships density / accumulation / temperature as separate
   % point collections keyed by lat/lon/timestamp; the parsing below is kept
   % deliberately small and is exercised against the real files once the
   % gitignored cache is populated. Until then fetchSumup(strict=true) upstream
   % guarantees this is only reached when a file exists.

   record = [];
   nc = dir(fullfile(source_dir, sprintf('*%s*.nc', variable)));
   csv = dir(fullfile(source_dir, sprintf('*%s*.csv', variable)));

   if ~isempty(nc)
      file = fullfile(nc(1).folder, nc(1).name);
      [record, note] = selectNearest( ...
         readSumupNetcdf(file), point, radius_km, startdate, enddate, variable);
   elseif ~isempty(csv)
      file = fullfile(csv(1).folder, csv(1).name);
      [record, note] = selectNearest( ...
         readtable(file), point, radius_km, startdate, enddate, variable);
   else
      note = sprintf('no SUMup %s file in cache', variable);
   end
end

function tbl = readSumupNetcdf(file)
   %READSUMUPNETCDF Read a SUMup NetCDF point collection into a table.
   %
   % SUMup point collections expose 1-D coordinate variables (latitude,
   % longitude, timestamp) and the measured variable(s). Read the common
   % coordinates and every numeric data variable into one flat table.
   info = ncinfo(file);
   names = string({info.Variables.name});
   tbl = table();
   for k = 1:numel(names)
      v = ncread(file, names(k));
      if isvector(v)
         tbl.(matlab.lang.makeValidName(char(names(k)))) = double(v(:));
      end
   end
end

function [record, note] = selectNearest(tbl, point, radius_km, ...
      startdate, enddate, variable)
   %SELECTNEAREST Keep table rows within radius_km of the point (optional window).
   record = [];
   names = string(tbl.Properties.VariableNames);
   lat_name = firstMatch(names, ["latitude", "lat"]);
   lon_name = firstMatch(names, ["longitude", "lon"]);
   if lat_name == "" || lon_name == ""
      note = sprintf('SUMup %s file has no lat/lon columns', variable);
      return
   end

   proj = icemodel.forcing.helpers.psnProjection();
   [px, py] = projfwd(proj, point(1), point(2));
   [rx, ry] = projfwd(proj, tbl.(lat_name), tbl.(lon_name));
   d_km = hypot(rx - px, ry - py) / 1000;
   keep = d_km <= radius_km;

   % Optional temporal window on a timestamp column when present.
   time_name = firstMatch(names, ["timestamp", "time", "date"]);
   if time_name ~= "" && ~strcmp(string(startdate), "")
      t = tbl.(time_name);
      if ~isdatetime(t)
         t = datetime(t, 'ConvertFrom', 'datenum');
      end
      t1 = icemodel.verification.setup.ensureUtc(startdate);
      t2 = icemodel.verification.setup.ensureUtc(enddate);
      if isempty(t.TimeZone)
         t.TimeZone = 'UTC';
      end
      keep = keep & t >= t1 & t <= t2;
   end

   if ~any(keep)
      note = sprintf('no SUMup %s record within %.1f km', variable, radius_km);
      return
   end
   record = tbl(keep, :);
   note = sprintf('%d SUMup %s record(s) within %.1f km', ...
      nnz(keep), variable, radius_km);
end

function name = firstMatch(names, candidates)
   %FIRSTMATCH First candidate that appears in names (case-insensitive), else "".
   name = "";
   lnames = lower(names);
   for c = candidates
      hit = find(lnames == lower(c), 1);
      if ~isempty(hit)
         name = names(hit);
         return
      end
   end
end
