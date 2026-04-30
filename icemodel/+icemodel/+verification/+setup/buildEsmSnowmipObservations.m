function [observations, metadata] = buildEsmSnowmipObservations(sitename, kwargs)
   %BUILDESMSNOWMIPOBSERVATIONS Convert ESM-SnowMIP obs NetCDF to verification targets.
   %
   %  [observations, metadata] = buildEsmSnowmipObservations(sitename)
   %  [observations, metadata] = buildEsmSnowmipObservations(sitename, ...
   %     source_dir=..., startdate=..., enddate=...)
   %
   %  Reads one ESM-SnowMIP obs_insitu_<sitename>_*.nc file and returns
   %  the observed snow / surface variables as a verification-target
   %  timetable. Variable mapping mirrors the existing schema used by
   %  the cdp / wfj smoke fixtures so the broader verification harness
   %  (comparecase, plotcase, runIcemodelSnowCandidate adapter) does
   %  not need a per-site branch.
   %
   %  Variable mapping (ESM-SnowMIP -> verification target):
   %     snd_auto / snd_man -> snow_depth_m   [m]
   %     snw_auto / snw_man -> swe_kg_m2      [kg m-2]
   %     ts                 -> surface_temp_C [C]
   %     tsl(:, k)          -> soil_temp_<k>_C [C] for each soil depth
   %
   %  When the optional automatic channel is missing for a site, the
   %  manual channel is used instead. The number of soil-temperature
   %  columns equals the number of `sdepth` levels in the obs file.
   %
   %  Inputs
   %    sitename : string  ESM-SnowMIP site code (e.g. "cdp").
   %
   %  Name-value
   %    source_dir : string (default data/verification/snow/esm_snowmip)
   %    startdate  : datetime or "" (default "")
   %    enddate    : datetime or "" (default "")
   %
   %  Returns
   %    observations : timetable (verification target schema)
   %    metadata     : struct with provenance + soil depth metadata
   %
   % See also: icemodel.verification.setup.buildEsmSnowmipForcing,
   %  icemodel.verification.setup.importEsmSnowmip

   arguments
      sitename (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate  = ""
      kwargs.enddate    = ""
   end

   if kwargs.source_dir == ""
      source_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'snow', 'esm_snowmip'));
   else
      source_dir = kwargs.source_dir;
   end

   % Locate the obs file (one per site).
   pattern = sprintf('obs_insitu_%s_*.nc', sitename);
   matches = dir(fullfile(source_dir, pattern));
   if isempty(matches)
      error('icemodel:verification:buildEsmSnowmipObservations:fileNotFound', ...
         'no file matching %s in %s', pattern, source_dir);
   end
   if numel(matches) > 1
      names = strjoin({matches.name}, ', ');
      error('icemodel:verification:buildEsmSnowmipObservations:ambiguousFile', ...
         'multiple files matching %s in %s: %s', pattern, source_dir, names);
   end
   obsfile = fullfile(matches.folder, matches.name);

   % --- Read obs channels ---------------------------------------------
   % ESM-SnowMIP sites are heterogeneous: boreal sites lack
   % snd_auto and report snd_gap_auto / snd_gap1_auto instead;
   % sap lacks tsl entirely; rme / sod lack albs; sap lacks sdepth;
   % sod lacks ts. Use site-aware / optional readers so the rest
   % of the builder is generic.
   obs_time = readNetcdfTime(obsfile, 'time');
   snd_auto = readBestSnowDepth(obsfile);
   snd_man  = optionalObs(obsfile, 'snd_man', numel(obs_time));
   swe_auto = optionalObs(obsfile, 'snw_auto', numel(obs_time));
   swe_man  = optionalObs(obsfile, 'snw_man',  numel(obs_time));
   ts       = optionalObs(obsfile, 'ts',       numel(obs_time));

   % Soil temperature is two-dimensional (sdepth x time). Sites
   % without tsl produce zero soil-temperature columns rather than
   % a placeholder NaN soil layer.
   info = ncinfo(obsfile);
   names = string({info.Variables.Name});
   if any(names == "tsl")
      tsl = readObs(obsfile, 'tsl');
   else
      tsl = zeros(0, numel(obs_time));
   end
   if any(names == "sdepth")
      sdepth = double(ncread(obsfile, 'sdepth'));
   else
      sdepth = zeros(0, 1);
   end

   % Merge primary / fallback observations.
   snow_depth_m = preferPrimary(snd_auto, snd_man);
   swe_kg_m2    = preferPrimary(swe_auto, swe_man);

   % --- Optional time-window subset -----------------------------------
   idx = trueMask(numel(obs_time));
   if ~strcmp(string(kwargs.startdate), "")
      idx = idx & obs_time >= ensureUtc(kwargs.startdate);
   end
   if ~strcmp(string(kwargs.enddate), "")
      idx = idx & obs_time <= ensureUtc(kwargs.enddate);
   end
   if ~any(idx)
      error('icemodel:verification:buildEsmSnowmipObservations:emptyWindow', ...
         'no observation samples in the requested window for site %s', sitename);
   end

   Time = obs_time(idx);

   % Soil-temperature columns are expanded into one named variable per
   % depth so comparison variables can pick the right column without
   % introspecting an ND array.
   tsl_window = tsl(:, idx);  % [n_soil x n_time]
   n_soil = size(tsl_window, 1);
   soil_names = "soil_temp_" + string(1:n_soil)' + "_C";
   variable_names = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"; soil_names];
   variable_values = [{snow_depth_m(idx), swe_kg_m2(idx), ts(idx)}, ...
      num2cell(tsl_window', 1)];

   observations = timetable(Time, variable_values{:}, ...
      'VariableNames', cellstr(variable_names));
   observations.Time.TimeZone = 'UTC';

   metadata = struct( ...
      'sitename',         sitename, ...
      'obs_file',         string(obsfile), ...
      'n_rows',           height(observations), ...
      'window_start',     observations.Time(1), ...
      'window_end',       observations.Time(end), ...
      'soil_depths_m',    sdepth(:), ...
      'snow_depth_source', "snd_auto with snd_man fallback", ...
      'swe_source',       "snw_auto with snw_man fallback when available");
end

% =====================================================================
% Local helpers (mirror buildEsmSnowmipForcing for symmetry).
% =====================================================================

function time = readNetcdfTime(pathname, varname)
   raw = double(ncread(pathname, varname));
   units = string(ncreadatt(pathname, varname, 'units'));
   tref = datetime(extractAfter(units, 'hours since '), ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss.S', 'TimeZone', 'UTC');
   time = tref + hours(raw);
   time = time(:);
end

function values = readObs(pathname, varname)
   values = double(ncread(pathname, varname));
   values(values <= -900) = NaN;
end

function values = optionalObs(pathname, varname, ntime)
   info = ncinfo(pathname);
   if any(string({info.Variables.Name}) == varname)
      values = readObs(pathname, varname);
   else
      values = nan(ntime, 1);
   end
end

function values = readBestSnowDepth(pathname)
   %READBESTSNOWDEPTH Site-aware snow-depth selector (mirrors the
   %  forcing builder so both produce the same column).
   info = ncinfo(pathname);
   names = string({info.Variables.Name});
   candidates = ["snd_auto", "snd_gap_auto", "snd_gap1_auto", "snd_man"];
   for c = candidates
      if any(names == c)
         values = readObs(pathname, char(c));
         return
      end
   end
   error('icemodel:verification:buildEsmSnowmipObservations:noSnowDepth', ...
      'no usable snow-depth channel in %s', pathname);
end

function merged = preferPrimary(primary, fallback)
   merged = primary;
   replace = ~isfinite(merged) & isfinite(fallback);
   merged(replace) = fallback(replace);
end

function mask = trueMask(n)
   mask = true(n, 1);
end

function t = ensureUtc(t)
   if isstring(t) || ischar(t)
      t = datetime(t, 'TimeZone', 'UTC');
   elseif isdatetime(t) && isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
end
