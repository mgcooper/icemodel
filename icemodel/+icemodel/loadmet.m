function [met, opts] = loadmet(opts, fileiter) %#codegen
   %LOADMET Load one or more icemodel met files as a timetable.
   %
   %  met = icemodel.loadmet(opts) loads and concatenates all met files named
   %  in opts.metfname. This is the canonical path for multi-year runs.
   %
   %  met = icemodel.loadmet(opts, fileiter) loads only the requested met file
   %  index/indices from opts.metfname.

   if nargin < 1
      opts = icemodel.setopts('icemodel', 'behar', 2016, 'kanm');
   end
   if nargin < 2 || isempty(fileiter)
      fileiter = 1:numel(opts.metfname);
   end

   % Load and post-process each requested file before concatenation so yearly
   % swap-source files remain aligned with the source met file year.
   metcell = cell(1, numel(fileiter));
   for n = 1:numel(fileiter)
      metcell{n} = loadOneMetFile(opts, fileiter(n));
   end

   if isscalar(metcell)
      met = metcell{1};
   else
      met = vertcat(metcell{:});
      met = sortrows(met);
   end

   met = addCanonicalSnowDepth(met);

   % Require the concatenated met data to cover every requested simulation
   % year, whether the source is one multi-year file or one file per year.
   if ~all(ismember(opts.simyears(:), unique(year(met.Time))))
      error('met data do not cover all requested simulation years')
   end
   
   met.De ...
      = icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      met.wspd, opts.z0_bulk, opts.z_tair, opts.z_wind);
end

%%
function met = loadOneMetFile(opts, fileiter)

   % Load the met file
   met = load(opts.metfname{fileiter}, 'met');
   met = met.met;

   met = prepareMetData(met, opts);
   if shouldSwapData(opts)
      met = swapMetData(met, opts);
   end
end

%%
function met = prepareMetData(met, opts)
   %PREPAREMETDATA remove leap inds, trim to simyears, check for bad data

   % remove leap inds if the met data is on a leap-year calendar
   if strcmp('noleap', opts.calendar_type)
      feb29 = month(met.Time) == 2 & day(met.Time) == 29;
      met = met(~feb29,:);
   end

   % subset the met file to the requested simyears
   met = met(ismember(year(met.Time), opts.simyears), :);

   met.Time.TimeZone = 'UTC';

   % Optional explicit datetime-window override. When opts.startdate
   % and/or opts.enddate are set, narrow the met data to the
   % requested window in addition to the year-granularity simyears
   % subset above. Either bound may be left as NaT to disable that
   % side of the window.
   if isfield(opts, 'startdate') && ~isnat(opts.startdate)
      met = met(met.Time >= opts.startdate, :);
   end
   if isfield(opts, 'enddate') && ~isnat(opts.enddate)
      met = met(met.Time <= opts.enddate, :);
   end
end

%%
function met = addCanonicalSnowDepth(met)
   %ADDCANONICALSNOWDEPTH Standardize optional forcing snow-depth aliases.

   if isvariable('snow_depth', met)
      snow_depth = met.snow_depth;
   else
      snow_depth = nan(height(met), 1);
   end

   aliases = {'snowd', 'snow'};
   for n = 1:numel(aliases)
      alias = aliases{n};
      if ~isvariable(alias, met)
         continue
      end
      values = met.(alias);
      replace = ~isfinite(snow_depth) & isfinite(values);
      snow_depth(replace) = values(replace);
   end

   if isvariable('snow_depth', met) || any(cellfun(@(v) isvariable(v, met), aliases))
      snow_depth(~isfinite(snow_depth)) = nan;
      met.snow_depth = snow_depth;
   end
end

%%
function met = swapMetData(met, opts)

   uservars = normalizeUservars(opts.uservars);
   simyear = year(met.Time);
   for thisyear = reshape(unique(simyear)', 1, [])
      ii = simyear == thisyear;
      Data = [];

      for n = 1:numel(uservars)
         targetvar = uservars{n};
         requireVariable(met, targetvar, 'met');

         sourcevar = findSourceVar(met.Properties.VariableNames, ...
            inlineSourceCandidates(opts.userdata, targetvar));
         if isempty(sourcevar)
            if isempty(Data)
               Data = loadExternalSwapData(opts, thisyear, met.Time(ii));
            end
            sourcevar = findSourceVar(Data.Properties.VariableNames, ...
               externalSourceCandidates(opts.userdata, targetvar));
            if isempty(sourcevar)
               error('swap variable for "%s" not found in %s source', ...
                  targetvar, opts.userdata);
            end
            swapdata = Data.(sourcevar);
         else
            swapdata = met.(sourcevar)(ii);
         end

         swapdata = sanitizeSwapData(swapdata, sourcevar, met.(targetvar)(ii));
         met.(targetvar)(ii) = swapdata;
      end
   end
end

%%
function tf = shouldSwapData(opts)

   tf = ~(isempty(opts.userdata) || isblanktext(opts.userdata)) ...
      && ~strcmpi(opts.userdata, 'none') ...
      && ~strcmpi(opts.forcings, opts.userdata);
end

%%
function uservars = normalizeUservars(uservars)

   if ischar(uservars)
      uservars = cellstr(uservars);
   elseif isstring(uservars)
      uservars = cellstr(uservars);
   end
end

%%
function Data = loadExternalSwapData(opts, thisyear, mettime)
   %LOADEXTERNALSWAPDATA Load a met-preferred external swap source.
   % The public option is still named USERDATA for backward compatibility, but
   % the selected value is a source label. Prefer met/<source>/ files because
   % forcing-channel swaps are meteorological data; fall back to legacy
   % userdata/<source>/ Data files only when no matching met file is staged.

   [filepath, kind] = resolveSwapSourceFile(opts, thisyear, mettime);
   if strlength(filepath) == 0 || ~isfile(filepath)
      error('\n swap source file does not exist: \n\n %s \n', filepath);
   end

   switch kind
      case "met"
         payload = load(filepath, 'met');
         if ~isfield(payload, 'met')
            error('swap met file does not contain timetable "met": %s', ...
               filepath);
         end
         Data = payload.met;

      case "userdata"
         payload = load(filepath, 'Data');
         if ~isfield(payload, 'Data')
            error('userdata file does not contain timetable "Data": %s', ...
               filepath);
         end
         Data = payload.Data;
         rejectDailyUserdata(Data, filepath, opts.userdata);

      otherwise
         error('unrecognized swap source kind: %s', kind);
   end

   if ~istimetable(Data)
      error('swap source file does not contain a timetable: %s', filepath);
   end
   Data.Time.TimeZone = mettime.TimeZone;
   Data = retime(Data, mettime, 'linear');
end

%%
function [filepath, kind] = resolveSwapSourceFile(opts, thisyear, mettime)
   %RESOLVESWAPSOURCEFILE Prefer met files, then legacy userdata files.

   % An explicit manifest-selected Data artifact is authoritative. This is how
   % native 30-minute and default-hourly variants with the same site/source are
   % kept distinct at runtime; legacy calls without the option remain met-first.
   if ~isempty(explicitUserdataFiles(opts))
      filepath = resolveUserdataFile(opts, thisyear, mettime);
      kind = "userdata";
      return
   end

   [filepath, found] = resolveSwapMetFile(opts, thisyear, mettime);
   if found
      kind = "met";
      return
   end

   filepath = resolveUserdataFile(opts, thisyear, mettime);
   kind = "userdata";
end

%%
function [filepath, found] = resolveSwapMetFile(opts, thisyear, mettime)
   %RESOLVESWAPMETFILE Locate a met file for the selected swap source.
   % Runtime met may be 15 min while reusable RCM swap sources are hourly, so
   % search the active run cadence first and then the other supported cadence.

   found = false;
   filepath = "";
   source = char(opts.userdata);
   met_base = fullfile(opts.pathinput, 'met');
   metname = swapMetName(opts.sitename, source);
   base = ['met_' metname '_' source];
   suffixes = swapMetSuffixes(opts);

   for s = 1:numel(suffixes)
      suffix = suffixes{s};
      for d = icemodel.forcing.helpers.sourceSearchDirs(met_base, source)
         met_dir = d{1};
         enclosing = icemodel.forcing.helpers.findEnclosingWindowFile( ...
            met_dir, base, ['_' suffix '.mat'], min(mettime), max(mettime));
         if strlength(enclosing) > 0
            filepath = fullfile(met_dir, char(enclosing));
            found = true;
            return
         end

         peryear_name = [base '_' int2str(thisyear) '_' suffix '.mat'];
         candidate = fullfile(met_dir, peryear_name);
         if isfile(candidate)
            filepath = candidate;
            found = true;
            return
         end
      end
   end
end

%%
function metname = swapMetName(sitename, source)
   %SWAPMETNAME Match the forcing-station aliasing used by met-file loading.

   if strcmpi(source, 'kanl') && ismember(sitename, {'ak4', 'upperbasin'})
      metname = 'kanl';
   elseif strcmpi(source, 'kanm') && ismember(sitename, {'slv1', 'slv2', 'behar'})
      metname = 'kanm';
   else
      metname = char(sitename);
   end
end

%%
function suffixes = swapMetSuffixes(opts)
   %SWAPMETSUFFIXES Prefer run cadence, then established compatible fallbacks.

   switch opts.dt
      case 900
         suffixes = {'15m', '1hr'};
      case 1800
         % Native 30-minute runs prefer exact swap met, then the finer default
         % model-met artifact, then the legacy hourly source.
         suffixes = {'30m', '15m', '1hr'};
      case 3600
         suffixes = {'1hr', '15m'};
      otherwise
         error('unsupported dt for swap met file naming: %g', opts.dt)
   end
end

%%
function rejectDailyUserdata(Data, filepath, source)
   %REJECTDAILYUSERDATA Keep daily observations out of met-channel swapping.

   if ~istimetable(Data) || height(Data) < 2
      return
   end
   dt_seconds = seconds(diff(Data.Time));
   dt_seconds = dt_seconds(isfinite(dt_seconds));
   if isempty(dt_seconds)
      return
   end
   if median(dt_seconds) >= 23 * 3600
      error('icemodel:loadmet:dailySwapData', ...
         ['swap source "%s" resolved only to daily userdata: %s. ' ...
         'Stage a met file or hourly userdata for met-channel swapping.'], ...
         source, filepath);
   end
end

%%
function filepath = resolveUserdataFile(opts, thisyear, mettime)
   %RESOLVEUSERDATAFILE Locate the userdata file covering this run year's met.
   % Prefers a full-period window file <site>_<source>_<YYYYMMDD>_<YYYYMMDD>.mat
   % whose encoded period brackets the met samples being swapped (the
   % writeuserdata naming="window" form); falls back to the legacy per-year
   % <site>_<source>_<YYYY>.mat. The caller retimes whichever file onto METTIME,
   % so a single full-period file serves every run year. The window lookup is the
   % shared icemodel.forcing.helpers.findEnclosingWindowFile (same primitive
   % icemodel.createMetFileNames uses for met files), bracketed by the actual
   % met time span rather than the whole calendar year.

   % Manifest-selected paths take precedence over cadence-blind legacy name
   % discovery. Select the widest explicit artifact that actually brackets this
   % year's met support so multi-year lists remain deterministic.
   explicit = explicitUserdataFiles(opts);
   if ~isempty(explicit)
      filepath = selectExplicitUserdataFile(explicit, mettime);
      return
   end

   % Look in the per-source subfolder userdata/<source>/ first, then the flat
   % userdata dir, in each preferring the window file then the per-year file.
   % The subfolder-first ordering is the shared sourceSearchDirs primitive.
   base = [opts.sitename '_' char(opts.userdata)];
   peryear_name = [base '_' int2str(thisyear) '.mat'];
   for u = icemodel.forcing.helpers.sourceSearchDirs( ...
         opts.pathuserdata, opts.userdata)
      udir = u{1};
      enclosing = icemodel.forcing.helpers.findEnclosingWindowFile( ...
         udir, base, '.mat', min(mettime), max(mettime));
      if strlength(enclosing) > 0
         filepath = fullfile(udir, char(enclosing));
         return
      end
      if isfile(fullfile(udir, peryear_name))
         filepath = fullfile(udir, peryear_name);
         return
      end
   end
   % Not found in either layout: return the flat per-year path so the caller
   % raises a clean "userdata file does not exist" error.
   filepath = fullfile(opts.pathuserdata, peryear_name);
end

function files = explicitUserdataFiles(opts)
   %EXPLICITUSERDATAFILES Normalize optional caller/manifest-selected paths.
   files = strings(0, 1);
   if ~isfield(opts, 'userdatafname') || isempty(opts.userdatafname)
      return
   end
   files = reshape(strtrim(string(opts.userdatafname)), [], 1);
   files = files(strlength(files) > 0);
   for n = 1:numel(files)
      if ~isAbsolutePath(files(n))
         files(n) = fullfile(string(opts.pathuserdata), files(n));
      end
   end
end

function tf = isAbsolutePath(filename)
   %ISABSOLUTEPATH Recognize POSIX and drive-qualified absolute paths.
   tf = startsWith(filename, string(filesep)) ...
      || ~isempty(regexp(char(filename), '^[A-Za-z]:[\\/]', 'once'));
end

function filepath = selectExplicitUserdataFile(files, mettime)
   %SELECTEXPLICITUSERDATAFILE Select an existing file enclosing METTIME.
   exists = isfile(files);
   if ~any(exists)
      % Preserve the normal missing-file diagnostic with the first exact path.
      filepath = files(1);
      return
   end

   candidates = files(exists);
   enclosing = strings(numel(candidates), 1);
   durations = nan(numel(candidates), 1);
   n_enclosing = 0;
   for n = 1:numel(candidates)
      % Explicit manifest paths are authoritative; corrupt referenced files
      % should surface their load error instead of silently selecting a sibling.
      saved = load(candidates(n), 'Data');
      if ~isfield(saved, 'Data') || ~istimetable(saved.Data) ...
            || isempty(saved.Data)
         continue
      end
      candidate_time = saved.Data.Time;
      candidate_time.TimeZone = mettime.TimeZone;
      if min(candidate_time) <= min(mettime) ...
            && max(candidate_time) >= max(mettime)
         n_enclosing = n_enclosing + 1;
         enclosing(n_enclosing) = candidates(n);
         durations(n_enclosing) = seconds( ...
            max(candidate_time) - min(candidate_time));
      end
   end
   if n_enclosing == 0
      error('icemodel:loadmet:explicitUserdataCoverage', ...
         ['Explicit userdata artifacts do not cover the requested met support ' ...
         '%s to %s: %s'], string(min(mettime)), string(max(mettime)), ...
         strjoin(candidates, ', '))
   end

   % Match legacy enclosing-window selection: widest support, then lexical path.
   enclosing = enclosing(1:n_enclosing);
   durations = durations(1:n_enclosing);
   rank = table(-durations, enclosing, ...
      'VariableNames', {'negative_duration', 'filename'});
   [~, order] = sortrows(rank, {'negative_duration', 'filename'});
   filepath = enclosing(order(1));
end

%%
function candidates = inlineSourceCandidates(userdata, targetvar)

   if strcmpi(userdata, 'modis') && strcmpi(targetvar, 'albedo')
      candidates = {'modis', 'MODIS'};
   else
      candidates = {};
   end
end

%%
function candidates = externalSourceCandidates(userdata, targetvar)

   candidates = {targetvar};
   if strcmpi(userdata, 'modis') && strcmpi(targetvar, 'albedo')
      candidates = [{'modis', 'MODIS'}, candidates];
   end
end

%%
function varname = findSourceVar(varnames, candidates)

   varname = '';
   for n = 1:numel(candidates)
      idx = find(strcmpi(varnames, candidates{n}), 1);
      if ~isempty(idx)
         varname = varnames{idx};
         return
      end
   end
end

%%
function requireVariable(T, varname, label)

   if ~isvariable(varname, T)
      error('%s does not contain variable "%s"', label, varname);
   end
end

%%
function values = sanitizeSwapData(values, sourcevar, fallback)

   if any(strcmpi(sourcevar, {'modis'}))
      bad = values <= 0 | values >= 1;
      if any(bad)
         values(bad) = fallback(bad);
         warning('bad albedo')
      end
   end
end
