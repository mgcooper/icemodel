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
   % userdata swaps remain aligned with the source met file year.
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
   if shouldSwapUserdata(opts)
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
               Data = loadExternalUserdata(opts, thisyear, met.Time(ii));
            end
            sourcevar = findSourceVar(Data.Properties.VariableNames, ...
               externalSourceCandidates(opts.userdata, targetvar));
            if isempty(sourcevar)
               error('userdata variable for "%s" not found in %s source', ...
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
function tf = shouldSwapUserdata(opts)

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
function Data = loadExternalUserdata(opts, thisyear, mettime)

   filepath = resolveUserdataFile(opts, thisyear, mettime);
   if ~isfile(filepath)
      error('\n userdata file does not exist: \n\n %s \n', filepath);
   end

   Data = load(filepath, 'Data');
   if ~isfield(Data, 'Data')
      error('userdata file does not contain timetable "Data": %s', filepath);
   end
   Data = Data.Data;
   Data.Time.TimeZone = mettime.TimeZone;
   Data = retime(Data, mettime, 'linear');
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

   base = [opts.sitename '_' char(opts.userdata)];
   enclosing = icemodel.forcing.helpers.findEnclosingWindowFile( ...
      opts.pathuserdata, base, '.mat', min(mettime), max(mettime));
   if strlength(enclosing) > 0
      filepath = fullfile(opts.pathuserdata, char(enclosing));
   else
      filepath = fullfile(opts.pathuserdata, ...
         [base '_' int2str(thisyear) '.mat']);
   end
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
