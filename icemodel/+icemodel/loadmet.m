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

   % Require the concatenated met data to cover every requested simulation
   % year, whether the source is one multi-year file or one file per year.
   if ~all(ismember(opts.simyears(:), unique(year(met.Time))))
      error('met data do not cover all requested simulation years')
   end

   % Compute the wind transfer coefficient for the processed met data.
   met.De = WINDCOEF(met.wspd, opts.z_0, opts.z_tair, opts.z_wind);
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

   userfile = [opts.sitename '_' opts.userdata '_' int2str(thisyear) '.mat'];
   filepath = fullfile(opts.pathuserdata, userfile);
   if exist(filepath, 'file') ~= 2
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
