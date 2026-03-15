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

   if ~strcmp('none', opts.userdata) && ~strcmp(opts.forcings, opts.userdata)
      met = swapMetData(met, opts);
   end
   met = prepareMetData(met, opts);
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

   % check for bad modis data before swapping out default albedo for modis
   varmodis = {'modis', 'MODIS'};
   hasmodis = ismember(varmodis, met.Properties.VariableNames);
   if any(hasmodis)
      bi = met.(varmodis{hasmodis}) <= 0 | met.(varmodis{hasmodis}) >= 1;
      if sum(bi) > 0
         met.(varmodis{hasmodis})(bi) = met.albedo(bi);
         warning('bad albedo')
      end
   end
   met.Time.TimeZone = 'UTC';
end

%%
function met = swapMetData(met, opts)

   % convert uservars to a cellstr for compatibility with multiple uservars
   if ischar(opts.uservars)
      opts.uservars = cellstr(opts.uservars);
   end

   is_grid_run = strcmp(opts.sitename, 'sector');

   if is_grid_run

      % For grid runs, the met file already carries the grid-cell forcing and
      % userdata swaps are currently limited to replacing forcing albedo with
      % the embedded MODIS albedo field.

      % Activate this and move the swap loop below outside the else to swap
      % generic variables for gridded runs.
      % Data = met;
      % if strcmp('modis', opts.userdata)
      %    Data.modis = met.modis;
      % end

      if strcmp('modis', opts.userdata)
         met.albedo = met.modis;
      end

   else

      % Most met files should have a MODIS column
      if strcmp('modis', opts.userdata)
         if isvariable('MODIS', met)
            met.albedo = met.MODIS;

         elseif isvariable('modis', met)
            met.albedo = met.modis;
         end
         return
      end

      % Load the userdata by year and swap it into the matching met rows.
      for thisyear = reshape(unique(year(met.Time))', 1, [])
         userfile = [opts.sitename '_' opts.userdata '_' int2str(thisyear) '.mat'];
         filepath = fullfile(opts.pathuserdata, userfile);
         if exist(filepath, 'file') ~= 2
            error('\n userdata file does not exist: \n\n %s \n', filepath);
         end

         Data = load(filepath, 'Data');
         Data = Data.Data;
         ii = year(met.Time) == thisyear;
         Data.Time.TimeZone = met.Time.TimeZone;
         if Data.Time(2)-Data.Time(1) ~= met.Time(2)-met.Time(1)
            Data = retime(Data, met.Time(ii), 'linear');
         else
            Data = Data(isbetween(Data.Time, met.Time(find(ii, 1)), ...
               met.Time(find(ii, 1, 'last'))), :);
         end

         % Swap out the data in the metfile for the user data
         for n = 1:numel(opts.uservars)
            thisvar  = opts.uservars{n};
            swapvar  = thisvar;

            % MODIS requires changing the variable name to albedo
            if strcmpi('modis', opts.userdata) && strcmpi('albedo', thisvar)
               thisvar = 'albedo';
               swapvar = 'modis';
            end
            met.(thisvar)(ii) = Data.(swapvar);
         end
      end
   end
end
