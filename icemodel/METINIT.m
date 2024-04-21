function [tair, swd, lwd, albedo, wspd, rh, psfc, rain, tppt, ...
      De, S, time] = METINIT(opts, ii)
   %METINIT initialize the met file
   %
   %#codegen

   % The 2nd input is the index into the metfile name list generated in setopts
   if nargin < 2
      ii = 1;
   end

   % Load the met file
   met = load(opts.metfname{ii}, 'met');
   met = met.met;

   if strcmp('none', opts.userdata) == false
      met = swapMetData(met, opts);
   end
   met = prepareMetData(met, opts);

   % Transfer the met data to vectors
   rh = met.rh;
   swd = met.swd;
   lwd = met.lwd;
   tair = met.tair;
   wspd = met.wspd;
   psfc = met.psfc;
   time = met.Time;
   albedo = met.albedo;

   if isvariable('rain', met)
      rain = met.rain / 3600; % opts.dt - the mar rain/snow is mWE / hr not / dt
      rain(isnan(rain)) = 0.0;
      rain = 0 * tair;
   else
      rain = 0 * tair;
   end

   % Solve for wet bulb
   [Ls, cp_air] = icemodel.physicalConstant('Ls', 'cp_air');
   tppt = nan(size(rh));
   for n = 1:numel(rh)
      tppt(n) = SOLVEWB(tair(n), rh(n), Ls, cp_air, psfc(n));
   end

   % Compute the wind transfer and stability coefficients
   [De, S] = WINDCOEF(wspd, opts.z_0, opts.z_obs, opts.z_wind);
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

   if strcmp(opts.sitename, 'sector')

      % Swap out the forcing data albedo for modis albedo

      % Activate this and move the swap loop below outside the else to swap
      % generic variables for gridded (sector) runs.
      % Data = met;
      % if strcmp('modis', opts.userdata)
      %    Data.modis = met.modis;
      % end

      if strcmp('modis', opts.userdata)
         met.albedo = met.modis;
      end

   else
      % Load the userdata and retime from hourly to 15 m if the met data is 15 m
      simyears = num2str(opts.simyears(1));
      userfile = [opts.userdata '_' opts.sitename '_' simyears '.mat'];

      if isfile(fullfile(opts.pathinput, 'userdata', userfile))
         Data = load(fullfile(opts.pathinput, 'userdata', userfile), 'Data');
         Data = Data.Data;
         Data.Time.TimeZone = met.Time.TimeZone;
         if Data.Time(2)-Data.Time(1) ~= met.Time(2)-met.Time(1)
            Data = retime(Data, met.Time, 'linear');
         end
      else
         error('\n userdata file does not exist: \n\n %s \n', ...
            fullfile(opts.pathinput, 'userdata', userfile));
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
         met.(thisvar) = Data.(swapvar);
      end
   end
end
