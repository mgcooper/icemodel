function [opts, tair, swdn, lwdn, refl, wspd, relh, psfc, De, Time] = METINIT( ...
   opts, fileiter) %#codegen

% The second input is the index into the metfile name list generated in setopts
if nargin < 2
   fileiter = 1;
end

% load the met file
met = load(opts.metfname{fileiter}, 'met');
met = met.met;

if strcmp('sector', opts.sitename)

   [met, opts] = processMetData(met, opts);
else
   % for point/catchment-scale simulations, option to swap out a variable from
   % userdata in place of the default met data
   if ~strcmp('none', opts.userdata)
      met = swapMetData(met , opts);
   end
   [met, opts] = processMetData(met, opts);
end

% Transfer the met data to vectors
swdn = met.swd;
lwdn = met.lwd;
tair = met.tair;
relh = met.rh;
wspd = met.wspd;
psfc = met.psfc;
refl = met.albedo;
Time = met.Time;

% compute the wind speed transfer coefficient
De = WINDCOEF(wspd, opts.z_0, opts.z_obs);

%%
function [met, opts] = processMetData(met, opts)

% remove leap inds if the met data is on a leap-year calendar
if strcmp('noleap', opts.calendar_type)
   feb29 = month(met.Time) == 2 & day(met.Time) == 29;
   met = met(~feb29,:);
end

% subset the met file to the requested simyears
met = met(ismember(year(met.Time), opts.simyears), :);

% compute the total number of model timesteps and the timestep
opts.maxiter = size(met, 1)/opts.numyears;
opts.dt = seconds(met.Time(2)-met.Time(1));

% check for bad modis data before swapping out default albedo for modis
varmodis = {'modis', 'MODIS'};
hasmodis = ismember(varmodis, met.Properties.VariableNames);
if any(hasmodis)
   bi = met.(varmodis{hasmodis}) <= 0 | met.(varmodis{hasmodis}) >= 1;
   if sum(bi) > 0
      met.(varmodis{hasmodis})(bi) = met.albedo(bi);
      warning('bad albedo')
   end

   % % THIS IS WHERE THE SECTOR-RUNS WERE SWAPPED, REACTIVATE IF RECONCILING
   % DOESN'T WORK
   % % swap out chosen forcing data albedo for modis albedo
   % if strcmp('modis', opts.userdata)
   %    met.albedo = met.modis;
   % end
end

met.Time.TimeZone = 'UTC';

%%
function met = swapMetData(met , opts)

% convert uservars to a cellstr for compatibility with multiple uservars
if ischar(opts.uservars)
   opts.uservars = cellstr(opts.uservars);
end

% load the forcingUserData file and retime to match the metfile
simyears = num2str(opts.simyears(1));
userfile = [opts.userdata '_' opts.sitename '_' simyears '.mat'];

if isfile(fullfile(opts.userpath, userfile))
   Data = load(fullfile(opts.userpath, userfile), 'Data');
   Data = Data.Data;
else
   error('userdata file does not exist');
end

% the userdata is hourly, retime to 15 m if the met data is 15 m
if Data.Time(2)-Data.Time(1) ~= met.Time(2)-met.Time(1)
   Data = retime(Data, met.Time, 'linear');
end

% swap out the data in the metfile for the user data
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