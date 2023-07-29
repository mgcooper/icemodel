function [met,opts] = METINIT(opts, fileiter)

% The second input is the index into the metfile name list generated in setopts
if nargin < 2
   fileiter = 1;
end

% load the met file
load(opts.metfname{fileiter}, 'met')

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

% check for bad albedo data
bi = met.modis <= 0 | met.modis >= 1;
if sum(bi) > 0
   met.modis(bi) = met.albedo(bi);
   warning('bad albedo')
end

% swap out chosen forcing data albedo for modis albedo
if strcmp('modis', opts.userdata)
   met.albedo = met.modis;
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
   load(fullfile(opts.userpath, userfile), 'Data');
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