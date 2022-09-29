function [met,opts] = METINIT(opts)
   
%------------------------------------------------------------------------------
% this section is for gridded (sector-scale) simulations
%------------------------------------------------------------------------------
if opts.sitename == "sector"

   % load the met file
   load(opts.metfname,'met')

   % remove leap inds if the met data is on a leap-year calendar
   if opts.calendar_type == "noleap"
      met = rmleapinds(met);
   end
   
   % compute the total number of model timesteps
   opts.maxiter = height(met)/opts.numyears;
   opts.dt  = seconds(met.Time(2)-met.Time(1));

   % swap out chosen forcing data albedo for modis albedo
   if string(opts.userdata) == "modis"
      met.albedo = met.modis;
   end
   
   return;
   % for a sector-scale run, we are finished
end

%------------------------------------------------------------------------------
% this section is for point-scale or catchment-scale simulations
%------------------------------------------------------------------------------

% point- or catchment-scale forcing file
   load(opts.metfname,'met')

% next swaps out a variable from userdata in place of the default met data
if opts.userdata ~= "none"
    
   % for point/catchment-scale runs, multi-year simulations are not presently
   % supported, so the simyear should equal the start/endyear
   simyear     =   num2str(opts.simyears(1));
   sitename    =   opts.sitename;
   userdata    =   opts.userdata;
   uservars    =   opts.uservars;
   userpath    =   opts.userpath;
   
   % convert uservars to a cellstr for compatibility with multiple uservars
   if ischar(uservars)
      uservars = cellstr(uservars);
   end
   
   numuservars =   numel(uservars);
    
% load the forcingUserData file and retime to match the metfile
   userfile    =   [userpath userdata '_' sitename '_' simyear '.mat'];
   
   if ~exist(userfile,'file')
      error('userdata file does not exist');
   end
   
   % the userdata is hourly, retime to 15 m if the met data is 15 m
   load(userfile,'Data');
   if Data.Time(2)-Data.Time(1) ~= met.Time(2)-met.Time(1)
      Data  = retime(Data,met.Time,'linear');
   end

   % swap out the data in the metfile for the user data
   for n = 1:numuservars

      thisvar  = uservars{n};
      swapvar  = thisvar;

      % MODIS requires changing the variable name to albedo
      if lower(userdata) == "modis" && lower(thisvar) == "albedo"
         thisvar = 'albedo';
         swapvar = 'modis';
      end

      met.(thisvar) = Data.(swapvar);
   end
end

if isempty(met.Time.TimeZone)
   met.Time.TimeZone = 'UTC';
end

%%%% test
% met = retime(met,'hourly','mean');
%%%% test

% reset the timestep using the built-in time functions to ensure precision
opts.maxiter   = size(met,1);
opts.dt        = seconds(met.Time(2)-met.Time(1));

