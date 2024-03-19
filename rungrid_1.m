clean

savedata = true;
sitename = 'sector';
forcings = 'mar';
userdata = 'mar';
uservars = 'albedo';
simmodel = 'icemodel';
simyears = 2008:2018;

gridnums = loadZobsGridNums("round2", makeplot=false);

% Set the path where the output will be saved.
pathdata = icemodel.setpath('output', sitename, simmodel, userdata);

% Divide the gridpoints into jobs for running in chunks.
[~, E] = icemodel.chunkgridcell(1, numel(gridnums), 2, 1);

% Overrule S if restarting.
S = 257;

%% config

% Initialize the model options.
opts = icemodel.setopts(simmodel, sitename, simyears, forcings, ...
   userdata, uservars, savedata, [], 'zobs');

% Display the run information.
disp([simmodel ' ' userdata ' ' int2str(simyears(1)) ':'   ...
   int2str(simyears(end)) ' ' int2str(S) ':' int2str(E)])

% Run the model for each point.
for n = S:E

   thispoint = gridnums(n);

   opts.casename = int2str(thispoint);
   opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
      ['met_' int2str(thispoint) '.mat'])};

   % run the model
   switch simmodel
      case 'icemodel'
         icemodel(opts);
      case 'skinmodel'
         skinmodel(opts);
   end
end
