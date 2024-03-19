clean

% Run grid points. See second section to run grid points for each catchment.

savedata = true;
sitename = 'sector';
forcings = 'mar';
userdata = 'mar'; % don't use "none" - set equal to forcings
uservars = 'albedo';
simmodel = 'skinmodel';
simyears = 2008:2018;
% gridnums = 671:744;

% 298, 671, 1046,  1420
% 372, 744, 1116, 1487

% gridnums = [368, 2337]; % bad mar pixels
% gridnums = [2320]; % bad modis pixel

% Set the path where the output will be saved.
pathdata = icemodel.setpath('output', sitename, simmodel, userdata, [], testname);

filename = fullfile(getenv('ICEMODELPATH'), 'testbed','data','zobs_vars');
load(filename, 'wspd', 'tair', 'relh')
idx = wspd > 6.5 | tair > 274 | relh < 70;
gridnums = find(idx);

% TEMP - run the ones I didn't run
gridnums = find(~idx & E < 1200);

figure
scatter(X, Y, 40); hold on
scatter(X(gridnums), Y(gridnums), 40, E(gridnums), 'filled')

%% config

% initialize the model options
opts = icemodel.setopts(simmodel, sitename, simyears, forcings, ...
   userdata, uservars, savedata, [], 'zobs');

% display the run information
disp([simmodel ' ' userdata ' ' int2str(simyears(1)) ':'   ...
   int2str(simyears(end)) ' ' int2str(gridnums(1)) ':' int2str(gridnums(end))])

% run the model for each point (thispoint = 1415, thispoint = 1503 rio behar) 
for n = 1:numel(gridnums)
   
   thispoint = gridnums(n);
   
   % if isfile(fullfile(opts.pathoutput, '2018', ['ice1_' num2str(thispoint) '.mat']))
   %    continue
   % end

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

% load the met data and run the post processing 
% [ice1, ice2, met] = icemodel.loadresults(opts);
% 
% % for testing
% opts.sitename = 'behar';

%% Override the full grid points with catchment points

% allsites = {'behar', 'slv1', 'slv2', 'ak4', 'upperbasin'};
% simyears = {2015:2016, 2015, 2015, 2009:2016, 2016};

% allsites = {'slv1', 'slv2', 'ak4', 'upperbasin'};
% simyears = {2015, 2015, 2009:2016, 2016};
% gridnums = getSiteMetIndex(allsites);

% run the model for each point (thispoint = 1415, thispoint = 1503 rio behar) 
for n = 1:numel(allsites)

   thissite = allsites{n};
   testname = thissite;
   thispoint = gridnums.(thissite);

   % initialize the model options
   opts = icemodel.setopts(simmodel, sitename, simyears{n}, forcings, ...
      userdata, uservars, savedata, [], testname);
   
   opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
      ['met_' int2str(thispoint) '.mat'])};

   % run the model
   tic; icemodel(opts); toc
end

