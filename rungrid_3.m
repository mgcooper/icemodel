clean

savedata = true;
sitename = 'sector';
forcings = 'mar';
userdata = 'modis';
uservars = 'albedo';
simmodel = 'icemodel';
simyears = 2008:2018;
% gridnums = 671:744;

% 298, 671, 1046,  1420
% 372, 744, 1116, 1487

pathdata = icemodel.setpath('output',sitename,simmodel,userdata);

%%

load('IceMask.mat')
M = IceMask.IceMaskSLA;
E = IceMask.Elev(M);
X = IceMask.X(M);
Y = IceMask.Y(M);
filename = fullfile(getenv('ICEMODELPATH'), 'testbed','data','zobs_vars');
load(filename, 'wspd', 'tair', 'relh')
idx = wspd > 6.5 | tair > 274 | relh < 70;
gridnums = find(idx);

% To run all the ones that weren't rerun
gridnums = find(~idx);

% figure
% scatter(X(gridnums), Y(gridnums), 40, E(gridnums), 'filled')

%% config

% initialize the model options
opts = icemodel.setopts(simmodel, sitename, simyears, forcings, ...
   userdata, uservars, savedata, [], 'zobs');

%
% N = 3;
% si = 1;
% ni = 154;
%
% S = si+ni*(N-1);
% E = si+ni*N-1;

% S = 1;
S = 105;
E = floor(numel(gridnums) / 2);

% display the run information
disp([simmodel ' ' userdata ' ' int2str(simyears(1)) ':'   ...
   int2str(simyears(end)) ' ' int2str(S) ':' int2str(E)])

for n = S:E

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

