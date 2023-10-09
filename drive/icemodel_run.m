clearvars
close
clc

% this demonstrates how to set up a one-year run at one site

%% set the run-specific model configuration

savedata = true;
sitename = 'behar';        % options: 'kanm', 'behar'
forcings = 'mar';          % options: 'mar','kanm'
userdata = 'modis';        % options: 'modis','racmo','merra','mar','kanm','none'
uservars = 'albedo';       % options: 'albedo', or any var in met
simmodel = 'icemodel';     % options: 'icemodel','skinmodel'
simyears = 2016;

%% set the input and output paths in icemodel_config.m

pathadd(getprojectfolder('runoff'))

cfg = icemodel_config();

%% build the 'opts' model configuration structure

opts = icemodel_opts(sitename,simmodel,simyears,forcings,userdata,uservars,  ...
   savedata);

%% run the model
switch simmodel
   case 'icemodel'
      tic; icemodel(opts); toc
   case 'skinmodel'
      tic; skinmodel(opts); toc
end

% load the met data and run the post processing 
[ice1, ice2, met] = icemodel.loadresults(opts);

%% prep the output for plotting

[Runoff,Discharge,Catchment] = prepRunoff(opts, ice1);
AblationHourly = prepAblation(opts, ice1, 'hourly');
AblationDaily = prepAblation(opts, ice1, 'daily');

% this controls the time period over which ablation and runoff are plotted
t1 = datetime(simyears(1),7,1,0,0,0,'TimeZone','UTC');
t2 = datetime(simyears(1),9,1,0,0,0,'TimeZone','UTC');

% plot the runoff
if opts.simmodel == "skinmodel"
   h1 = plotRunoff(Runoff,Discharge,Catchment,'plotsurf',true,'sitename',  ...
      sitename,'userdata',userdata,'forcingdata',forcings);
else
   h1 = plotRunoff(Runoff,Discharge,Catchment,'sitename',sitename,         ...
      'userdata',userdata,'forcingdata',forcings,'t1',t1,'t2',t2);
end

% plot ablation
t1 = datetime(simyears(1),7,1,0,0,0,'TimeZone','UTC');
plotPromice(AblationDaily,'refstart',t1);
plotPromice(AblationHourly,'refstart',t1);

% plot the energy balance
% plotEnbal(ice1,met);

%%

% [Lv, ro_liq] = icemodel.physicalConstant('Lv', 'ro_liq');
% evap = ice1.lhf ./ (Lv * ro_liq) .* opts.dt ./ opts.dz_thermal;
%
% figure;
% plot(ice1.Time, ice1.evap); hold on;
% plot(ice1.Time, cumsum(evap, 'omitnan'));

% postive lhf means energy into the surface
