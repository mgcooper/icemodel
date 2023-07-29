clearvars
close
clc

% this demonstrates how to set up a one-year run at one site

%% set the run-specific model configuration

savedata = false;
sitename = 'behar';        % options: 'kanm', 'behar'
forcings = 'mar';          % options: 'mar','kanm'
userdata = 'modis';        % options: 'modis','racmo','merra','mar','kanm','none'
uservars = 'albedo';       % options: 'albedo', or any var in met
simmodel = 'icemodel';     % options: 'icemodel','skinmodel'
simyears = 2016;

%% source the input and output paths set in icemodel_config.m

icemodel_config()

%% build the 'opts' model configuration structure

opts = icemodel_opts(sitename,simmodel,simyears,forcings,userdata,uservars,  ...
   savedata);

%% run the model

tic
switch simmodel
   case 'icemodel'
      [ice1,ice2,met,opts] = icemodel(opts);
   case 'skinmodel'
      [ice1,ice2,met,opts] = skinmodel(opts);
end
toc

%% save the data

if savedata == true
   if ~isfolder(opts.pathoutput); mkdir(opts.pathoutput); end
   save([opts.pathoutput opts.fsave],'ice1','ice2','opts');
end

%% prep the output for plotting

[Runoff,Discharge,Catchment]  = prepRunoff(opts, ice1);
AblationHourly                = prepAblation(opts, ice1, 'hourly');
AblationDaily                 = prepAblation(opts, ice1, 'daily');

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


