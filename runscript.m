clearvars
close
clc

% I need this to run a single point if requested. Ideally it would run all
% points inside a catchment and conservatively remap but for now, it at least
% needs to run a given point.
runpoint = true;

%% set the run-specific model configuration
savedata = false;
sitename = 'behar';        % options: 'kanm', 'behar'
forcings = 'mar';         % options: 'mar','kanm'
userdata = 'none';        % options: 'modis','racmo','merra','mar','kanm','none'
uservars = 'albedo';       % options: 'albedo', or any var in met
simmodel = 'icemodel';     % options: 'icemodel','skinmodel'
simyears = 2016;

pathadd(getprojectfolder('runoff'))

%% Set the project configuration and model options
cnfg = icemodel.config();
opts = icemodel.setopts(simmodel, sitename, simyears, forcings, ...
   userdata, uservars, savedata);

%% run the model
switch simmodel
   case 'icemodel'
      tic; [ice1, ice2] = icemodel(opts); toc
   case 'skinmodel'
      tic; [ice1, ice2] = skinmodel(opts); toc
end

% load the met data and run the post processing
if opts.savedata
   [ice1, ice2, met] = icemodel.loadresults(opts);
else
   met = icemodel.loadmet(opts, numel(simyears));
   [ice1, ice2] = POSTPROC(ice1, ice2, opts, ...
      met.swd, met.lwd, met.albedo, met.Time);
end

% [mean(itercount) std(itercount) min(itercount) max(itercount) numel(itercount)]
% 
% figure; histogram(itercount)

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
   [h1, data] = plotRunoff(Runoff,Discharge,Catchment,'sitename',sitename,         ...
      'userdata',userdata,'forcingdata',forcings,'t1',t1,'t2',t2);
end

% plot ablation
t1 = datetime(simyears(1),7,1,0,0,0,'TimeZone','UTC');
plotPromice(AblationDaily,'refstart',t1);
plotPromice(AblationHourly,'refstart',t1);

% plot the energy balance
% plotEnbal(ice1,met);

diffs = (data{end, 2:end} - data.ADCP(end)) ./ data.ADCP(end);
diffs = array2table(diffs, 'VariableNames', {'Racmo', 'Mar', 'Merra', 'IceModel'})
[ice1.runoff(end) ice1.melt(end)]

%%

% [Lv, ro_liq] = icemodel.physicalConstant('Lv', 'ro_liq');
% evap = ice1.lhf ./ (Lv * ro_liq) .* opts.dt ./ opts.dz_thermal;
%
% figure;
% plot(ice1.Time, ice1.evap); hold on;
% plot(ice1.Time, cumsum(evap, 'omitnan'));

% postive lhf means energy into the surface
