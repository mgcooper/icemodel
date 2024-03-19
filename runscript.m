clearvars
close all
clc

%% set the run-specific model configuration
savedata = false;
sitename = 'sector';        % options: 'kanm', 'behar'
forcings = 'mar';         % options: 'mar','kanm'
userdata = 'mar';       % options: 'modis','racmo','merra','mar','kanm','none'
uservars = 'albedo';       % options: 'albedo', or any var in met
simmodel = 'icemodel';    % options: 'icemodel','skinmodel'
simyears = 2018:2018;

% If sitename is 'sector', specify which grid point to run
runpoint = 2195;

% % These were the settings I had to rerun 
% sitename = 'sector';        % options: 'kanm', 'behar'
% forcings = 'mar';         % options: 'mar','kanm'
% userdata = 'modis';          % options: 'modis','racmo','merra','mar','kanm','none'
% simyears = 2013:2013;
% runpoint = 1794;

%% Set the model options
opts = icemodel.setopts(simmodel, sitename, simyears, forcings, ...
   userdata, uservars, savedata);

if ~isempty(runpoint)
   opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
      ['met_' int2str(runpoint) '.mat'])};
end

% run the model
switch simmodel
   case 'icemodel'
      tic; [ice1, ice2, numfail] = icemodel(opts); toc
   case 'skinmodel'
      tic; [ice1, ice2] = skinmodel(opts); toc
end

% load the met data and run the post processing
if opts.savedata
   [ice1, ice2, met] = icemodel.loadresults(opts);
else
   [ice1, ice2, met] = POSTPROC(ice1, ice2, opts, simyears);
end


%% prep the output for plotting
setzero = false;
[Runoff, Discharge, Catchment, Melt] = prepRunoff(opts, ice1, 'set_negative_runoff_zero', setzero);
AblationHourly = prepAblation(opts, ice1, 'hourly', setzero);
AblationDaily = prepAblation(opts, ice1, 'daily', setzero);

% this controls the time period over which ablation and runoff are plotted
t1 = datetime(simyears(1),6,1,0,0,0,'TimeZone','UTC');
t2 = datetime(simyears(1),9,1,0,0,0,'TimeZone','UTC');

%% TEST
melt = Melt.icemodel;
roff = Runoff.icemodel;
roff(roff<0) = 0;

melt0 = cumsum(melt, 'omitnan');
roff1 = cumsum(Runoff.icemodel, 'omitnan');
roff2 = cumsum(roff, 'omitnan');

% Plot before reassigning
figure;
plot(roff1); hold on
plot(roff2, ':')
plot(melt0); 
legend('roff', 'roff2', 'melt')

% compute the percent diff cumulative
(roff2(end) - melt0(end)) / melt0(end)
(roff1(end) - melt0(end)) / melt0(end)

% Reassign runoff
Runoff.icemodel = roff;

%% Plot runoff
if opts.simmodel == "skinmodel"
   [h1, data] = plotRunoff(Runoff, Discharge, Catchment, ...
      'plotsurf', true, 'sitename', sitename, 'userdata', userdata, ...
      'forcingdata', forcings);
else
   [h1, data] = plotRunoff(Runoff, Discharge, Catchment, ...
      'sitename', sitename, 'userdata', userdata, 'forcingdata', forcings, ...
      't1', t1, 't2', t2, 'refstart', false);
end

% plot ablation
h2 = plotPromice(AblationDaily,'refstart',t1);
h3 = plotPromice(AblationHourly,'refstart',t1);

% Repeat with July
t1 = datetime(simyears(1),7,1,0,0,0,'TimeZone','UTC');
h2 = plotPromice(AblationDaily,'refstart',t1);
h3 = plotPromice(AblationHourly,'refstart',t1);

%%

% The standard method above uses prepRunoff which uses prepice1 which loads the
% Data files in getenv('ICEMODELINPUTPATH'). Those don't exsit for racmo
% surface, so they are subsurface. NO CATCHMENT scaling is done there, that's
% done in plotRunoff -> prep_runoff, but there, catchment min/med/max is used
% 
% Below, load_icemodel loads the runoff_ files in getenv('ICEMODELDATAPATH')
% 
L = getlegend;
LString = [L.String, {'Surf', 'Subsurf'}];

% if sitename == "slv1"
%    load('/Users/mattcooper/myprojects/matlab/runoff/data/icemodel/eval/racmo_runoff_slv1_2015.mat')
% elseif sitename == "slv2"
%    load('/Users/mattcooper/myprojects/matlab/runoff/data/icemodel/eval/racmo_runoff_slv2_2015.mat')
% end


[~, RacmoSrf] = loadRunoff(sitename, {'mar', 'racmo'}, data.Time, ...
   'racmo', 'surface');
[~, RacmoSub] = loadRunoff(sitename, {'mar', 'racmo'}, data.Time, ...
   'racmo', 'subsurface');
RsrfL = catchmentRunoff(RacmoSrf.min, Catchment.min.ease.area, 'kg/m2/s', false, true);
RsrfH = catchmentRunoff(RacmoSrf.max, Catchment.max.ease.area, 'kg/m2/s', false, true);
RsubL = catchmentRunoff(RacmoSub.min, Catchment.min.ease.area, 'kg/m2/s', false, true);
RsubH = catchmentRunoff(RacmoSub.max, Catchment.max.ease.area, 'kg/m2/s', false, true);


hold on
plot(datenum(data.Time), 1e9 * (RsrfL + RsrfH) / 2, 'b:')
plot(datenum(data.Time), 1e9 * (RsubL + RsubH) / 2, 'b--')
legend(LString, 'location', 'nw')


% RracM = catchment_runoff(Racmo.med, Catchment.med.ease.area,'kg/m2/s',inan,setzero);


% A = (Catchment.min.ease.area + Catchment.max.ease.area) / 2;
% RracL = catchment_runoff(RacmoSrf.min, A, 'kg/m2/s', false, true);
% RracH = catchment_runoff(RacmoSrf.max, A, 'kg/m2/s', false, true);

%%
% % plot the energy balance
% plotEnbal(ice1, met);
% 
% diffs = (data{end, 2:end} - data.ADCP(end)) ./ data.ADCP(end);
% diffs = array2table(diffs, 'VariableNames', {'Racmo', 'Mar', 'Merra', 'IceModel'})
% [ice1.runoff(end) ice1.melt(end)]
% 
% figure
% plot(ice1.balance)
% ylabel('balance')
% 
% figure
% plot(ice1.chf)
% ylabel('Q_c')
% 
% figure
% histogram(ice1.Tice_numiter, 'Normalization', 'probability')
% ylabel('ICEENBAL numiter')
% 
% fprintf('number successful Tsfc: %.0f\n', sum(ice1.Tsfc_converged))
% fprintf('percent icenbal < 5 iterations: %.2f\n', ...
%    100 * sum(ice1.Tice_numiter < 5) / numel(ice1.Tice_numiter))
% 
% 
% % figure;
% % plot(ice1.Time, ice1.melt); hold on
% % plot(met.Time, cumsum(met.melt))
% % plot(ice1.Time, ice1.runoff);
% % legend('total melt', 'met melt', 'runoff')
% 
% % zD = ice1.zD;
% % [min(zD(zD > 0)) mean(zD(zD > 0)) max(zD(zD > 0))]
% % 
% % figure; plot(ice1.Time, zD); hold on
% % horzline(opts.dz_thermal)
% 
% 
% %%
% % icemodel.diags.stability(met.Time, met.tair, ice1.tsfc, met.wspd, ...
% %       opts.z_0, opts.z_obs, opts.z_wind)
% 
% %%
% 
% % [Lv, ro_liq] = icemodel.physicalConstant('Lv', 'ro_liq');
% % evap = ice1.lhf ./ (Lv * ro_liq) .* opts.dt ./ opts.dz_thermal;
% %
% % figure;
% % plot(ice1.Time, ice1.evap); hold on;
% % plot(ice1.Time, cumsum(evap, 'omitnan'));
% 
% % postive lhf means energy into the surface
