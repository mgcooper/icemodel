clearvars
close all
clc

% I need this to run a single point if requested. Ideally it would run all
% points inside a catchment and conservatively remap but for now, it at least
% needs to run a given point.
runpoint = true;

%% set the run-specific model configuration
savedata = false;
sitename = 'upperbasin';        % options: 'kanm', 'behar'
forcings = 'mar';         % options: 'mar','kanm'
userdata = 'modis';          % options: 'modis','racmo','merra','mar','kanm','none'
uservars = 'albedo';       % options: 'albedo', or any var in met
simmodel = 'icemodel';    % options: 'icemodel','skinmodel'
simyears = 2016;

%% Set the model options
opts = icemodel.setopts(simmodel, sitename, simyears, forcings, ...
   userdata, uservars, savedata);

% opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
%       ['met_' int2str(333) '.mat'])};

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

[Runoff,Discharge,Catchment] = prepRunoff(opts, ice1);
AblationHourly = prepAblation(opts, ice1, 'hourly');
AblationDaily = prepAblation(opts, ice1, 'daily');

% this controls the time period over which ablation and runoff are plotted
t1 = datetime(simyears(1),6,1,0,0,0,'TimeZone','UTC');
t2 = datetime(simyears(1),9,1,0,0,0,'TimeZone','UTC');

%% TEST
melt = Runoff.icemodelMelt;
roff = Runoff.icemodelRunoff;
roff(roff<0) = 0;

melt0 = cumsum(melt, 'omitnan');
roff1 = cumsum(Runoff.icemodelRunoff, 'omitnan');
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
Runoff.icemodelRunoff = roff;

%% Plot runoff
if opts.simmodel == "skinmodel"
   [h1, data] = plotRunoff(Runoff, Discharge, Catchment, 'plotsurf', true, ...
      'sitename', sitename, 'userdata', userdata, 'forcingdata', forcings);
else
   [h1, data] = plotRunoff(Runoff, Discharge, Catchment, 'sitename', sitename, ...
      'userdata', userdata, 'forcingdata', forcings, 't1', t1, 't2', t2);
end

%% Plot ablation
t1 = datetime(simyears(1),7,1,0,0,0,'TimeZone','UTC');
h2 = plotPromice(AblationDaily,'refstart',t1);
h3 = plotPromice(AblationHourly,'refstart',t1);

%% Plot the energy balance
plotEnbal(ice1, met);

diffs = (data{end, 2:end} - data.ADCP(end)) ./ data.ADCP(end);
diffs = array2table(diffs, 'VariableNames', {'Racmo', 'Mar', 'Merra', 'IceModel'})
[ice1.runoff(end) ice1.melt(end)]

figure
plot(ice1.balance)
ylabel('balance')

figure
plot(ice1.chf)
ylabel('Q_c')

figure
histogram(ice1.Tice_numiter, 'Normalization', 'probability')
ylabel('ICEENBAL numiter')

fprintf('number successful Tsfc: %.0f\n', sum(ice1.Tsfc_converged))
fprintf('percent icenbal < 5 iterations: %.2f\n', ...
   100 * sum(ice1.Tice_numiter < 5) / numel(ice1.Tice_numiter))


% figure;
% plot(ice1.Time, ice1.melt); hold on
% plot(met.Time, cumsum(met.melt))
% plot(ice1.Time, ice1.runoff);
% legend('total melt', 'met melt', 'runoff')

% zD = ice1.zD;
% [min(zD(zD > 0)) mean(zD(zD > 0)) max(zD(zD > 0))]
%
% figure; plot(ice1.Time, zD); hold on
% horzline(opts.dz_thermal)


%%
% icemodel.diags.stability(met.Time, met.tair, ice1.tsfc, met.wspd, ...
%       opts.z_0, opts.z_obs, opts.z_wind)

%%

% [Lv, ro_liq] = icemodel.physicalConstant('Lv', 'ro_liq');
% evap = ice1.lhf ./ (Lv * ro_liq) .* opts.dt ./ opts.dz_thermal;
%
% figure;
% plot(ice1.Time, ice1.evap); hold on;
% plot(ice1.Time, cumsum(evap, 'omitnan'));

% postive lhf means energy into the surface
