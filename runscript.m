
% Apr 2024 - converted runscript to run.point
% copied the settings in runscript below.

%% runscript settings
clearvars
close all
clc

savedata = false;
sitename = 'sector';
forcings = 'mar';
userdata = 'mar';
uservars = 'albedo';
simmodel = 'icemodel';
simyears = 2018:2018;

% If sitename is 'sector', specify which grid point to run
% runpoint = 239; % I think this was upperbasin
runpoint = 2195; % see notes at end about this point

% % These were the settings I had to rerun
% sitename = 'sector';        % options: 'kanm', 'behar'
% forcings = 'mar';           % options: 'mar','kanm'
% userdata = 'modis';         % options: 'modis','racmo','merra','mar','kanm','none'
% simyears = 2013:2013;
% runpoint = 1794;

%% rungrid settings

testname = 'zobs2';
% Option to run specific grid points:
gridnums = 2195;
% gridnums = [368, 2337]; % bad mar pixels
% gridnums = [2320]; % bad modis pixel

icemodel.run.grid(gridnums=gridnums, testname=testname)


%% Example of reading from saved nc output to compare with this run

% This compares the missing skinmodel/mar/zobs point 2195 for year 2018 using
% the 'sector' plus 'runpoint' option in runscript instead of the rungrid
% option. The ice1 files were already deleted for those runs so the nc file was
% the only way to compare so this reads a single point from them.
%
% I ran skinmodel/modis and it was identical to the saved data in the nc file.
% Note that this run did not spinup, but it was still identical. Since mar is
% the missing one, I did not compare, but I think the point was to see if there
% was something wrong with that point in particular.

filename = icemodel.netcdf.setfilename( ...
   'ice1', simmodel, forcings, userdata, sitename, simyears, ...
   fullfile(opts.pathoutput, 'zobs'));

varnames = {'Tsfc', 'runoff'};
data = icemodel.netcdf.ncread(filename, varnames, [runpoint 1], [1 8760]);

figure;
tiledlayout(1, 2)

nexttile; hold on
plot(ice1.Time, ice1.Tsfc)
plot(ice1.Time, data.Tsfc, ':')
ylabel('Tsfc');
legend('ice1 saved', 'this run')

nexttile; hold on
plot(ice1.Time, ice1.runoff)
plot(ice1.Time, data.runoff, ':')
ylabel('runoff');
legend('ice1 saved', 'this run')

figure;
tiledlayout(1, 2)

nexttile
plot(ice1.Tsfc, data.Tsfc, 'o')
legend('Tsfc');
xylabel('ice1 saved', 'this run');
formatPlotMarkers
addOnetoOne

nexttile
plot(ice1.runoff, data.runoff, 'o')
legend('runoff');
xylabel('ice1 saved', 'this run')
formatPlotMarkers
addOnetoOne

