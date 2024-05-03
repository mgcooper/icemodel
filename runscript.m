clearvars
close all
clc

% This script demonstrates how to run icemodel (or skinmodel) at a point,
% or for an individual catchment with catchment-scale average forcings

saveflag = false;
sitename = 'behar';
forcings = 'kanm';
userdata = 'kanm';
uservars = 'albedo';
smbmodel = 'icemodel';
simyears = 2016:2016;
gridcell = [];

% To run an individual grid cell, specify sitename = 'sector' and the grid ID:
% gridcell = [];

% Run the model
[ice1, ice2, met, opts] = icemodel.run.point(...
   saveflag=saveflag, ...
   sitename=sitename, ...
   forcings=forcings, ...
   smbmodel="skinmodel", ...
   userdata=userdata, ...
   uservars=uservars, ...
   simyears=simyears, ...
   gridcell=gridcell);


% % These were the settings I had to rerun
% sitename = 'sector';
% forcings = 'mar';
% userdata = 'modis';
% simyears = 2013:2013;
% gridcell = 1794;

% % mar albedo:
% forcings = 'mar';
% userdata = 'mar';
% uservars = 'albedo';
% smbmodel = 'icemodel';
% simyears = 2018:2018;
% gridcell = 239; % I think this was upperbasin
% gridcell = 2195; % see notes at end about this point


%% rungrid settings

% This demonstrates how to run icemodel (or skinmodel) for a set of grid cells,
% using a list of grid cell IDs.

testname = 'zobs2';
% Option to run specific grid points:
gridnums = 2195;
% gridnums = [368, 2337]; % bad mar pixels
% gridnums = [2320]; % bad modis pixel

icemodel.run.grid(gridnums=gridnums, testname=testname)


%% Example of reading from saved nc output to compare with this run

% This compares the missing skinmodel/mar/zobs point 2195 for year 2018 using
% the 'sector' plus 'gridcell' option in runscript instead of the rungrid
% option. The ice1 files were already deleted for those runs so the nc file was
% the only way to compare so this reads a single point from them.
%
% I ran skinmodel/modis and it was identical to the saved data in the nc file.
% Note that this run did not spinup, but it was still identical. Since mar is
% the missing one, I did not compare, but I think the point was to see if there
% was something wrong with that point in particular.

filename = icemodel.netcdf.setfilename( ...
   'ice1', smbmodel, forcings, userdata, sitename, simyears, ...
   fullfile(opts.pathoutput, 'zobs'));

varnames = {'Tsfc', 'runoff'};
data = icemodel.netcdf.ncread(filename, varnames, [gridcell 1], [1 8760]);

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

