%% For sector nc files

% TLDR: this is the start of how I will need to read in from the nc4 files, and
% also could be a function used to call from runscript when I run a gridpoint to
% read the grid point and compare with the simulation output in ice1

% This was to compare a sim using runscript which reran the missing
% skinmodel/mar/zobs point 2195 for year 2018 i.e. used the 'sector' plus
% 'runpoint' option in runscript instead of the rungrid option. 

% I wanted to compare the simulation to the skinmodel/mar one, or the
% skinmodel/modis/zobs one, but I already deleted the ice1 files for both of
% those, so I needed a function to load a point from the nc file, which meant I
% needed a function to set the filename.

% i ran skinmodel/modis and compared the one-year run here, which doesn't
% spinup, to the saved data in the nc file and they were essentially identical. 

filename = icemodel.netcdf.setfilename( ...
   'ice1', simmodel, forcings, userdata, sitename, simyears, fullfile(opts.pathoutput, 'zobs'));

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
