clean

% There was a version of this script called _test that tested only writing the
% 2d vars (and maybe only the 0-d vars based on the note below but the note may
% be a mistake, it seems it tested only writing 2-d vars). In the test version,
% comprsz was set to 9, so I added that to the if test_compression == true
% section, but I think that was just a one-off setting, separate from the actual
% test of including 0-d, 1-d, or 2-d vars.

savedata = true;
sitename = 'sector';
simyears = 2008:2018;
simmodel = 'icemodel'; % {'icemodel', 'skinmodel'};
userdata = 'mar'; % {'mar', 'modis'};
siteopts = setBasinOpts('sitename', sitename, 'simmodel', simmodel, 'userdata', userdata);

make_backups = false;

% % I also tested if compression is better w/o the 1-d and 2-d vars - it isn't
% test_compression = true;
% if test_compression == true
%    deflateLevel = 9;           % compression size
% else
%    deflateLevel = 1;           % compression size
% end
% file size in mb/cell for min/med/max compression
% 1: 10.9 (10.8 w/o 1-d or 0-d vars)
% 5: 9.4
% 9: 9.2 (9.1 w/o 1-d or 0-d vars)


% Set path to data
pathdata = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename, simmodel, userdata);
pathsave = pathdata;

%%

% [info, data] = icemodel.netcdf.makencfile(pathdata, pathsave, simyears);
% icemodel.netcdf.makencfile(pathdata, pathsave, simyears(2:end), ...
%    deflateLevel=9, test_write=false);

icemodel.netcdf.makencfile(pathdata, pathsave, 2012, ...
   deflateLevel=9, test_write=true);



%%
info = ncinfo(fullfile(pathdata, 'icemodel_2008.nc4'));

tic
f_ice = ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'f_ice');
toc

numhrs = 8784;
numgrids = 1487;

% Takes about 0.92 seconds to load the data
tic
data = load(fullfile(fileparts(fileparts(pathdata)), 'icemodel_mar_2008.mat'));
toc

% Takes about 0.88 seconds to create the equivalent data from the nc4 files
tic
data.M = nan(numhrs, numgrids);
data.F = nan(numhrs, numgrids);
data.R = nan(numhrs, numgrids);
data.Tsfc = nan(numhrs, numgrids);

data.M = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'depth_melt'));
data.F = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'depth_freeze'));
data.R = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'runoff'));
data.Tsfc = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'Tsfc'));
toc


%    test = ncread(outfilename,'f_liq');
%    test1 = squeeze(test(1,:,:));
%    test2 = squeeze(test(2,:,:));
%    test3 = squeeze(test(3,:,:));
%
%    figure; plot(test1(1,:)); hold on; plot(test1(2,:)); plot(test1(3,:));


% % % % % % % % % % % % % % % % % % % % % % % % % %
% % % this shows that runoff cannot be computed from f_liq after the fact,
% % whereas it can be computed from df_liq, so I need to save df_liq. I
% % should be able to compute ro_sno, cp_sno, etc from f_liq and f_ice, so I
% % still need to save them along with Tice
%
%
% vzero = zeros(size(ice2.f_liq(:,1)));
% df_test = [vzero diff(ice2.f_liq,1,2)];
%
% dz = 0.04;
% runoff_test = tocolumn(cumsum( sum( 4*dz.*df_test ) ));
% % runoff_test = tocolumn(cumsum( sum( dz.*ice2.df_liq ) ));
%
% figure;
% plot(ice1.Time,ice1.runoff); hold on;
% plot(ice1.Time,4.*runoff_test,':');
