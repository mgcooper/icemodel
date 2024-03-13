clean

% to finish the makencfile
% - add grid_mapping
% - add a gridcell variable
% - figure out the auxiliary coordinate variable thing

% For the new runs, test on skinmodel, mar, 2008. The raw data is 17.36 GB.

savedata = true;
sitename = 'sector';
simyears = 2009:2018;
simmodel = 'skinmodel'; % {'icemodel', 'skinmodel'};
userdata = 'mar'; % {'mar', 'modis'};
siteopts = setBasinOpts('sitename', sitename, 'simmodel', simmodel, 'userdata', userdata);

make_backups = false;

simyears = 2009:2011;

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

icemodel.netcdf.makencfile(pathdata, pathsave, simyears, deflateLevel=9, test_write=false);

%%
% info_transpose = ncinfo(fullfile(pathdata, 'icemodel_2009.nc4'));
% info_no_transpose = ncinfo(fullfile(pathdata, 'icemodel_2009_no_transpose.nc4'));
%
% Tice_transpose = ncread(fullfile(pathdata, 'icemodel_2009.nc4'), 'Tice');
% Tice_no_transpose = ncread(fullfile(pathdata, 'icemodel_2009_no_transpose.nc4'), 'Tice');
%
% isequal(Tice_transpose, Tice_no_transpose)
%
% onecell_transpose = squeeze(Tice_transpose(1, :, :));
% onecell_no_transpose = squeeze(Tice_no_transpose(1, :, :));
%
% figure; hold on
% plot(mean(onecell_transpose, 2), 1:500)
% plot(mean(onecell_no_transpose, 2), 1:500)
% set(gca, 'YDir', 'reverse')
% legend('transpose', 'no transpose')
%
% % now permute as one would
% % onecell_transpose = flipud(permute(onecell_transpose, [2 1 3]));
% % onecell_no_transpose = flipud(permute(onecell_no_transpose, [2 1 3]));
% onecell_transpose = flipud(permute(onecell_transpose, [2 1 3]));
% onecell_no_transpose = flipud(permute(onecell_no_transpose, [2 1 3]));
%
% figure; hold on
% plot(mean(onecell_transpose, 1), 1:500)
% plot(mean(onecell_no_transpose, 1), 1:500)
% set(gca, 'YDir', 'reverse')
% legend('transpose', 'no transpose')
%
%
%
% % tic
% % f_ice = ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'f_ice');
% % toc
% %
% % numhrs = 8784;
% % numgrids = 1487;
% %
% % % Takes about 0.92 seconds to load the data
% % tic
% % data = load(fullfile(fileparts(fileparts(pathdata)), 'icemodel_mar_2008.mat'));
% % toc
% %
% % % Takes about 0.88 seconds to create the equivalent data from the nc4 files
% % tic
% % data.M = nan(numhrs, numgrids);
% % data.F = nan(numhrs, numgrids);
% % data.R = nan(numhrs, numgrids);
% % data.Tsfc = nan(numhrs, numgrids);
% %
% % data.M = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'depth_melt'));
% % data.F = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'depth_freeze'));
% % data.R = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'runoff'));
% % data.Tsfc = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'Tsfc'));
% % toc
% %
% % %%
% %
% %
% % % skinmodel is not sensitive to the domain thickness, 12 versus 20, but the
% % % shape of
% %
% % ice1_12 = ice1;
% % ice2_12 = ice2;
% %
% % % ice1_20 = ice1;
% % % ice2_20 = ice2;
% %
% % Z_12 = 0.02:0.04:12;
% % Z_20 = 0.02:0.04:20;
% %
% % tidx = month(ice1.Time) == 3; % & hour(ice1.Time) == 15;
% %
% % figure; hold on
% % plot(mean(ice2_12.Tice(:, tidx), 2), Z_12);
% % plot(mean(ice2_20.Tice(:, tidx), 2), Z_20);
% % legend('12', '20')
% % set(gca, 'YDir', 'reverse')
% %
% % figure
% % scatterfit(ice1_12.chf, ice1_20.chf)
% %
% % % the ones from december are 500 (20 m depth)
% % % the ones from february are 300 (12 m depth)
% % % they both have the same grid resolution, so just keep the top 300 cells
% %
% %
% %
% % ice2_1 = load('/Volumes/Samsung_T5b/icemodel/output/v10b/sector/skinmodel/mar/zobs/2009/ice2_1.mat').('ice2');
% % ice2_2 = load('/Volumes/Samsung_T5b/icemodel/output/v10b/sector/skinmodel/mar/zobs/2009/ice2_2.mat').('ice2');
% %
% % size(ice2_1.Tice)
% % size(ice2_2.Tice)
% %
% % 12/300
% % 20/500
% %
% %
% %
% % %%
% % %    test = ncread(outfilename,'f_liq');
% % %    test1 = squeeze(test(1,:,:));
% % %    test2 = squeeze(test(2,:,:));
% % %    test3 = squeeze(test(3,:,:));
% % %
% % %    figure; plot(test1(1,:)); hold on; plot(test1(2,:)); plot(test1(3,:));
% %
% %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % this shows that runoff cannot be computed from f_liq after the fact,
% % % % whereas it can be computed from df_liq, so I need to save df_liq. I
% % % % should be able to compute ro_sno, cp_sno, etc from f_liq and f_ice, so I
% % % % still need to save them along with Tice
% % %
% % %
% % % vzero = zeros(size(ice2.f_liq(:,1)));
% % % df_test = [vzero diff(ice2.f_liq,1,2)];
% % %
% % % dz = 0.04;
% % % runoff_test = tocolumn(cumsum( sum( 4*dz.*df_test ) ));
% % % % runoff_test = tocolumn(cumsum( sum( dz.*ice2.df_liq ) ));
% % %
% % % figure;
% % % plot(ice1.Time,ice1.runoff); hold on;
% % % plot(ice1.Time,4.*runoff_test,':');
