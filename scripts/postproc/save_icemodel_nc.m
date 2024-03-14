clean

% to finish the makencfile
% - add grid_mapping
% - decide how to deal with zobs
% - separate into ice1 and ice2 files, latter on monthly timestep?
% - DONE add a gridcell variable
% - DONE figure out the auxiliary coordinate variable thing

% To deal with zobs, pick N random points from icemask X, Y, then run them with
% runscript runpoint functionality using no rounding, compare with the saved
% data in zobs ... if they're equivalent then there's no getting around it, I
% should just delete the old data and use the new data. BUT, might check
% sensitivity to numiter, bc, etc. 

% file size in mb/cell for deflateLevel = 1,5,9 w/wo 1d or 0d vars
% 1: 10.9 (10.8 w/o 1-d or 0-d vars)
% 5: 9.4
% 9: 9.2 (9.1 w/o 1-d or 0-d vars)
% % I also tested if compression is better w/o the 1-d and 2-d vars - it isn't

savedata = true;
sitename = 'sector';
simyears = 2009:2018;
simmodel = 'skinmodel'; % {'icemodel', 'skinmodel'};
forcings = 'mar';
userdata = 'mar'; % {'mar', 'modis'};
siteopts = setBasinOpts('sitename', sitename, 'simmodel', simmodel, 'userdata', userdata);

make_backups = false;

simyears = 2009;

% Set path to data
pathdata = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename, simmodel, userdata);
pathsave = pathdata;

%%

% [info, data] = icemodel.netcdf.makencfile(pathdata, pathsave, simyears);
% icemodel.netcdf.makencfile(pathdata, pathsave, simyears(2:end), ...
%    deflateLevel=9, test_write=false);

icemodel.netcdf.makencfile(pathdata, pathsave, simmodel, ...
   forcings, userdata, simyears, whichdata="ice1", ...
   deflateLevel=9, test_write=false);

%% Check x, y, lat, lon, grid cell, time

filelist = listfiles(pathdata, pattern="nc4", aslist=true, fullpath=true);
f = filelist{1};

info = ncparse(f);
data = ncreaddata(f);
x = data.x_easting;
y = data.y_northing;

% Plot all 1d vars
vars = {'Tsfc', 'melt', 'freeze', 'runoff', 'cond', 'subl'};
figontop
for n = 1:numel(vars)
   % v = ncread(f, vars{n});
   v = data.(vars{n});
   scatter(x, y, 20, mean(v, 2), 'filled')
   colorbar
   title(vars{n})
   pause
   clf
end

% %% Check Tice
%
% % Can't read past 3896
% % tice = ncread(f, 'Tice', [1 1 1], [1 500 8760]);
% tice = squeeze(ncread(f, 'Tice', [1 1 1], [1 500 3896]));
%
% % This confirms that up to 3896, the data is oriented correctly
% figure
% plot(mean(tice, 2), depth)
% set(gca, 'YDir', 'Reverse')
%
% % Pick some random indices
% figure; hold on
% plot(tice(:, 1:500:3500), depth)
% set(gca, 'YDir', 'Reverse')
% legend(string((1:500:3500).'))
%
% %%
% tice = squeeze(ncread(f, 'Tice', [1 1 1], [2479 500 1]));
%
% % the data cannot be loaded beyond timestep 3896, so check the file
% 3896*40
%
% % if the first timestep is read, then the data comes in up to cell 40:
% 40 * 500
%
% % try to read beyond grid cell 40
% tice = squeeze(ncread(f, 'Tice', [40 1 1], [1 500 3896]));
%
% % the chunk size is [138,28,487]
% % the dimensions are [2479,500,8760]
% 3896/487 % 8
% 2479/138 % not even
% 500/28 % not even
%
% % so 8 chunks along the time dimension were read
%
% %% This confirmed that the transpose is WRONG
% % info_transpose = ncinfo(fullfile(pathdata, 'icemodel_2009.nc4'));
% % info_no_transpose = ncinfo(fullfile(pathdata, 'icemodel_2009_no_transpose.nc4'));
% %
% % Tice_transpose = ncread(fullfile(pathdata, 'icemodel_2009.nc4'), 'Tice');
% % Tice_no_transpose = ncread(fullfile(pathdata, 'icemodel_2009_no_transpose.nc4'), 'Tice');
% %
% % isequal(Tice_transpose, Tice_no_transpose)
% %
% % onecell_transpose = squeeze(Tice_transpose(1, :, :));
% % onecell_no_transpose = squeeze(Tice_no_transpose(1, :, :));
% %
% % figure; hold on
% % plot(mean(onecell_transpose, 2), 1:500)
% % plot(mean(onecell_no_transpose, 2), 1:500)
% % set(gca, 'YDir', 'reverse')
% % legend('transpose', 'no transpose')
% %
% % % now permute as one would
% % % onecell_transpose = flipud(permute(onecell_transpose, [2 1 3]));
% % % onecell_no_transpose = flipud(permute(onecell_no_transpose, [2 1 3]));
% % onecell_transpose = flipud(permute(onecell_transpose, [2 1 3]));
% % onecell_no_transpose = flipud(permute(onecell_no_transpose, [2 1 3]));
% %
% % figure; hold on
% % plot(mean(onecell_transpose, 1), 1:500)
% % plot(mean(onecell_no_transpose, 1), 1:500)
% % set(gca, 'YDir', 'reverse')
% % legend('transpose', 'no transpose')
% %
% %
% % % <model>.<forcings>.<albedo>.<sitename>.<yyyy>.nc4
% % % skinmodel.mar3p11.mar3p11.sw.2009.nc4
% % % skinmodel.mar3p11.modis.sw.2009.nc4
% % %
% % % skinmodel.2009.sw.mar3p11.mar3p11.nc4
% % % skinmodel.2009.sw.mar3p11.modis.nc4
% % %
% % % skinmodel.mar_forcings.modis_albedo.sw.2009.nc4
% %
% % % tic
% % % f_ice = ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'f_ice');
% % % toc
% % %
% % % numhrs = 8784;
% % % numgrids = 1487;
% % %
% % % % Takes about 0.92 seconds to load the data
% % % tic
% % % data = load(fullfile(fileparts(fileparts(pathdata)), 'icemodel_mar_2008.mat'));
% % % toc
% % %
% % % % Takes about 0.88 seconds to create the equivalent data from the nc4 files
% % % tic
% % % data.M = nan(numhrs, numgrids);
% % % data.F = nan(numhrs, numgrids);
% % % data.R = nan(numhrs, numgrids);
% % % data.Tsfc = nan(numhrs, numgrids);
% % %
% % % data.M = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'depth_melt'));
% % % data.F = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'depth_freeze'));
% % % data.R = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'runoff'));
% % % data.Tsfc = transpose(ncread(fullfile(pathdata, 'icemodel_2008.nc4'), 'Tsfc'));
% % % toc
% % %
% % % %%
% % %
% % %
% % % % skinmodel is not sensitive to the domain thickness, 12 versus 20, but the
% % % % shape of
% % %
% % % ice1_12 = ice1;
% % % ice2_12 = ice2;
% % %
% % % % ice1_20 = ice1;
% % % % ice2_20 = ice2;
% % %
% % % Z_12 = 0.02:0.04:12;
% % % Z_20 = 0.02:0.04:20;
% % %
% % % tidx = month(ice1.Time) == 3; % & hour(ice1.Time) == 15;
% % %
% % % figure; hold on
% % % plot(mean(ice2_12.Tice(:, tidx), 2), Z_12);
% % % plot(mean(ice2_20.Tice(:, tidx), 2), Z_20);
% % % legend('12', '20')
% % % set(gca, 'YDir', 'reverse')
% % %
% % % figure
% % % scatterfit(ice1_12.chf, ice1_20.chf)
% % %
% % % % the ones from december are 500 (20 m depth)
% % % % the ones from february are 300 (12 m depth)
% % % % they both have the same grid resolution, so just keep the top 300 cells
% % %
% % %
% % %
% % % ice2_1 = load('/Volumes/Samsung_T5b/icemodel/output/v10b/sector/skinmodel/mar/zobs/2009/ice2_1.mat').('ice2');
% % % ice2_2 = load('/Volumes/Samsung_T5b/icemodel/output/v10b/sector/skinmodel/mar/zobs/2009/ice2_2.mat').('ice2');
% % %
% % % size(ice2_1.Tice)
% % % size(ice2_2.Tice)
% % %
% % % 12/300
% % % 20/500
% % %
% % %
% % %
% % % %%
% % % %    test = ncread(outfilename,'f_liq');
% % % %    test1 = squeeze(test(1,:,:));
% % % %    test2 = squeeze(test(2,:,:));
% % % %    test3 = squeeze(test(3,:,:));
% % % %
% % % %    figure; plot(test1(1,:)); hold on; plot(test1(2,:)); plot(test1(3,:));
% % %
% % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % this shows that runoff cannot be computed from f_liq after the fact,
% % % % % whereas it can be computed from df_liq, so I need to save df_liq. I
% % % % % should be able to compute ro_sno, cp_sno, etc from f_liq and f_ice, so I
% % % % % still need to save them along with Tice
% % % %
% % % %
% % % % vzero = zeros(size(ice2.f_liq(:,1)));
% % % % df_test = [vzero diff(ice2.f_liq,1,2)];
% % % %
% % % % dz = 0.04;
% % % % runoff_test = tocolumn(cumsum( sum( 4*dz.*df_test ) ));
% % % % % runoff_test = tocolumn(cumsum( sum( dz.*ice2.df_liq ) ));
% % % %
% % % % figure;
% % % % plot(ice1.Time,ice1.runoff); hold on;
% % % % plot(ice1.Time,4.*runoff_test,':');
