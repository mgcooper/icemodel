clean

% TODO: Save to nc and delete ice1 files:
% icemodel, mar, zobs
% icemodel, modis, zobs
% icemodel, mar
% icemodel, modis
% skinmodel, mar, zobs
% skinmodel, modis, zobs
% DONE skinmodel, mar
% DONE skinmodel, modis

% to finish the makencfile
% - add grid_mapping
% - decide how to deal with zobs
% - remove GetDimsFromData stuff esp getdimid it was compact and clear
% - DONE separate into ice1 and ice2 files, latter on monthly timestep?
% - DONE add a gridcell variable
% - DONE figure out the auxiliary coordinate variable thing
% - DONE add runoff to the verify script
% - DONE move ncells, nhrs, nlayrs into getdimdata


% FINAL DECISION ON PRECISION
% - Save all nc files to double precision, it will still save a lot of disk
% space plus eliminate unneccesary variables e.g. skinmodel freeze
% - Then I can decide how to handle zobs, and swap out data as needed
% - For the final archive, if needed I can save as single but by then, I won't
% need to
% - The tricky decision is whether to save Qm and Qe, but since I can regenerate
% skinmodel data quickly, just save the core data. RECALL that all variables can
% be generated for individual sites, this is just the runoff for the sector
% grids

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
userdata = 'modis'; % {'mar', 'modis'};
siteopts = setBasinOpts('sitename', sitename, 'simmodel', simmodel, 'userdata', userdata);

testname = 'zobs';

% simyears = 2009;

% Set path to data
pathdata = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename, simmodel, userdata, testname);
pathsave = pathdata;

%%

% [info, data] = icemodel.netcdf.makencfile(pathdata, pathsave, simyears);
% icemodel.netcdf.makencfile(pathdata, pathsave, simyears(2:end), ...
%    deflateLevel=9, test_write=false);

% ice 1
icemodel.netcdf.makencfile('ice1', pathdata, pathsave, simmodel, ...
   forcings, userdata, 'sw', simyears, xtype='NC_DOUBLE', ...
   deflateLevel=9, testwrite=false);

% ice 2
% icemodel.netcdf.makencfile('ice2', pathdata, pathsave, simmodel, ...
%    forcings, userdata, 'sw', simyears, xtype='NC_DOUBLE', ...
%    deflateLevel=9, testwrite=false, Z=12, dz=0.04);


%% Check x, y, lat, lon, grid cell, time

filelist = listfiles(pathdata, pattern="ice2", aslist=true, fullpath=true);
f = filelist{2};
data = ncreaddata(f);
% icemodel.netcdf.plotncfile(f)

verify_ncfile(pathdata, 'ice1')

%%

% <model>.<forcings>.<albedo>.<sitename>.<yyyy>.nc4
% skinmodel.mar3p11.mar3p11.sw.2009.nc4
% skinmodel.mar3p11.modis.sw.2009.nc4
%
% skinmodel.2009.sw.mar3p11.mar3p11.nc4
% skinmodel.2009.sw.mar3p11.modis.nc4
%
% skinmodel.mar_forcings.modis_albedo.sw.2009.nc4


%%

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
