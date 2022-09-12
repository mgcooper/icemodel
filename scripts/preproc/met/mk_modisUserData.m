clean

% this saves modis data for the catchments, so a situation like
% opts.metfname = 'met_kanm_modis_2016.mat' will swap out the
% catchment-scale modis albedo for the point=scale kanm albedo, whereas the
% setup prior to this would have just used the point-scale modis albedo
% saved in the kanm metfile, but upon inspection this is not good b/c at
% the single pixel scale, the modis data can flatline for extended periods
% of time, whereas averagig over the catchment reduces the impact of this

% need files like this:
% input/userData/behar/KANM_behar_2016.mat

srcpath = '/Users/coop558/MATLAB/GREENLAND/icemodel/input/met/';
% dstpath = '/Users/coop558/MATLAB/GREENLAND/icemodel/model/input/met/';
dstpath = '/Users/coop558/MATLAB/GREENLAND/icemodel/input/userData/';
cd(srcpath)

% want to copy: 
% met_KANM_KANM_YYYY to met_behar_KANM_

% DONE: slv1, slv2, behar, upperBasin, ak4, KANL, KANM
years    = 2014:2016;
siteName = 'hills'; %'KANL';
dataName = 'MAR'; %'KANL';
fileCopy = ['met_' siteName '_' dataName '_'];
fileDest = ['modis_' siteName '_'];

for n = 1:numel(years)
   thisYear = num2str(years(1)+n-1);
   srcfile  = [srcpath fileCopy thisYear '_1hr.mat'];
   dstfile  = [dstpath fileDest thisYear '.mat'];
   
   % actually need to rename from 'met' to 'Data', not just copy
   load(srcfile,'met'); Time=met.Time; modis=met.MODIS;
   Data = timetable(Time,modis);
   save(dstfile,'Data');
   %copyfile(srcfile,dstfile);
end


% below is the original usage of this script

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Behar, 2015, KANM albedo
% srcfile = 'met_KANM_KANMforcing_KANMalbedo_2015_';
% dstfile = 'met_behar_KANMforcing_KANMalbedo_2015_';
% 
% copyfile([srcfile '1hr.mat'],[dstfile '1hr.mat']);
% copyfile([srcfile '15m.mat'],[dstfile '15m.mat']);
% 
% % Behar, 2016, KANM albedo
% srcfile = 'met_KANM_KANMforcing_KANMalbedo_2016_';
% dstfile = 'met_behar_KANMforcing_KANMalbedo_2016_';
% 
% copyfile([srcfile '1hr.mat'],[dstfile '1hr.mat']);
% copyfile([srcfile '15m.mat'],[dstfile '15m.mat']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Behar, 2015, MODIS albedo
% srcfile = 'met_KANM_KANMforcing_MODISalbedo_2015_';
% dstfile = 'met_behar_KANMforcing_MODISalbedo_2015_';
% 
% copyfile([srcfile '1hr.mat'],[dstfile '1hr.mat']);
% copyfile([srcfile '15m.mat'],[dstfile '15m.mat']);
% 
% % Behar, 2016, KANM albedo
% srcfile = 'met_KANM_KANMforcing_MODISalbedo_2016_';
% dstfile = 'met_behar_KANMforcing_MODISalbedo_2016_';
% 
% copyfile([srcfile '1hr.mat'],[dstfile '1hr.mat']);
% copyfile([srcfile '15m.mat'],[dstfile '15m.mat']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Upper Basin, 2016, KANL albedo
% srcfile = 'met_KANL_KANLforcing_KANLalbedo_2016_';
% dstfile = 'met_upperBasin_KANLforcing_KANLalbedo_2016_';
% 
% copyfile([srcfile '1hr.mat'],[dstfile '1hr.mat']);
% copyfile([srcfile '15m.mat'],[dstfile '15m.mat']);
% 
% % Upper Basin, 2016, MODIS albedo
% srcfile = 'met_KANL_KANLforcing_MODISalbedo_2016_';
% dstfile = 'met_upperBasin_KANLforcing_MODISalbedo_2016_';
% 
% copyfile([srcfile '1hr.mat'],[dstfile '1hr.mat']);
% copyfile([srcfile '15m.mat'],[dstfile '15m.mat']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
