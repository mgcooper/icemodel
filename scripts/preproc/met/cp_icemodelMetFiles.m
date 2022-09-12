clean

% in addition to renaming certain files such as KANL_ to upperBasin, since
% they're identical, this also copies the KAN met files to the 'user data'
% folders so the forcingUserData method works to swap out met station vars
% e.g. if we want to force with MAR but swap in KanM albedo

% need files like this:
% input/userData/behar/KANM_behar_2016.mat

srcpath = '/Users/coop558/MATLAB/GREENLAND/icemodel/input/met/';
% dstpath = '/Users/coop558/MATLAB/GREENLAND/icemodel/model/input/met/';
dstpath = '/Users/coop558/MATLAB/GREENLAND/icemodel/input/userData/';
cd(srcpath)


% DONE: Upper Basin, behar, 
years    = 2015;
siteName = 'slv2';
copyName = 'KANM';
fileCopy = ['met_' copyName '_' copyName 'forcing_' copyName 'albedo_'];
fileDest = [siteName '/' copyName '_' siteName '_'];

for n = 1:numel(years)
   thisYear = num2str(years(1)+n-1);
   srcfile  = [srcpath fileCopy thisYear '_1hr.mat'];
   dstfile  = [dstpath fileDest thisYear '.mat'];
   
   % actually need to rename from 'met' to 'Data', not just copy
   load(srcfile,'met'); Data = met;
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
