clean

% send mk_snowMask to my other computer
funcname = 'mk_snowMask';
pathsave = '/Volumes/GoogleDrive/My Drive/TEMP/';

info     = buildsandbox(funcname,'pathsave',pathsave,'dryrun',false);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % this is the oriignal setup which was then copied to personal mac:
% funcpath = setpath('GREENLAND/icemodel/model/experiment/v10/');
% funcname = [funcpath 'a_icemodel_drive.m'];
% pathsave = '/Users/coop558/myprojects/matlab/icemodel/v10/';
% 
% info     = buildsandbox(funcname,'pathsave',pathsave,'strexclude','Cupid');
% % info     = buildsandbox(funcname,'pathsave',pathsave);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % this is to check for diff's b/w experiment/v10 and versions/v10
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % first build a sandbox based on experiment/v10
% funcpath = setpath('GREENLAND/icemodel/model/experiment/v10/');
% oldpath  = addpath(genpath(funcpath));
% funcname = [funcpath 'b_drive_region.m'];
% pathsave = '/Users/coop558/myprojects/matlab/icemodel/v10_a/';
% info     = buildsandbox(funcname,'pathsave',pathsave,'strexclude','Cupid');
% rmpath(genpath(funcpath)); % remove from path before building the next one
% 
% % now build a sandbox based on versions/v10
% funcpath = setpath('GREENLAND/icemodel/model/versions/v10/');
% oldpath  = addpath(genpath(funcpath));
% funcname = [funcpath 'b_drive.m'];
% pathsave = '/Users/coop558/myprojects/matlab/icemodel/v10_b/';
% info     = buildsandbox(funcname,'pathsave',pathsave,'strexclude','Cupid');
% rmpath(genpath(funcpath)); % remove from path before building the next one
% 
