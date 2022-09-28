clean

%------------------------------------------------------------------------------
% send the latest icemodel dependencies to maxim
%------------------------------------------------------------------------------
funcname = '/Users/coop558/myprojects/matlab/icemodel/drive/b_drive.m';
pathsave = '/Users/coop558/myprojects/matlab/sandbox/icemodel/';

workoff skinmodel
workon icemodel
workon runoff

info     = buildsandbox(funcname,'pathsave',pathsave,'dryrun',true);

%------------------------------------------------------------------------------
% repeat for skinmdoel
%------------------------------------------------------------------------------
funcname = '/Users/coop558/myprojects/matlab/skinmodel/drive/b_drive.m';
pathsave = '/Users/coop558/myprojects/matlab/sandbox/skinmodel/';

workon skinmodel
workon runoff
workoff icemodel
info     = buildsandbox(funcname,'pathsave',pathsave,'dryrun',false,'strexclude','Cupid');


%------------------------------------------------------------------------------
% use this to move unneccesary files out of functions/ and into private/
%------------------------------------------------------------------------------
srcdir1  = '/Users/coop558/myprojects/matlab/sandbox/icemodel/functions/';
srcdir2  = '/Users/coop558/myprojects/matlab/sandbox/skinmodel/functions/';
cpydir   = '/Users/coop558/myprojects/matlab/icemodel/functions/';
dstdir   = '/Users/coop558/myprojects/matlab/icemodel/private/functions/';

% these are lists of functions in the sandboxes
srclist1 = getlist(srcdir1,'.m');
srclist2 = getlist(srcdir2,'.m');
srclist1 = {srclist1.name}';
srclist2 = {srclist2.name}';

% this is the combined list of functions in the sandboxes
srclist  = unique([srclist1;srclist2]);

% this is the list of functions in the main icemodel/functions dir
cpylist  = getlist(cpydir,'.m');
cpylist  = {cpylist.name}';

% this is the list of functions that are already in icemodel/private/functions 
dstlist  = getlist(dstdir,'.m');
dstlist  = {dstlist.name}';

for n = 1:numel(cpylist)
   
   % if the file in the main icemodel function dir isn't in the sandbox, it
   % means we don't need it to run either icemodel or skinmodel, so move it to
   % the private/functions dir
   if ~ismember(cpylist(n),srclist)
      %fprintf(['moving ' [cpydir cpylist{n}] ' to ' [dstdir cpylist{n}] '\n\n']);
      if ~ismember(cpylist(n),dstlist)
         %system(['git mv ' [cpydir cpylist{n}] ' ' [dstdir cpylist{n}]])
         movefile([cpydir cpylist{n}],[dstdir cpylist{n}]);
      else
         fprintf(['file ' [cpydir cpylist{n}] ' already present in private/']);
      end
   end
end

% % send mk_snowMask to my other computer
% funcname = 'mk_snowMask';
% pathsave = '/Volumes/GoogleDrive/My Drive/TEMP/';
% 
% info     = buildsandbox(funcname,'pathsave',pathsave,'dryrun',false);

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
