function varargout = setpath(pathstr,varargin)
%SETPATH return full path string or (if no output request) add the full path
% 
%  pathstr = setpath(pathstr,'project') returns the full path to input
%  pathstr located in MATLABPROJECTPATH
% 
%  string (char) 
% 
% See also


%------------------------------------------------------------------------------
% input parsing
validpathtypes = @(x)any(validatestring(x,{'matlab','project','data','user'}));
validpostset   = @(x)any(validatestring(x,{'goto','none'}));
p              = magicParser;
p.FunctionName = 'setpath';

p.addRequired( 'pathstr',              @(x)ischar(x)|isstruct(x)  );
p.addOptional( 'pathtype', 'matlab',   validpathtypes             );
p.addOptional( 'postset',  'none',     validpostset               );

p.parseMagically('caller');

pathtype = p.Results.pathtype;
postset  = p.Results.postset;
%------------------------------------------------------------------------------

% first determine what type of path (matlab, project, or data path). NOTE: the
% 'otherwise' option was supposed to negate the need to pass 'user' for a full
% path, but it doesn't work b/c default option is 'matlab' for convenience
switch pathtype
   case 'matlab'
      pathroot = getenv('MATLABUSERPATH');
   case 'project'
      pathroot = getenv('MATLABPROJECTPATH');
   case 'data'
      pathroot = getenv('USERDATAPATH');
   case 'user'
      pathroot = '';       % assume the full path is passed in
   otherwise
      pathroot = pathtype; % assume the root is passed in
end

if isstruct(pathstr)
   fields = fieldnames(pathstr);

   for n = 1:length(fields)
      pathname = [pathroot pathstr.(fields{n})];
      if ~contains(pathstr,'.')
         if pathname(end) ~= "/"
            pathname = strcat(pathname,'/');
         end
      end
      pathstr.(fields{n}) = pathname;
   end

elseif ischar(pathstr)
   pathstr = [pathroot pathstr];

   if ~contains(pathstr,'.')
      if pathstr(end) ~= "/"
         pathstr = strcat(pathstr,'/');
      end
   end
end

% if postset = 'goto', cd to the path
if strcmp(postset,'goto')
   cd(pathstr)
end

% parse outputs
switch nargout
   case 1
      % return the full path
      varargout{1} = pathstr;
   case 0
      % add the path, don't return it
      addpath(genpath(pathstr));
end




% the old way that required knowing which computer I was on:

% homepath = pwd;
%
% if isstruct(path)
%     fields = fieldnames(path);
%
%     for n = 1:length(fields)
%
%         pathname = path.(fields{n});
%
%         if strcmp(homepath(8:14),'coop558')
%             pathname = ['/Users/coop558/Dropbox/MATLAB/' pathname];
%         elseif strcmp(homepath(8:17),'mattcooper')
%             pathname = ['/Users/mattcooper/Dropbox/MATLAB/' pathname];
%         end
%
%         path.(fields{n}) = pathname;
%     end
%
% elseif ischar(path)
%
%     if strcmp(homepath(8:14),'coop558')
%             path = ['/Users/coop558/Dropbox/MATLAB/' path];
%     elseif strcmp(homepath(8:17),'mattcooper')
%             path = ['/Users/mattcooper/Dropbox/MATLAB/' path];
%     end
%
% end
%
%
%
