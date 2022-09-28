function varargout = setpath(pathstr,varargin)
   
%------------------------------------------------------------------------------
% input parsing
p = MipInputParser;
p.FunctionName='setpath';
p.addRequired('path',@(x)ischar(x)|isstruct(x));
p.addOptional('pathtype','matlab',@(x)ischar(x));
p.parseMagically('caller');
pathtype=p.Results.pathtype;
%------------------------------------------------------------------------------

% first determine what type of path (matlab, project, or data path)
switch pathtype
   case 'matlab'
      pathroot = getenv('MATLABUSERPATH');
   case 'project'
      pathroot = getenv('MATLABPROJECTPATH');
   case 'data'
      pathroot = getenv('USERDATAPATH');
   otherwise
      pathroot = pathtype; % assume the root is passed in
end

if isstruct(pathstr)
   fields = fieldnames(pathstr);

   for n = 1:length(fields)
      pathname = [pathroot pathstr.(fields{n})];
      if pathname(end) ~= "/"
         pathname = strcat(pathname,'/');
      end
      pathstr.(fields{n}) = pathname;
   end

elseif ischar(pathstr)
   pathstr = [pathroot pathstr];

   if pathstr(end) ~= "/"
      pathstr = strcat(pathstr,'/');
   end
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
