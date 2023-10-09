function fullpath = projectpath(varargin)
%PROJECTPATH return the full path to the top-levl project directory. 
% 
% 
% 
% See also

fullpath = fileparts(fileparts(fileparts(mfilename('fullpath'))));

if nargin == 1
   fullpath = fullfile(fullpath, varargin{:});
end

