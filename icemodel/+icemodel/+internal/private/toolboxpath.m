function fullpath = toolboxpath(varargin)
   %TOOLBOXPATH Return the full path to the toolbox directory.
   %
   %  FULLPATH = TOOLBOXPATH(VARARGIN)
   %
   % See also: projectpath

   fullpath = fileparts(fileparts(fileparts(fileparts( ...
      mfilename('fullpath')))));

   if nargin == 1
      fullpath = fullfile(fullpath, varargin{:});
   end
end
