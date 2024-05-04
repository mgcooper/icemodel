function pathname = projectpath(varargin)
   %PROJECTPATH Return the full path to the top-level project directory.
   %
   %  PATHNAME = PROJECTPATH(VARARGIN)
   %
   % See also

   pathname = fileparts(toolboxpath());

   if nargin == 1
      pathname = fullfile(pathname, varargin{:});
   end
end
