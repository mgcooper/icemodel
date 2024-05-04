function pathname = fullpath(varargin)
   %PATHNAME Build full path to toolbox folder or file.
   %
   %  PATHNAME = FULLPATH() Returns the full path to the top level folder.
   %  PATHNAME = FULLPATH(FOLDERNAME)
   %  PATHNAME = FULLPATH(FOLDERNAME, FILENAME)
   %
   % Description
   %  PATHNAME = FULLPATH() Returns the full path to the top level (project)
   %  folder.
   %
   %  PATHNAME = FULLPATH(FOLDERNAME) Returns the full path to the
   %  toolbox/FOLDERNAME folder.
   %
   %  PATHNAME = FULLPATH(FOLDERNAME, FILENAME) Returns the full path to the
   %  toolbox/FOLDERNAME/FILENAME file.
   %
   %  PATHNAME = FULLPATH(varargin) Recursively appends the contents of varargin
   %  to the toolbox/ folder.
   %
   % See also:

   % Ensure inputs are chars
   [varargin{:}] = convertStringsToChars(varargin{:});
   cellfun(@(arg) validateattributes(arg, {'char'}, {'row'}), varargin)

   % Top level path
   pathname = projectpath();

   % Recursively validate sub-folders and append them to the path string.
   for subpath = varargin(:)'
      pathname = fullfile(pathname, subpath{:});
   end
end
