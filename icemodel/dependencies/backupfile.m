function varargout = backupfile(filename, makecopy, makezip)
   %BACKUPFILE Create a backup file name or folder name and (optionally) a copy.
   %
   % [fullpath_bk, filename_bk] = backupfile(filename) returns the name and
   % full path of the backup file or folder without actually making a copy.
   %
   % [fullpath_bk, filename_bk] = backupfile(filename, true) makes a copy.
   %
   % [fullpath_bk, filename_bk] = backupfile(filename, true, true) makes a copy
   % and moves the copy to a zip file.
   %
   % Note: Provide the full path to avoid unexpected behavior especially when
   % working with private directories or paths not in the MATLAB search space.
   %
   % Examples:
   % 1. backupfile('/Users/user/test.m')
   % 2. backupfile('/Users/user/test_folder', true)
   %
   % See also: tempdir, tempfile
   %
   %#codegen

   if nargin < 2
      makecopy = false;
   end
   if nargin < 3
      makezip = false;
   end

   if ~isoctave && verLessThan('matlab', '9.3') % R2017b
      folderexists = @(folder) exist(folder, 'dir') == 7;
      fileexists = @(file) exist(file, 'file') == 2;
      stringsToChars = @(str) char(str);
   else
      folderexists = @(folder) isfolder(folder);
      fileexists = @(file) isfile(file);
      stringsToChars = @(str) convertStringsToChars(str);
   end

   filename = stringsToChars(filename);
   validateattributes(filename, {'char'}, {'row'}, mfilename, 'FILENAME', 1)
   validateattributes(makecopy, {'logical'}, {'scalar'}, mfilename, 'MAKECOPY', 2)
   validateattributes(makezip, {'logical'}, {'scalar'}, mfilename, 'MAKEZIP', 3)

   % Generate a date-time string to append to the file name.
   filedate = mkfiledate('dd-MMM-yyyy_hh-mm-ss');

   % Remove trailing filesep if it exists.
   filename = rmtrailingsep(filename);

   % Get the parent folder path, filename, and file extension.
   [filepath, filename, fileext] = fileparts(filename);
   fullpath = fullfile(filepath, [filename fileext]);

   % Create a back up file (or folder) name and full path.
   filename_bk = [filename '_' filedate fileext];
   fullpath_bk = fullfile(filepath, filename_bk);

   if makecopy
      if ~fileexists(fullpath) && ~folderexists(fullpath)
         % This lets backupfile be called without if isfile() in the caller
         warning('No backup made. File not found: %s', fullpath)
         return
      end
      if fileexists(fullpath_bk) || folderexists(fullpath_bk)
         warning('No backup made. Backup already exists: %s', fullpath_bk);
      else
         copyfile(fullpath, fullpath_bk);

         if makezip
            zip(fullpath_bk, fullpath_bk);
            if folderexists(fullpath_bk)
               rmdir(fullpath_bk, "s");
            elseif fileexists(fullpath_bk)
               delete(fullpath_bk);
            end
            filename_bk = [filename_bk '.zip'];
            fullpath_bk = [fullpath_bk '.zip'];
         end
         fprintf('Backup created: %s\n', fullpath_bk);
      end
   end

   switch nargout
      case 1
         varargout{1} = fullpath_bk;
         varargout{2} = filename_bk;
   end
end

function filename = rmtrailingsep(filename)
   if strcmp(filename(end), filesep)
      filename(end) = [];
   end
end
