function varargout = listfiles(folderlist, opts)
   %LISTFILES List all files in folder and (optionally) subfolders.
   %
   % FILELIST = LISTFILES(FOLDERNAME)
   % FILELIST = LISTFILES(FOLDERLIST)
   % FILELIST = LISTFILES(_, 'ASLIST', TRUE)
   % FILELIST = LISTFILES(_, 'ASSTRUCT', TRUE)
   % FILELIST = LISTFILES(_, 'ASSTRING', TRUE)
   % FILELIST = LISTFILES(_, 'SUBFOLDERS', TRUE)
   % FILELIST = LISTFILES(_, 'FULLPATH', TRUE)
   % FILELIST = LISTFILES(_, 'MFILES', TRUE)
   % FILELIST = LISTFILES(_, 'RMPATTERN', RMPATTERN)
   % [FILELIST, NUMFILES] = LISTFILES(_)
   %
   % Description:
   %
   % FILELIST = LISTFILES(FOLDERNAME) Returns a directory structure in the
   % same format returned by DIR containing all files in directory FOLDERNAME.
   %
   % [FILELIST, NUMFILES] = LISTFILES(FOLDERNAME) Also returns the number of
   % files.
   %
   % FILELIST = LISTFILES(_, 'SUBFOLDERS', TRUE) Returns all files in all
   % subfolders of FOLDERNAME.
   %
   % FILELIST = LISTFILES(_, 'MFILES', TRUE) Only returns m-files.
   % FILELIST = LISTFILES(_, 'MFILES', TRUE, 'MATFILES', TRUE) Also returns
   % mat-files.
   %
   % FILELIST = LISTFILES(_, 'ASLIST', TRUE) Returns a list of
   % filenames in a cell array rather than the default directory struct.
   %
   % FILELIST = LISTFILES(_, 'ASLIST', TRUE, 'ASSTRING', TRUE) Returns
   % a string array of filenames rather than the default directory struct.
   %
   % FILELIST = LISTFILES(_, 'ASLIST', TRUE, 'FULLPATHS', TRUE) Returns
   % fullpaths to each file in a cell array rather than the default directory
   % struct.
   %
   % Copyright (c) 2023, Matt Cooper, BSD 3-Clause License, github.com/mgcooper.
   %
   % See also: mfilename, mcallername, mfoldername, ismfile

   arguments
      folderlist (:,1) string = pwd()
      opts.subfolders (1,1) logical = false
      opts.asstruct (1,1) logical = true
      opts.aslist (1,1) logical = false
      opts.asstring (1,1) logical = false
      opts.fullpath (1,1) logical = false
      opts.pattern (1,1) string = "*"
      opts.mfiles (1,1) logical = false
      opts.matfiles (1,1) logical = false
      opts.rmpatterns (1,:) string = ".git"
      opts.ignoredfolders (1, :) string = []
   end

   opts.pattern = parseFilePattern(opts);

   if opts.aslist || opts.asstring
      opts.asstruct = false;
   end

   % Create the list of files
   list = arraymap(@(folder) processOneFolder(folder, opts), folderlist);

   % I think this was added to stack dir structs, but it failed when list was a
   % cell array of paths i.e. opts.aslist == true, so I added opts.asstruct.
   % Update: When folder is non-scalar, list is a cell array where each element
   % is a cell array of files in a subfolder. I think that's why I had vertcat,
   % and it failed when list was scalar b/c it tried to vertcat all the
   % individual paths into one char array. So I added the numel > 1 check.
   % But I am not sure if other use cases depend on the files being in separate
   % cell elements, so I commented that out and added the aslist check.
   if iscell(list)
      if opts.asstruct % || numel(list) > 1
         list = vertcat(list{:});
      elseif numel(list) > 1 && opts.aslist
         list = vertcat(list{:});
      end
   end

   % Parse outputs
   switch nargout
      case 1
         varargout{1} = list;
      case 2
         varargout{1} = list;
         varargout{2} = numel(list);
   end
end

%% Local Functions

function list = processOneFolder(folder, opts)

   % Get all files in main folder and if requested, sub folders
   list = dir(fullfile(folder, opts.pattern));
   list(strncmp({list.name}, '.', 1)) = [];
   list = list(~[list.isdir]);

   if opts.mfiles
      if opts.matfiles
         list = list(...
            strncmp(reverse({list.name}), 'm.', 2) | ...
            strncmp(reverse({list.name}), 'tam.', 4));
      else
         list = list(strncmp(reverse({list.name}), 'm.', 2));
      end
   end

   % Remove files containing the "RMPATTERNS"
   list = trimfiles(list, opts);

   if opts.aslist

      if opts.fullpath
         list = transpose(fullfile({list.folder}, {list.name}));
      else
         list = transpose(fullfile({list.name}));
      end

      % convert to string array if requested
      if opts.asstring
         list = string(list);
      end
   end
end

function pattern = parseFilePattern(opts)
   %PARSEFILEPATTERN Create a wildcard pattern from the user supplied one

   pattern = opts.pattern;

   % Wildcard * is appended to the pattern by default, so remove if supplied
   if startsWith(pattern, "*") || endsWith(pattern, "*")
      pattern = erase(pattern, "*");
   end
   if startsWith(pattern, ".")
      pattern = erase(pattern, ".");
   end
   if opts.subfolders
      if pattern ~= ""
         pattern = strcat("**/*", pattern, "*");
      else
         pattern = "**/*"; % Prevent strcat from creating **/**
      end
   else
      if pattern ~= ""
         pattern = strcat("*", pattern, "*");
      else
         pattern = "*"; % Prevent strcat from creating ** which returns subfolders
      end
   end
end

function list = trimfiles(list, opts)

   % Note: call trimfiles prior to converting from dir struct to list to
   % respect the requested output type, but perform the same "fullpath"
   % check that is performed in the "aslist" section to respect the
   % requested output type meaning the rmpattern must apply to the filename
   % OR the fullpath. note - could be preferable to add a "nameonly" option
   % to distinguish whether the rmpatterns apply to the fullpath or filename.

   % Append default patterns that should be removed
   rmpatterns = [opts.rmpatterns, "~$"]; % add more as needed

   if opts.fullpath
      files = transpose(fullfile({list.folder}, {list.name}));
   else
      files = transpose(fullfile({list.name}));
   end

   list = list(~contains(files, rmpatterns));
end
