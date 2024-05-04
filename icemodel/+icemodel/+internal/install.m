function [requirementsList, urlList] = install(kwargs)
   %INSTALL Install requirements from Github.
   %
   % INSTALL(INSTALLMISSING=TRUE)
   % INSTALL(REQUIREDFILES=REQUIREDFILES)
   % INSTALL(PROJECTFOLDER=PATHNAME)
   % INSTALL(REQUIREMENTSFILE=FILENAME)
   % INSTALL(SAVEREQUIREMENTSFILE=TRUE)
   % INSTALL(LOCALSOURCEFOLDER=PATHNAME)
   % INSTALL(REMOTEREPONAME=REPONAME)
   % INSTALL(REMOTEBRANCH=BRANCHNAME)
   % INSTALL(GITHUBUSERNAME=USERNAME)
   %
   % Input Arguments (optional name-value parameters)
   %
   %  requiredFiles - a list of required functions.
   %  projectPath - folder that needs the required files.
   %  localSourcePath - folder with local versions of the dependencies.
   %  remoteRepoName - Github repo for the localSourcePath.
   %  remoteBranch - Branch to use when downloading from the remote Github repo.
   %
   % Note: This function generates a list of required files to ship with a
   % toolbox, or install them into a folder within the toolbox. Third party
   % users then run an install script which reads the requirements list and
   % installs them from GitHub. This function does not allow the third party to
   % generate the requirements list, because the requirements must be on the
   % local path where this function runs.
   %
   % See also: getRequiredFiles

   arguments
      kwargs.installMissing (1, 1) logical {mustBeNumericOrLogical} ...
         = false

      kwargs.requiredFiles (:, :) string ...
         = []

      kwargs.projectPath (1, :) string {mustBeFolder} ...
         = projectpath()

      kwargs.installPath (1, :) string {mustBeFolder} ...
         = fullfile(toolboxpath(), "dependencies");

      kwargs.ignoreFolder (1, :) string ...
         = "testbed"

      kwargs.requirementsFile (1, :) string ...
         = "requirements.mat"

      kwargs.saveRequirementsFile (1, 1) logical {mustBeNumericOrLogical} ...
         = false

      kwargs.localSourcePath (1, :) {mustBeFolder} ...
         = getenv('MATLABFUNCTIONPATH')

      kwargs.remoteRepoName (1, :) string ...
         = "matfunclib"

      kwargs.remoteBranch (1, :) string ...
         = "dev"

      kwargs.GitHubUserName (1, :) string ...
         = getenv('GITHUB_USER_NAME')
   end

   [projectPath, ignoreFolder, localSourcePath, remoteSourcePath] = ...
      parseargs(kwargs);

   % Remember current folder and go to folder for external dependencies
   job = withcd(projectPath);

   % Find the required files
   if isempty(kwargs.requiredFiles)
      requiredFiles = getRequiredFiles(projectPath, "ignore", ignoreFolder);
      requiredFiles = requiredFiles.missingFiles;
   else
      requiredFiles = kwargs.requiredFiles;
   end

   % Build a url list for the remote files
   [requirementsList, urlList] = remoteDependencyList( ...
      requiredFiles, projectPath, localSourcePath, remoteSourcePath);

   % Option to install the missing requirement locally
   if kwargs.installMissing

      fileList = kwargs.installPath + filesep + requirementsList;
      if ~isfolder(kwargs.installPath)
         mkdir(savePath)
      end
      outFileList = strings(size(requirementsList));
      for n = 1:numel(requirementsList)
         try
            outFileList(n) = websave(fileList(n), urlList(n));
         catch ME
            warning('Failed to download file: %s\nReason: %s', ...
               requirementsList(n), ME.message);
         end
      end
   end
end

%% Local Functions
function [projectPath, ignoreFolder, localSourcePath, ...
      remoteSourcePath, requirementsFile] = parseargs(kwargs)

   % Retrieve the Github user name
   if isempty(kwargs.GitHubUserName)
      error('Set "GitHubUserName" or environment variable "GITHUB_USER_NAME"')
   else
      GITHUB_USER_NAME = kwargs.GitHubUserName;
   end

   % Note: for general use, this should be userpath or MATLABPATH, I think.
   if isempty(kwargs.localSourcePath)
      localSourcePath = userpath();
   else
      localSourcePath = kwargs.localSourcePath;
   end

   if isempty(kwargs.remoteRepoName)
      error('Set "remotesource" or environment variable "MATLABFUNCTIONPATH"')
   else
      GITHUB_URL = 'https://raw.githubusercontent.com/';
      remoteSourcePath = strcat(GITHUB_URL, GITHUB_USER_NAME, '/', ...
         kwargs.remoteRepoName, '/', kwargs.remoteBranch);

      % This works too:
      %GITHUB_URL = 'https://github.com/';
      %remotesource = strcat(GITHUB_URL, GITHUB_USER_NAME, '/', ...
      %   Opts.remoteRepoName, '/raw/', Opts.remotebranch);
   end

   % Pull out required args and remaining optional args
   projectPath = kwargs.projectPath;
   requirementsFile = kwargs.requirementsFile;

   % Full path to ignore folder
   ignoreFolder = fullfile(projectPath, kwargs.ignoreFolder);
end

function [requirementsList, urlList] = remoteDependencyList( ...
      requiredFiles, projectPath, localsource, remotesource)
   %REMOTEDEPENDENCYLIST Get a list of remote url's to function dependencies.

   % This operates on one file at a time

   [requirementsList, urlList] = deal(strings(length(requiredFiles), 1));

   % For each dependency
   for ifile = 1:length(requiredFiles)

      % Get the file name with extension
      [requiredFilePath, requiredFileName, ext] = fileparts(requiredFiles{ifile});
      requiredFileName = strcat(requiredFileName, ext);

      if skipfile(requiredFileName, requiredFilePath, ...
            requirementsList, projectPath)
         continue
      end

      % If the required file exists in the local source repo, add it to the
      % requirementsList and build a full path to the remote file.
      if contains(requiredFilePath, localsource)

         % Add file names to list of external depencies
         requirementsList(ifile) = requiredFileName;

         % Get the subfolder path relative to the top-level source repo
         relativePath = erase(requiredFilePath, localsource);
         relativePath = strrep(relativePath, filesep , '/');
         if relativePath(1) == filesep
            relativePath = relativePath(2:end);
         end

         % Use '/' not fullfile b/c fullfile is platform specific
         urlList(ifile) = remotesource + '/' + relativePath + '/' ...
            + requirementsList(ifile);
      end
   end
   requirementsList(requirementsList == "") = [];
   urlList(urlList == "") = [];
   assert(all(endsWith(urlList, requirementsList)))
end

function tf = skipfile(requiredFileName, requiredFilePath, ...
      requirementsList, projectPath)

   [~, ~, ext] = fileparts(requiredFileName);

   % skip this file if it is the target function, a mex file, already found,
   % or already satisfied b/c it exists in the projectPath.
   tf = ...
      strcmp(ext, '.mex') | ...
      any(strcmpi(requiredFileName, requirementsList)) | ...
      contains(requiredFilePath, projectPath);
end

%%
function Requirements = getRequiredFiles(target, kwargs)
   %GETREQUIREDFILES Get a list of functions required by a function.
   %
   %  REQUIREMENTS = GETREQUIREDFILES(TARGET)
   %  REQUIREMENTS = GETREQUIREDFILES(_, PROJECTPATH=PATHNAME)
   %  REQUIREMENTS = GETREQUIREDFILES(_, REQUIREMENTSFILENAME=FILENAME)
   %  REQUIREMENTS = GETREQUIREDFILES(_, SAVEREQUIREMENTSFILE=TRUE)
   %  REQUIREMENTS = GETREQUIREDFILES(_, IGNOREFILELIST=IGNORELIST)
   %
   %
   % Description
   %
   %  REQUIREMENTS = GETREQUIREDFILES(TARGET) Returns REQUIREMENTS, a struct
   %  containing a list of functions and products required to run the function
   %  or functions specified by the function(s) or folder of functions TARGET.
   %
   %  REQUIREMENTS = GETREQUIREDFILES(TARGET, PROJECTPATH=PATHNAME) Compares
   %  the required functions in flist to functions in the folder specified by
   %  PATHNAME to indicate if required files are missing.
   %
   % Inputs
   %
   %  TARGET - the target function, list of functions, or folder containing
   %  functions that other functions depend upon. TARGET can be a single char or
   %  string, or a cell array or string array of targets. If TARGET is a
   %  filename or list of filenames, they can be fullly qualified paths or file
   %  names, but if they are not full paths, the files must be on the matlab
   %  path or this function must be called from within the folder they are in
   %  for expected behavior. If TARGET is a folder name, it should be a single
   %  folder.
   %
   %  PROJECTPATH - The project top-level path, or path to folder which contains
   %  the files supplied by the project. Since these files are included with the
   %  project, they do not need to be acquired externally.
   %
   % Example
   %
   %
   % Matt Cooper, 23-Dec-2022, https://github.com/mgcooper
   %
   % See also: resolveDependencies, getFunctionConflicts

   arguments
      target (1, :) string
      kwargs.projectPath (1, :) string {mustBeFolder} = projectpath()
      kwargs.requirementsFileName (1, :) string = "requirements.mat"
      kwargs.saveRequirementsFile (1, 1) logical = false
      kwargs.ignore (1, :) string = []
   end

   ignore = kwargs.ignore;
   projectPath = kwargs.projectPath;
   saveRequirementsFile = kwargs.saveRequirementsFile;
   requirementsFileName = kwargs.requirementsFileName;

   % validate each member of the target file / folder list
   target = cellfun(@validateFileList, target, 'Uniform', false);
   ignore = cellfun(@validateFileList, ignore, 'Uniform', false);

   % If target is a folder, convert to file list
   [targetFiles, projectFiles] = prepareFileLists(...
      target, ignore, projectPath);

   % call requiredFilesAndProducts on the file list
   [requiredFiles, requiredProducts] = processFileList(targetFiles);

   % Return the file lists in a struct
   Requirements.requiredFiles = requiredFiles;
   Requirements.requiredProducts = requiredProducts;

   % if a project path is provided, find missing files (required files which
   % are not supplied with the project itself)
   if isempty(projectFiles)
      missingFiles = 'unknown';
   else
      missingFiles = setdiff(requiredFiles, projectFiles);
   end

   Requirements.missingFiles = missingFiles;

   if saveRequirementsFile
      save(requirementsFileName, 'Requirements')
   end
end

%% Local Function
function [requiredFiles, requiredProducts] = processFileList(fileList)

   % requiredFilesAndProducts will error if a file has syntax errors. Likely
   % in other cases too. So use try-catch and return empty if it fails.
   try
      [requiredFiles, requiredProducts] = ...
         matlab.codetools.requiredFilesAndProducts(fileList);
      requiredFiles = unique(transpose(requiredFiles));
   catch e
      requiredFiles = {};
      requiredProducts = {};
   end
   % requiredFiles = cell2table(requiredFiles,'VariableNames',{'RequiredFiles'});
end
%%
function [targetFiles, projectFiles, ignoreFiles] = prepareFileLists(...
      target, ignore, projectPath)

   % Generate a list of target filenames
   if all(isfolder(target))
      if numel(target) > 1
         error('operate one folder at a time')
      end
      targetFiles = listfiles(target, "subfolders", true, ...
         "mfiles", true, "aslist", true, "fullpath", true);
   else
      targetFiles = target;
   end

   % Generate a list of filenames included with the project
   if isempty(projectPath)
      projectFiles = projectPath;
   else
      projectFiles = listfiles(projectPath, "subfolders", true, ...
         "mfiles", true, "aslist", true, "fullpath", true);
   end

   % Generate a list of filenames to ignore
   if all(isfolder(ignore))
      if numel(ignore) > 1
         error('operate one folder at a time')
      end
      ignoreFiles = listfiles(ignore, "subfolders", true, ...
         "mfiles", true, "aslist", true, "fullpath", true);

   elseif all(isfile(ignore))
      ignoreFiles = ignore;
   else
      ignoreFiles = "";
   end

   % Remove the ignored file list from the target list
   targetFiles = setdiff(targetFiles, ignoreFiles);
end

%%
function target = validateFileList(target)
   % If a full file or folder path is provided and exists, use it directly.
   % Otherwise, if target is not a full path to an existing file or folder:
   if ~isfile(target) && ~isfolder(target)

      if isfullfile(target)
         % If a full file path was passed in but doesn't exist
         error('File does not exist or cannot be found.');

      elseif ispathlike(target)
         % If a full folder path was passed in but doesn't exist
         error('Folder does not exist or cannot be found.');

      else
         % Try to find it as a function on the MATLAB path
         target = which(target);
         if isempty(target)
            target = which(strcat(target, '.m')); % One last attempt
         end
         if isempty(target)
            error('Function does not exist on the MATLAB path.');
         end
      end
   end
end
