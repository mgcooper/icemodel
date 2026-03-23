function [requirementsList, urlList] = installRequiredFiles(requiredFiles, kwargs)
   %INSTALLREQUIREDFILES Install required files from Github.
   %
   %  INSTALLREQUIREDFILES(REQUIREDFILES)
   %  INSTALLREQUIREDFILES(PROJECTPATH=PATHNAME)
   %  INSTALLREQUIREDFILES(REQUIREMENTSFILE=FILENAME)
   %
   %  INSTALLREQUIREDFILES(_, INSTALLPATH=PATHNAME)
   %  INSTALLREQUIREDFILES(_, LOCALSOURCEPATH=PATHNAME)
   %  INSTALLREQUIREDFILES(_, IGNOREFOLDER=FOLDERNAME)
   %  INSTALLREQUIREDFILES(_, REMOTEREPONAME=REPONAME)
   %  INSTALLREQUIREDFILES(_, REMOTEBRANCH=BRANCHNAME)
   %  INSTALLREQUIREDFILES(_, GITHUBUSERNAME=USERNAME)
   %  INSTALLREQUIREDFILES(_, INSTALLMISSING=TRUE)
   %
   % Description
   %
   %  The use case for this function is to install a list of required files
   %  from GitHub. The list could be shipped with a toolbox, and third party
   %  users run an install script which reads the requirements list and installs
   %  them from GitHub. Alternatively, the toolbox maintainer can use this
   %  function to package the requirements with the toolbox.
   %
   % Input Arguments
   %
   %  The following arguments control what files get installed: Either a
   %  pre-existing list of requirements (REQUIREDFILES or REQUIREMENTSFILE), or
   %  a list of requirements generated internally by this function for a
   %  specific project folder (PROJECTPATH), optionally ignoring any
   %  requirements for files contained in IGNOREFOLDER.
   %
   %  REQUIREDFILES - (optional, positional) a list of required functions. If
   %  not provided or if empty, PROJECTPATH becomes required, and the
   %  requirements for all files in the PROJECTPATH folder are installed.
   %
   %  PROJECTPATH - (optional, name-value) a full path (scalar text) to a
   %  folder. Requirements for all files within this folder are generated and
   %  installed into the folder. The default install location is a subfolder
   %  named "dependencies" at the top level of PROJECTPATH. Specify the optional
   %  INSTALLPATH argument to control where dependencies are installed.
   %
   %  REQUIREMENTSFILE - (optional, name-value) a full path to a file containing
   %  a list of required files. NOT CURRENTLY IMPLEMENTED.
   %
   %  NOTE: If none of the three arguments above are supplied, the default
   %  behavior uses the current working directory as the PROJECTPATH parameter,
   %  and the function proceeds as though PROJECTPATH were supplied.
   %
   %  IGNOREFOLDER - A folder or array of folder names to be ignored when
   %  generating the list of requirements for the PROJECTPATH folder.
   %  IGNOREFOLDER should contain a single folder name or array of folder names
   %  which are subfolders of PROJECTPATH. Use this option to ignore a scratch/
   %  or testbed/ or sandbox/ or examples/ folder which is not under source
   %  control and is not distributed with the toolbox or project.
   %
   %  These arguments control how the requirements are found and/or resolved:
   %
   %  LOCALSOURCEPATH - folder with local versions of the required files.
   %  REMOTEREPONAME - remote (Github) repo for the localSourcePath.
   %  REMOTEBRANCH - branch to use when downloading from the remote Github repo.
   %  GITHUBUSERNAME - GitHub username for the REMOTEREPONAME.
   %
   %  These arguments control if and where files are installed:
   %
   %  INSTALLPATH - full path to location where files are installed. The default
   %  value is a folder named "dependencies" in the PROJECTPATH.
   %  DRYRUN - logical flag indicating if the files are installed. If false, the
   %  list of files are returned and printed to the screen. The default value is
   %  true.
   %
   % See also: getRequiredFiles

   arguments
      %%% The following arguments control what gets installed:
      requiredFiles (:, :) string {mustBeText} ...
         = []

      kwargs.requirementsFile (1, :) string {mustBeTextScalar} ...
         = "" % Note: requirementsFile is currently not implemented.

      kwargs.projectPath (1, :) string {mustBeFolder} ...
         = projectpath()

      kwargs.ignoreFolder (1, :) string ...
         = "testbed"

      %%% The following arguments control how requirements are found:
      kwargs.localSourcePath (1, :) {mustBeFolder} ...
         = getenv('MATLAB_FUNCTION_PATH')

      kwargs.remoteRepoName (1, :) string {mustBeTextScalar} ...
         = "matfunclib"

      kwargs.remoteBranch (1, :) string {mustBeTextScalar} ...
         = "main"

      kwargs.GitHubUserName (1, :) string {mustBeTextScalar} ...
         = getenv('GITHUB_USER_NAME')

      %%% The following arguments control where and if files get installed:

      kwargs.installPath (1, :) string {mustBeTextScalar} ...
         = fullfile(toolboxpath(), "dependencies");

      kwargs.testRun (1, 1) logical {mustBeNumericOrLogical} ...
         = false
   end

   [projectPath, ignoreFolder, localSourcePath, remoteSourcePath] = ...
      parseargs(kwargs);

   % Remember current folder and go to folder for external dependencies
   job = withcd(projectPath);

   % Find the required files
   if all(isempty(requiredFiles))
      requiredFiles = getRequiredFiles(projectPath, "ignoreList", ignoreFolder);
      requiredFiles = requiredFiles.missingFiles;
   end

   % Build a url list for the remote files
   [requirementsList, urlList] = remoteDependencyList( ...
      requiredFiles, projectPath, localSourcePath, remoteSourcePath);

   % Option to install the missing requirement locally
   fileList = kwargs.installPath + filesep + requirementsList;
   if not(kwargs.testRun)

      if ~isfolder(kwargs.installPath)
         mkdir(kwargs.installPath)
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
   else
      fprintf("\n Files will be installed to: \n %s \n", kwargs.installPath)
      fprintf("\n The following files will be installed: \n")
      disp(urlList)
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
      error('Set "remotesource" or environment variable "MATLAB_FUNCTION_PATH"')
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

      % Note - this is problematic if the requiredFiles contain files which are
      % not in localSourcePath e.g. if one project depends on another. So this
      % needs to be refactored to work with localSourcePaths (plural).

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
