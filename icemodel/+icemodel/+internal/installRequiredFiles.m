function [requirementsList, urlList] = installRequiredFiles(requiredFiles, kwargs)
   %INSTALLREQUIREDFILES Install required files from Github.
   %
   %  INSTALLREQUIREDFILES(REQUIREDFILES)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, REQUIREMENTSFILE=FILENAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, PROJECTPATH=PATHNAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, INSTALLPATH=PATHNAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, LOCALSOURCEPATH=PATHNAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, IGNOREFOLDER=FOLDERNAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, REMOTEREPONAME=REPONAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, REMOTEBRANCH=BRANCHNAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, GITHUBUSERNAME=USERNAME)
   %  INSTALLREQUIREDFILES(REQUIREDFILES, INSTALLMISSING=TRUE)
   %
   % Description
   %  The use case for this function is to install a list of required files
   %  from GitHub. The list could be shipped with a toolbox, and third party
   %  users run an install script which reads the requirements list and installs
   %  them from GitHub. Alternatively, the toolbox maintainer can use this
   %  function to package the requirements with the toolbox.
   %
   % Input Arguments
   %  REQUIREDFILES - (optional) a list of required functions. If not provided
   %  or if empty, the requirements for the PROJECTPATH parameter or the
   %  supplied REQUIREMENTSFILE are installed.
   %
   % Optional name-value parameters
   %  projectPath - folder that needs the required files.
   %  localSourcePath - folder with local versions of the dependencies.
   %  remoteRepoName - Github repo for the localSourcePath.
   %  remoteBranch - Branch to use when downloading from the remote Github repo.
   %
   % See also: getRequiredFiles

   arguments
      requiredFiles (:, :) string {mustBeText} ...
         = []

      % Note: requirementsFile is currently not implemented.
      kwargs.requirementsFile (1, :) string {mustBeTextScalar} ...
         = ""

      kwargs.projectPath (1, :) string {mustBeFolder} ...
         = projectpath()

      kwargs.installPath (1, :) string {mustBeTextScalar} ...
         = fullfile(toolboxpath(), "dependencies");

      kwargs.localSourcePath (1, :) {mustBeFolder} ...
         = getenv('MATLAB_FUNCTION_PATH')

      kwargs.ignoreFolder (1, :) string ...
         = "testbed"

      kwargs.remoteRepoName (1, :) string {mustBeTextScalar} ...
         = "matfunclib"

      kwargs.remoteBranch (1, :) string {mustBeTextScalar} ...
         = "main"

      kwargs.GitHubUserName (1, :) string {mustBeTextScalar} ...
         = getenv('GITHUB_USER_NAME')

      kwargs.installMissing (1, 1) logical {mustBeNumericOrLogical} ...
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
   if kwargs.installMissing

      fileList = kwargs.installPath + filesep + requirementsList;
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
      requiredFiles, projectPath, localSourcePath, remoteSourcePath)
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

      if contains(requiredFilePath, localSourcePath)

         % Add file names to list of external depencies
         requirementsList(ifile) = requiredFileName;

         % Get the subfolder path relative to the top-level source repo
         relativePath = erase(requiredFilePath, localSourcePath);
         relativePath = strrep(relativePath, filesep , '/');
         if relativePath(1) == filesep
            relativePath = relativePath(2:end);
         end

         % Use '/' not fullfile b/c fullfile is platform specific
         urlList(ifile) = remoteSourcePath + '/' + relativePath + '/' ...
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
