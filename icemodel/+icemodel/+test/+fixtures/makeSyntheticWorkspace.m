function workspace = makeSyntheticWorkspace(simyears, kwargs)
   %MAKESYNTHETICWORKSPACE Create an isolated icemodel test workspace.
   %
   %  workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyears)
   %  workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyears, ...
   %     configure=true, nsteps=24, dt_seconds=3600)
   %  workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyears, ...
   %     parentdir=fullfile(icemodel.getpath('test'), 'artifacts', 'tmp'))
   %
   % The returned struct records the created paths and the previous icemodel
   % environment-variable values so tests can restore them during teardown.

   arguments
      simyears (:, 1) double {mustBeInteger, mustBePositive}
      kwargs.configure (1, 1) logical = true
      kwargs.sitename (1, :) char = 'kanm'
      kwargs.forcings (1, :) char = 'kanm'
      kwargs.nsteps (1, 1) double {mustBeInteger, mustBePositive} = 24
      kwargs.dt_seconds (1, 1) double {mustBePositive} = 3600
      kwargs.include_modis (1, 1) logical = false
      kwargs.parentdir (1, :) char = ''
      kwargs.rootdir (1, :) char = ''
   end

   % Allow callers to pin the workspace under a known parent folder or to
   % provide a fully resolved root path when a test needs stable artifact
   % placement for inspection/debugging.
   rootdir = resolveRootdir(kwargs.rootdir, kwargs.parentdir);
   inputdir = fullfile(rootdir, 'input');
   metdir = fullfile(inputdir, 'met');
   spectraldir = fullfile(inputdir, 'spectral');
   userdatadir = fullfile(inputdir, 'userdata');
   evaldir = fullfile(rootdir, 'eval');
   outputdir = fullfile(rootdir, 'output');

   % Create the minimal directory tree expected by the model/tooling code.
   mkdir(rootdir);
   mkdir(inputdir);
   mkdir(metdir);
   mkdir(userdatadir);
   mkdir(evaldir);
   mkdir(outputdir);

   % Copy the canonical spectral inputs from the repo-local demo dataset.
   spectral_src = fullfile(icemodel.getpath('demo'), 'data', 'input', ...
      'spectral');
   copyfile(spectral_src, spectraldir);

   % Write one synthetic met file per requested simulation year.
   for thisyear = reshape(simyears, 1, [])
      icemodel.test.fixtures.writeSyntheticMetFile(metdir, thisyear, ...
         'sitename', kwargs.sitename, ...
         'forcings', kwargs.forcings, ...
         'nsteps', kwargs.nsteps, ...
         'dt_seconds', kwargs.dt_seconds, ...
         'include_modis', kwargs.include_modis);
   end

   % Snapshot the caller's current icemodel config so teardown can restore
   % whatever environment was in place before the fixture ran.
   oldenv = icemodel.config('getenv', true);

   % Optionally point the process at this isolated workspace immediately.
   if kwargs.configure
      icemodel.config('ICEMODEL_DATA_PATH', rootdir, ...
         'ICEMODEL_INPUT_PATH', inputdir, ...
         'ICEMODEL_OUTPUT_PATH', outputdir, ...
         'ICEMODEL_EVAL_PATH', evaldir, ...
         'ICEMODEL_USERDATA_PATH', userdatadir);
   end

   % Return the created paths and fixture metadata for later reuse/cleanup.
   workspace = struct();
   workspace.rootdir = rootdir;
   workspace.inputdir = inputdir;
   workspace.metdir = metdir;
   workspace.spectraldir = spectraldir;
   workspace.userdatadir = userdatadir;
   workspace.evaldir = evaldir;
   workspace.outputdir = outputdir;
   workspace.simyears = reshape(simyears, 1, []);
   workspace.sitename = kwargs.sitename;
   workspace.forcings = kwargs.forcings;
   workspace.nsteps = kwargs.nsteps;
   workspace.dt_seconds = kwargs.dt_seconds;
   workspace.oldenv = oldenv;
end

function rootdir = resolveRootdir(rootdir_in, parentdir)
   %RESOLVEROOTDIR Choose the root folder for a synthetic workspace.

   if ~isempty(rootdir_in)
      rootdir = rootdir_in;
      return
   end

   if ~isempty(parentdir)
      if ~isfolder(parentdir)
         mkdir(parentdir);
      end
      rootdir = tempname(parentdir);
   else
      rootdir = tempname;
   end
end
