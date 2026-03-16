function workspace = makeSyntheticWorkspace(simyears, kwargs)
%MAKESYNTHETICWORKSPACE Create an isolated icemodel test workspace.
%
%  workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyears)
%  workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyears, ...
%     configure=true, nsteps=24, dt_seconds=3600)
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
   end

   rootdir = tempname;
   inputdir = fullfile(rootdir, 'input');
   metdir = fullfile(inputdir, 'met');
   spectraldir = fullfile(inputdir, 'spectral');
   userdatadir = fullfile(inputdir, 'userdata');
   evaldir = fullfile(rootdir, 'eval');
   outputdir = fullfile(rootdir, 'output');

   mkdir(rootdir);
   mkdir(inputdir);
   mkdir(metdir);
   mkdir(userdatadir);
   mkdir(evaldir);
   mkdir(outputdir);

   spectral_src = icemodel.internal.fullpath( ...
      'demo', 'data', 'input', 'spectral');
   copyfile(spectral_src, spectraldir);

   for thisyear = reshape(simyears, 1, [])
      icemodel.test.fixtures.writeSyntheticMetFile(metdir, thisyear, ...
         'sitename', kwargs.sitename, ...
         'forcings', kwargs.forcings, ...
         'nsteps', kwargs.nsteps, ...
         'dt_seconds', kwargs.dt_seconds, ...
         'include_modis', kwargs.include_modis);
   end

   envnames = {'ICEMODEL_DATA_PATH', 'ICEMODEL_INPUT_PATH', ...
      'ICEMODEL_OUTPUT_PATH', 'ICEMODEL_EVAL_PATH', 'ICEMODEL_USERDATA_PATH'};
   oldenv = struct();
   for n = 1:numel(envnames)
      oldenv.(envnames{n}) = getenv(envnames{n});
   end

   if kwargs.configure
      icemodel.config('ICEMODEL_DATA_PATH', rootdir, ...
         'ICEMODEL_INPUT_PATH', inputdir, ...
         'ICEMODEL_OUTPUT_PATH', outputdir, ...
         'ICEMODEL_EVAL_PATH', evaldir, ...
         'ICEMODEL_USERDATA_PATH', userdatadir);
   end

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
