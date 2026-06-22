function [rootdir, input_path, output_path, eval_path, cleanup] = ...
      bootstrapTestEnvironment(kwargs)
   %BOOTSTRAPTESTENVIRONMENT Add test paths and install the canonical test config.
   %
   %  rootdir = icemodel.test.helpers.bootstrapTestEnvironment()
   %  [rootdir, input_path, output_path, eval_path, cleanup] = ...
   %     icemodel.test.helpers.bootstrapTestEnvironment(configure_paths=true)
   %
   % Use one centralized helper for formal-suite setup. This helper adds the
   % broad icemodel/ and test/ trees to the MATLAB path, optionally installs
   % the canonical test config via `icemodel.config('casename','test')`, and
   % restores the caller's previous environment values on cleanup.

   arguments
      kwargs.configure_paths (1, 1) logical = true
   end

   rootdir = icemodel.internal.fullpath();

   % Add the source tree so CLI and interactive runs see the same repo layout
   addpath(genpath(fullfile(rootdir, 'icemodel')))

   % Add only test folders that actually contain MATLAB source so generated
   % artifact, profiler, and archive directories stay off the MATLAB path.
   addCodeFolders(icemodel.getpath('test'))

   % Put the exactremap toolbox on the path so the conservative
   % (area-weighted) remap tests RUN instead of self-skipping. Resolved from
   % the sibling exactremap dev repo; a clean no-op when it is absent (the
   % affected tests then skip cleanly, as before).
   ensureExactremap(rootdir)

   % Put the activelayer toolbox (and the matfunclib helpers it depends on) on
   % the path so icemodel.verification reads the Obu permafrost zones through
   % activelayer.readobuzones instead of duplicating shaperead logic. Resolved
   % from the sibling dev repos; a clean no-op when they are absent.
   ensureActivelayer(rootdir)

   % Return a no-op cleanup by default so callers can always keep one handle
   % alive even when path/config installation is disabled.
   cleanup = onCleanup(@() []);
   if kwargs.configure_paths
      [cfg, cleanup] = installTestConfig();
   else
      cfg = icemodel.config('casename', 'test', 'setenv', false);
   end

   % Return the resolved canonical test paths for callers that need them.
   input_path = string(cfg.ICEMODEL_INPUT_PATH);
   output_path = string(cfg.ICEMODEL_OUTPUT_PATH);
   eval_path = string(cfg.ICEMODEL_EVAL_PATH);
end

function [cfg, cleanup] = installTestConfig()
   %INSTALLTESTCONFIG Install the canonical formal-suite config with restoration.

   % Snapshot the caller's current config before installing the canonical
   % formal-suite config.
   cfg = icemodel.config('getenv', true);
   [names, values] = deal(string(fieldnames(cfg)), struct2cell(cfg));

   % Install the canonical test config as the single source of truth for
   % the formal suite, then restore the caller's prior config on cleanup.
   cfg = icemodel.config('casename', 'test');

   cleanup = onCleanup(@() restoreConfig(names, values));
end

function restoreConfig(names, values)
   %RESTORECONFIG Restore the caller's original ICEMODEL_* env values.

   for n = 1:numel(names)
      setenv(names(n), values{n});
   end
end

function ensureExactremap(rootdir)
   %ENSUREEXACTREMAP Put the exactremap toolbox on the path if not already.
   % exactremap supplies the conservative (overlap-area-weighted) polygon
   % remap used by icemodel.forcing.helpers.remapPolygon. It lives in a
   % sibling dev repo (projects/exactremap, https://github.com/mgcooper/
   % exactremap), kept off the icemodel repo. When it is already on the path
   % (e.g. activated by the user's startup.m) this is a no-op; when the dev
   % repo is missing this returns quietly and the remap tests skip cleanly.

   if ~isempty(which('exactremap'))
      return
   end

   % projects/icemodel/.. -> projects/, then projects/exactremap/toolbox.
   projects_root = fileparts(rootdir);
   tbdir = fullfile(projects_root, 'exactremap', 'toolbox');
   if isfolder(tbdir)
      addpath(genpath(tbdir))
   end
end

function ensureActivelayer(rootdir)
   %ENSUREACTIVELAYER Put activelayer + matfunclib on the path if not already.
   % activelayer.readobuzones reads the Obu (UiO PEX) permafrost-zones
   % shapefile; it depends on shared helpers from matfunclib
   % (parseFileName, dealout, mustBePolygonOrFile). Both live in sibling dev
   % repos (projects/activelayer, projects/matfunclib), kept off the icemodel
   % repo. When activelayer.readobuzones already resolves (e.g. activated by
   % the user's startup.m) this is a no-op; when a dev repo is missing this
   % returns quietly and the Obu-dependent analysis tools skip cleanly.

   if ~isempty(which('activelayer.readobuzones'))
      return
   end

   % projects/icemodel/.. -> projects/, then the sibling toolbox trees.
   % matfunclib goes on first so activelayer's helper dependencies resolve.
   projects_root = fileparts(rootdir);
   mfdir = fullfile(projects_root, 'matfunclib');
   if isfolder(mfdir)
      addpath(genpath(mfdir))
   end
   tbdir = fullfile(projects_root, 'activelayer', 'toolbox');
   if isfolder(tbdir)
      addpath(genpath(tbdir))
   end
end

function addCodeFolders(rootdir)
   %ADDCODEFOLDERS Add only source-bearing folders under the test tree.

   % Discover the folders that contain MATLAB source files under the
   % requested root, then add each folder once in stable order.
   files = dir(fullfile(rootdir, '**', '*.m'));
   folders = unique(string({files.folder}), 'stable');
   for n = 1:numel(folders)
      addpath(char(folders(n)))
   end
end
