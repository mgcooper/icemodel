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
