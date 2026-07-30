function [rootdir, input_path, output_path, eval_path, cleanup] = ...
      bootstrapTestEnvironment(kwargs)
   %BOOTSTRAPTESTENVIRONMENT Add test paths and install one scoped data config.
   %
   %  [rootdir, input_path, output_path, eval_path, cleanup] = ...
   %     icemodel.test.helpers.bootstrapTestEnvironment(configure_paths=true)
   %  [...] = icemodel.test.helpers.bootstrapTestEnvironment( ...
   %     icemodel_config_casename="verification", data_root="")
   %
   % Use one centralized helper for formal-suite setup. This helper adds the
   % broad icemodel/ and test/ trees to the MATLAB path, optionally installs
   % the requested config case or data-root override, and restores the caller's
   % previous MATLAB path and environment values on cleanup. Retain CLEANUP for
   % the lifetime of the scoped run. DATA_ROOT takes precedence over
   % ICEMODEL_CONFIG_CASENAME. The default remains the formal test case.

   arguments
      kwargs.configure_paths (1, 1) logical = true
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.data_root (1, 1) string = ""
   end

   % Snapshot caller state and arm cleanup before any path or environment
   % mutation so setup failures unwind to the exact entry state.
   original_path = path;
   [config_names, config_values] = snapshotConfig();
   cleanup = onCleanup(@() restoreEnvironment( ...
      original_path, config_names, config_values));

   rootdir = icemodel.internal.fullpath();

   % Add the source tree so CLI and interactive runs see the same repo layout
   addpath(genpath(fullfile(rootdir, 'icemodel')))

   % Add only test folders that actually contain MATLAB source so generated
   % artifact, profiler, and archive directories stay off the MATLAB path.
   addCodeFolders(icemodel.getpath('test'))

   % Wire the external dev-repo dependencies (exactremap, activelayer, and the
   % matfunclib helpers activelayer needs) through the single central config
   % function so they are never bootstrapped from arbitrary locations. This
   % lets the conservative (area-weighted) remap tests and the Obu
   % permafrost-zone analysis RUN instead of self-skipping. Each dependency is
   % a clean no-op when already on the path or absent (the affected tests then
   % skip cleanly, as before).
   icemodel.dependencies();

   if kwargs.configure_paths
      cfg = resolveConfig(kwargs.icemodel_config_casename, ...
         kwargs.data_root, true);
   else
      cfg = resolveConfig(kwargs.icemodel_config_casename, ...
         kwargs.data_root, false);
   end

   % Add the selected data tree itself when it exists. Direct unit runs and
   % fixture loaders may resolve bundled files relative to this root, so keeping
   % the selected parent on path avoids per-file bootstrap duplication.
   addDataFolder(string(cfg.ICEMODEL_DATA_PATH))

   % Return the resolved canonical test paths for callers that need them.
   input_path = string(cfg.ICEMODEL_INPUT_PATH);
   output_path = string(cfg.ICEMODEL_OUTPUT_PATH);
   eval_path = string(cfg.ICEMODEL_EVAL_PATH);
end

function [names, values] = snapshotConfig()
   %SNAPSHOTCONFIG Capture every environment value written by icemodel.config.

   % Derive the complete schema from config itself, then read the caller's raw
   % environment instead of the resolved preview. Empty raw values remain empty
   % so cleanup unsets variables that were genuinely unset on entry.
   names = string(fieldnames(icemodel.config('setenv', false)));
   values = cell(numel(names), 1);
   for n = 1:numel(names)
      values{n} = getenv(names(n));
   end
end

function cfg = resolveConfig(icemodel_config_casename, data_root, do_setenv)
   %RESOLVECONFIG Resolve data-root precedence with optional installation.

   % A caller-supplied root owns every derived path and intentionally makes the
   % case selector irrelevant for this scoped run.
   if ~isblanktext(data_root)
      cfg = icemodel.config('ICEMODEL_DATA_PATH', data_root, ...
         'setenv', do_setenv);
      return
   end

   % Otherwise resolve the documented case through the central config API.
   cfg = icemodel.config('casename', icemodel_config_casename, ...
      'setenv', do_setenv);
end

function restoreConfig(names, values)
   %RESTORECONFIG Restore every caller value changed by icemodel.config.

   for n = 1:numel(names)
      setenv(names(n), values{n});
   end
end

function restoreEnvironment(original_path, names, values)
   %RESTOREENVIRONMENT Restore configuration and the exact caller MATLAB path.

   % Restore both scoped resources from one cleanup callback so normal return
   % and setup or suite errors follow the same unwind path.
   restoreConfig(names, values)
   path(original_path)
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

function addDataFolder(rootdir)
   %ADDDATAFOLDER Add the selected scoped data root when present.

   if isfolder(rootdir)
      addpath(char(rootdir))
   end
end
