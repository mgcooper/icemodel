function files = familyManifestFiles(kwargs)
   %FAMILYMANIFESTFILES List staged verification family manifest files.
   %
   %  files = icemodel.verification.helpers.familyManifestFiles()
   %  files = icemodel.verification.helpers.familyManifestFiles( ...
   %     data_root=path)
   %  files = icemodel.verification.helpers.familyManifestFiles( ...
   %     evaluation_data_root=path)
   %
   % Inputs
   %  data_root                  Whole data tree containing eval/ and input/.
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             repo-local data/eval tree is searched.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %
   % Outputs
   %  files   Sorted string column of `eval/<family>/manifest.json` paths.
   %
   % Role
   %  Operational helper used by listcases to discover staged verification
   %  families. It reads the filesystem but does not mutate setup artifacts.

   arguments
      kwargs.data_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   % Resolve exactly one paired scope through the shared precedence contract.
   % Discovery uses only its eval leaf and never composes missing families from
   % demo, test, or another active configuration.
   [root, ~] = icemodel.verification.setup.resolveStagingRoots( ...
      data_root=kwargs.data_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);

   % Find one manifest per family under the selected root. The wildcard is
   % intentionally
   % one level deep because case folders live below the family folder.
   entries = dir(fullfile(root, "*", "manifest.json"));
   if isempty(entries)
      files = strings(0, 1);
      return
   end

   files = fullfile(string({entries.folder}), string({entries.name}));
   files = sort(files(:));
end
