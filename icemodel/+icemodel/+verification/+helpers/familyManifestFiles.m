function files = familyManifestFiles(kwargs)
   %FAMILYMANIFESTFILES List staged verification family manifest files.
   %
   %  files = icemodel.verification.helpers.familyManifestFiles()
   %  files = icemodel.verification.helpers.familyManifestFiles( ...
   %     evaluation_data_root=path)
   %
   % Inputs
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             path is resolved from icemodel.config.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %
   % Outputs
   %  files   Sorted string column of `eval/snow/<family>/manifest.json` paths.
   %
   % Role
   %  Operational helper used by listcases to discover staged verification
   %  families. It reads the filesystem but does not mutate setup artifacts.

   arguments
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
   end

   % Resolve only the snow-verification root; dataset-family subfolders are
   % discovered dynamically so new families do not need code changes here.
   snow_data_root = icemodel.verification.helpers.snowDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Find one manifest per family. The wildcard is intentionally one level deep
   % because case folders live below the family folder.
   entries = dir(fullfile(snow_data_root, "*", "manifest.json"));
   if isempty(entries)
      files = strings(0, 1);
      return
   end

   files = fullfile(string({entries.folder}), string({entries.name}))';
   files = sort(files); %#ok<TRSRT>
end
