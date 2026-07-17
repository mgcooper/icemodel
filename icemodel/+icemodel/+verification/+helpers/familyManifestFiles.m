function files = familyManifestFiles(kwargs)
   %FAMILYMANIFESTFILES List staged verification family manifest files.
   %
   %  files = icemodel.verification.helpers.familyManifestFiles()
   %  files = icemodel.verification.helpers.familyManifestFiles( ...
   %     evaluation_data_root=path)
   %
   % Inputs
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             repo-local data/eval tree is searched first,
   %                             then committed demo/test fallback families.
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
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   % Resolve the staged root first. Default discovery also includes the committed
   % demo/test root as a per-family fallback so clean checkouts can still list
   % families that have not yet been staged under top-level data/eval.
   roots = manifestRoots(kwargs.evaluation_data_root, ...
      kwargs.icemodel_config_casename);

   % Find one manifest per family under each root. The wildcard is intentionally
   % one level deep because case folders live below the family folder.
   entry_groups = cell(numel(roots), 1);
   n_groups = 0;
   for root = reshape(roots, 1, [])
      found = dir(fullfile(root, "*", "manifest.json"));
      if isempty(found)
         continue
      end
      n_groups = n_groups + 1;
      entry_groups{n_groups} = found(:);
   end
   if n_groups == 0
      files = strings(0, 1);
      return
   end

   entries = vertcat(entry_groups{1:n_groups});
   if isempty(entries)
      files = strings(0, 1);
      return
   end

   files = fullfile(string({entries.folder}), string({entries.name}))';
   files = deduplicateFamilyFiles(files);
   files = sort(files);
end

function roots = manifestRoots(evaluation_data_root, icemodel_config_casename)
   %MANIFESTROOTS Resolve default and explicit manifest discovery roots.
   primary = icemodel.verification.helpers.evaluationDataRoot( ...
      evaluation_data_root=evaluation_data_root, ...
      icemodel_config_casename=icemodel_config_casename);
   roots = primary;
   if isblanktext(evaluation_data_root) && isblanktext(icemodel_config_casename)
      fallback = [
         string(fullfile(icemodel.internal.fullpath('data', 'test', 'eval')))
         string(fullfile(icemodel.internal.fullpath('demo', 'data', 'eval')))];
      roots = unique([primary; fallback], 'stable');
   end
end

function files = deduplicateFamilyFiles(files)
   %DEDUPLICATEFAMILYFILES Keep the first manifest for each family folder.
   families = strings(numel(files), 1);
   for k = 1:numel(files)
      manifest_folder = fileparts(files(k));
      [~, families(k)] = fileparts(manifest_folder);
   end
   [~, keep] = unique(families, 'stable');
   files = files(sort(keep));
end
