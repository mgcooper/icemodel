function [data_root, met_root] = selectedDataRoot(met_dir)
   %SELECTEDDATAROOT Resolve one selected met path to its data and met roots.
   %
   %  [data_root, met_root] = ...
   %     icemodel.forcing.reconstruct.selectedDataRoot(met_dir)
   %
   % Supports the standard <data>/input/met[/<source>] layout, the compact
   % <data>/met[/<source>] layout, and caller-owned custom source
   % directories.
   %
   % See also: icemodel.forcing.reconstruct.acceptanceWindow,
   %  icemodel.forcing.reconstruct.verifyPromiceFilledReadiness,
   %  icemodel.verification.report.buildGapFillReport

   arguments
      met_dir (1, 1) string
   end

   % A path ending in met is the flat met root; otherwise its parent is the
   % root that can contain per-source folders and flat fallback artifacts.
   [parent, name] = fileparts(met_dir);
   if string(name) == "met"
      met_root = char(met_dir);
   else
      met_root = parent;
   end

   % Strip the standard input/met suffix or the compact met suffix while
   % preserving custom roots.
   [met_parent, met_name] = fileparts(met_root);
   if string(met_name) ~= "met"
      data_root = met_root;
      return
   end
   [canonical_root, parent_name] = fileparts(met_parent);
   if string(parent_name) == "input"
      data_root = canonical_root;
   else
      data_root = met_parent;
   end
end
