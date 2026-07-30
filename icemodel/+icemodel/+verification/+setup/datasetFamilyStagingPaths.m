function [family_root, manifest_file, met_outdir, userdata_outdir] = ...
      datasetFamilyStagingPaths(evaluation_data_root, input_root, dataset_family)
   %DATASETFAMILYSTAGINGPATHS Build the shared dataset-family output paths.
   %
   %  [family_root, manifest_file, met_outdir, userdata_outdir] = ...
   %     icemodel.verification.setup.datasetFamilyStagingPaths( ...
   %     evaluation_data_root, input_root, dataset_family)
   %
   % Importers resolve the two configurable roots separately, then use this
   % helper to name the family manifest and the standard runtime input folders.

   arguments
      evaluation_data_root (1, 1) string
      input_root (1, 1) string
      dataset_family (1, 1) string
   end

   % Keep the four shared path relationships in one place so family adapters
   % cannot drift between eval, manifest, met, and userdata roots.
   family_root = fullfile(evaluation_data_root, dataset_family);
   manifest_file = fullfile(family_root, "manifest.json");
   met_outdir = fullfile(input_root, "met");
   userdata_outdir = fullfile(input_root, "userdata");
end
