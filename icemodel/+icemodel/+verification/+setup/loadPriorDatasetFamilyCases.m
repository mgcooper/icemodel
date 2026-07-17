function cases = loadPriorDatasetFamilyCases(manifest_file, kwargs)
   %LOADPRIORDATASETFAMILYCASES Read cases needed by an additive native refresh.
   %
   %  cases = icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
   %     manifest_file, build_forcing=false, overwrite_family=false)
   %
   % Inputs
   %  manifest_file     Existing family manifest path.
   %  build_forcing     True when native artifacts are rebuilt in this call.
   %  overwrite_family  True at the explicit whole-family replacement boundary.
   %
   % Outputs
   %  cases   Existing manifest cases only when an additive native observation
   %          refresh must preserve prior artifact legs; otherwise empty.

   arguments
      manifest_file (1, 1) string
      kwargs.build_forcing (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
   end

   % A forcing rebuild or whole-family replacement owns the fresh native state.
   % Missing manifests represent new family roots and therefore have no prior
   % cases to preserve.
   cases = struct([]);
   if kwargs.build_forcing || kwargs.overwrite_family || ~isfile(manifest_file)
      return
   end

   % Use the canonical manifest reader so text and case fields have the same
   % normalized types as every other verification workflow.
   manifest = icemodel.verification.helpers.readFamilyManifest(manifest_file);
   cases = manifest.cases;
end
