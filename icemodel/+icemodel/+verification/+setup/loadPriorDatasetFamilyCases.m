function cases = loadPriorDatasetFamilyCases(manifest_file, kwargs)
   %LOADPRIORDATASETFAMILYCASES Read cases needed by an additive native refresh.
   %
   %  cases = icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
   %     manifest_file, overwrite_family=false, build_observations=true)
   %
   % Inputs
   %  manifest_file     Existing family manifest path.
   %  overwrite_family  True at the explicit whole-family replacement boundary.
   %  build_observations True when this call rebuilds observation contracts.
   %
   % Outputs
   %  cases   Existing cases for an additive refresh or a forcing-only
   %          replacement. A whole-family replacement that also builds
   %          observations returns empty.

   arguments
      manifest_file (1, 1) string
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.build_observations (1, 1) logical = true
   end

   % A whole-family replacement that also builds observations regenerates every
   % source, so it needs no prior cases. A forcing-only replacement still needs
   % the existing observation contract.
   cases = struct([]);
   if (kwargs.overwrite_family && kwargs.build_observations) ...
         || ~isfile(manifest_file)
      return
   end

   % Use the canonical manifest reader so text and case fields have the same
   % normalized types as every other verification workflow.
   manifest = icemodel.verification.helpers.readFamilyManifest(manifest_file);
   cases = manifest.cases;
end
