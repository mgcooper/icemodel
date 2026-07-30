function manifest = refreshManifestSourceLists(manifest_file)
   %REFRESHMANIFESTSOURCELISTS Recompute source lists without rebuilding data.
   %
   %  manifest = icemodel.verification.setup.refreshManifestSourceLists( ...
   %     manifest_file)
   %
   % Inputs
   %  manifest_file   Path to one eval/<family>/manifest.json file.
   %
   % Outputs
   %  manifest        Updated manifest struct, also written back to
   %                  manifest_file.
   %
   % Role
   %  Source-list contract changes are manifest-only when colocation metadata
   %  and data artifacts are already current. This helper updates
   %  forcing_sources/eval_sources directly from colocation so importers do not
   %  need to rebuild observations.mat just to refresh source availability.
   %  It does not inspect or modify MAT files, observations.mat, variable
   %  schemas, embedded metadata, numeric payloads, or time axes. Use
   %  repairRcmArtifactMetadata for supported RCM-artifact migrations and a
   %  canonical importer/restage when source payload or observations changed.
   %
   % See also: icemodel.verification.setup.repairRcmArtifactMetadata

   arguments
      manifest_file (1, 1) string
   end

   manifest = jsondecode(fileread(manifest_file));
   if ~isfield(manifest, 'cases') || isempty(manifest.cases)
      icemodel.verification.setup.writeManifest(manifest_file, manifest);
      return
   end

   % Each case owns its colocation graph. Source lists are a derived view of
   % staged met/Data files and should not carry stale readiness-gating policy.
   for k = 1:numel(manifest.cases)
      if ~isfield(manifest.cases(k), 'colocation')
         continue
      end
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists( ...
         manifest.cases(k).colocation);
      manifest.cases(k).forcing_sources = cellstr(forcing_sources(:));
      manifest.cases(k).eval_sources = cellstr(eval_sources(:));
   end

   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end
