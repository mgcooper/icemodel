function entry = makeFirnCaseManifestEntry(values)
   %MAKEFIRNCASEMANIFESTENTRY Build one firn case manifest entry.
   %
   %  entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values)
   %
   % Inputs
   %  values   Cell array matching the canonical firn case-manifest field
   %           order (icemodel.verification.setup.firnCaseManifestFieldNames).
   %
   % Outputs
   %  entry    Struct with the canonical forcing-agnostic firn case-manifest
   %           schema.
   %
   % Role
   %  Setup helper used by the firn staging driver to make the forcing-agnostic
   %  firn case schema explicit. A driver that adds or drops a field fails early
   %  here rather than writing a shifted JSON manifest. The entry references the
   %  bundled data-only observations.mat eval target via evaluation_file and
   %  records WHICH forcing/eval sources are available (by id, INFORMATIONAL) and
   %  the colocation regime; the forcing itself is not bundled or stipulated.
   %
   % See also: icemodel.verification.setup.makeCaseManifestEntry,
   %  icemodel.verification.setup.firnCaseManifestFieldNames

   names = icemodel.verification.setup.firnCaseManifestFieldNames();

   if numel(values) ~= numel(names)
      error('firn case manifest entry expects %d values', numel(names))
   end

   entry = cell2struct(values(:), names, 1);

   % Validate the two case descriptors against the canonical vocabularies. An
   % empty surface_zone ("") / eval_target is permitted where the regime or
   % capability is not meaningful.
   icemodel.verification.setup.validateSurfaceZone(entry.surface_zone);
   icemodel.verification.setup.validateEvalTarget(entry.eval_target);
   icemodel.verification.setup.validatePermafrostZone(entry.permafrost_zone);
end
