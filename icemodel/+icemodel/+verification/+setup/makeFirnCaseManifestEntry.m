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
   %  Setup helper that the firn staging driver uses to state the
   %  forcing-agnostic firn case schema. A driver that adds or drops a field
   %  fails here instead of writing a shifted JSON manifest. The entry points
   %  at the bundled data-only observations.mat eval target through
   %  evaluation_file. It records WHICH forcing and eval sources are available
   %  (by id, INFORMATIONAL) and the colocation regime. It does not bundle or
   %  require the forcing itself.
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
