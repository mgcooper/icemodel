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
   %  entry    Struct with the canonical metadata-only firn case-manifest schema.
   %
   % Role
   %  Setup helper used by the firn staging driver to make the metadata-only firn
   %  case schema explicit. A driver that adds or drops a field fails early here
   %  rather than writing a shifted JSON manifest. The entry records WHICH
   %  forcing/eval sources are available (by id) and the colocation regime; it
   %  does NOT bundle evaluation.mat/reference.mat data.
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
