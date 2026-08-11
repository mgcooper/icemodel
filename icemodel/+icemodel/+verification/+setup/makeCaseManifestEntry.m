function entry = makeCaseManifestEntry(values)
   %MAKECASEMANIFESTENTRY Build one case manifest entry from canonical fields.
   %
   %  entry = icemodel.verification.setup.makeCaseManifestEntry(values)
   %
   % Inputs
   %  values   Cell array matching the canonical case manifest field order.
   %
   % Outputs
   %  entry    Struct with the canonical case manifest schema.
   %
   % Role
   %  Setup helper used by dataset importers to make manifest schemas explicit
   %  and shared across verification families.

   % Get the canonical case manifest field names.
   names = icemodel.verification.setup.caseManifestFieldNames();

   % Fail early if an importer adds or removes a field without updating the
   % shared schema. That check keeps every family manifest on one schema.
   if numel(values) ~= numel(names)
      error('case manifest entry expects %d values', numel(names))
   end

   % Build the manifest entry in the canonical field order.
   entry = cell2struct(values(:), names, 1);

   % Validate the surface_zone and eval_target descriptors against the canonical
   % vocabularies. An empty value is allowed where the regime or capability has
   % no meaning, for example an analytical Laugh-Tests benchmark.
   icemodel.verification.setup.validateSurfaceZone(entry.surface_zone);
   icemodel.verification.setup.validateEvalTarget(entry.eval_target);
   icemodel.verification.setup.validatePermafrostZone(entry.permafrost_zone);
end
