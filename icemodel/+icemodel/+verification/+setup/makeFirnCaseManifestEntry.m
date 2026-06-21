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
   %  entry    Struct with the canonical firn case-manifest schema.
   %
   % Role
   %  Setup helper used by the firn staging driver to make the firn case
   %  schema explicit. It mirrors makeCaseManifestEntry but binds the
   %  firn-only field set (site_location, colocated_forcing) so a driver
   %  that adds or drops a field fails early rather than writing a shifted
   %  JSON manifest.
   %
   % See also: icemodel.verification.setup.makeCaseManifestEntry,
   %  icemodel.verification.setup.firnCaseManifestFieldNames

   names = icemodel.verification.setup.firnCaseManifestFieldNames();

   if numel(values) ~= numel(names)
      error('firn case manifest entry expects %d values', numel(names))
   end

   entry = cell2struct(values(:), names, 1);

   % Validate surface_zone against the canonical vocabulary. An empty value
   % ("") is permitted for cases where the surface regime is not meaningful.
   icemodel.verification.setup.validateSurfaceZone(entry.surface_zone);
end
